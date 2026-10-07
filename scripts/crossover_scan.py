#!/usr/bin/env python3
"""Crossover calls for every Athilafinder intact element (k-mers only).

Per element (internal region between the LTRs):
  * family-specific k-mer markers -> independent sites (islands)
  * changepoint likelihood ratio (llr_switch) for one A->B switch; a crossover is
    called above --llr-threshold (calibrated on held-out clean elements: 0.5% FPR)
  * nesting / merged-element flags:
      ltr_inside   >= --ltr-sites islands of LTR-specific k-mers (from Athilafinder
                   solo LTRs, absent from internal regions) inside the internal region
      duplicated / mixed_strand / out_of_order   from the element's TEsorter domains
  * breakpoint placed relative to the element's own domains

Usage: crossover_scan.py --markers M --ltr-markers ltr_only.txt --internal DIR \
          --summary-dir DIR --tesorter DIR --out-dir OUT [--llr-threshold 26.1]
"""
import argparse, glob, os, re
from collections import Counter, defaultdict
import lib_ids as L
from athila_recomb_nt import read_fasta, revcomp
from island_scan_lib import islands
from tract_scan import load_markers
from llr_lib import llr_switch

ap = argparse.ArgumentParser()
for a in ("--markers", "--ltr-markers", "--internal", "--summary-dir", "--tesorter", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--llr-threshold", type=float, default=26.1)
ap.add_argument("--ltr-sites", type=int, default=5)
ap.add_argument("--k", type=int, default=31)
args, _ = ap.parse_known_args()
os.makedirs(args.out_dir, exist_ok=True)
DOM = ["GAG", "PROT", "RT", "RH", "INT"]

marker = load_markers(args.markers)
ltr = set()
for line in open(args.ltr_markers):
    km = line.split()[0].upper()
    ltr.add(km); ltr.add(revcomp(km))
l5 = {}
for p in glob.glob(os.path.join(args.summary_dir, "*SUMMARY_TABLE.txt")):
    acc = L.accession_from_path(p)
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    h = rows[0]
    for r in rows[1:]:
        d = dict(zip(h, r))
        if d.get("quality") == "intact":
            try:
                l5[L.nsid(acc, d["TE_ID"])] = int(d["LTR5_length"])
            except ValueError:
                pass
doms = defaultdict(list)
for p in glob.glob(os.path.join(args.tesorter, "*.dom.gff3")):
    acc = L.accession_from_path(p)
    for line in open(p):
        f = line.split("\t")
        if len(f) < 9:
            continue
        g = re.search(r"gene=([A-Za-z]+)", f[8])
        if g and g.group(1) in DOM:
            doms[L.nsid(acc, f[0])].append((g.group(1), int(f[3]), int(f[4]), f[6]))

out_rows, primary = [], {}
for p in sorted(glob.glob(os.path.join(args.internal, "*.fasta"))):
    acc = L.accession_from_path(p)
    for el, s in read_fasta(p).items():
        ns = L.nsid(acc, el); s = s.upper()
        hits = [(i, marker[s[i:i + args.k]]) for i in range(len(s) - args.k + 1) if s[i:i + args.k] in marker]
        isl = islands(hits)
        if isl:
            primary[ns] = Counter(f for f, a, b in isl).most_common(1)[0][0]
        llr, fa, fb, bp = llr_switch(s, marker, args.k)
        if llr < args.llr_threshold or bp is None:
            continue
        flags = []
        lh = [(i, "L") for i in range(len(s) - args.k + 1) if s[i:i + args.k] in ltr]
        n_ltr = len(islands(lh))
        if n_ltr >= args.ltr_sites:
            flags.append("ltr_inside")
        d = sorted(doms.get(ns, []), key=lambda x: x[1])
        cnt = Counter(x[0] for x in d)
        if any(v > 1 for v in cnt.values()):
            flags.append("duplicated")
        st = {x[3] for x in d}
        if len(st) > 1:
            flags.append("mixed_strand")
        elif cnt and max(cnt.values()) == 1:
            order = [x[0] for x in (d if "+" in st else d[::-1])]
            if order != [kk for kk in DOM if kk in order]:
                flags.append("out_of_order")
        pos = bp + l5.get(ns, 0)
        where = "outside annotated domains"
        for x in d:
            if x[1] <= pos <= x[2]:
                where = f"within {x[0]}"
        for x, y in zip(d, d[1:]):
            if x[2] < pos < y[1]:
                where = f"between {x[0]} and {y[0]}"
        cls = "nested_or_merged" if flags else "crossover"
        co = L.element_coords(el) or {}
        out_rows.append([ns, acc, co.get("chrom", ""), co.get("start", ""), fa, fb, f"{llr:.1f}",
                         cls, ";".join(flags), n_ltr, where, ",".join(x[0] for x in d)])

with open(os.path.join(args.out_dir, "crossover_calls.tsv"), "w") as f:
    f.write("ns_id\taccession\tchrom\tstart\tfamily_5p\tfamily_3p\tllr\tclass\tnesting_flags\t"
            "ltr_sites_inside\tbreakpoint\tdomains\n")
    for r in out_rows:
        f.write("\t".join(map(str, r)) + "\n")
with open(os.path.join(args.out_dir, "primary_family.tsv"), "w") as f:
    f.write("ns_id\tprimary_family\n")
    for ns, fam in primary.items():
        f.write(f"{ns}\tATHILA{fam}\n")
c = Counter(r[7] for r in out_rows)
print(f"elements scanned: {len(primary)}; switch calls: {len(out_rows)} -> " + ", ".join(f"{k} {v}" for k, v in c.items()))
print("nesting flags:", dict(Counter(x for r in out_rows for x in r[8].split(";") if x)))
pairs = defaultdict(Counter)
for r in out_rows:
    pairs["-".join(sorted((r[4], r[5])))][r[7]] += 1
print("\nby pair (crossover / nested):")
for p, cc in sorted(pairs.items(), key=lambda x: -x[1]["crossover"]):
    print(f"  {p:>7}: crossover {cc['crossover']:>4}  nested {cc['nested_or_merged']:>3}")
print("\nbreakpoints of crossovers:", Counter(r[10] for r in out_rows if r[7] == "crossover").most_common())
