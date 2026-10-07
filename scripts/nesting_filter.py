#!/usr/bin/env python3
"""Crossover vs nesting for crossover-like (A>B) mosaic elements.

Uses each element's own TEsorter domains (dom.gff3, element coordinates), the
family-specific k-mer markers (to place the switch), and the split alignment to the
TAIR12 LTR exemplars (to find LTR sequence inside the element).

Nesting / merged-element evidence:
  internal_ltr     an LTR-derived piece > --end-margin bp from both element ends
  duplicated       a core domain (GAG..INT) occurs more than once
  mixed_strand     core domains on both strands
  out_of_order     core domains not in GAG<PROT<RT<RH<INT order along the strand
  oversized        element > --oversize x the median length of its host family
Crossover: none of the above, >= --min-domains single-copy core domains present,
  and the switch (last host marker -> first donor marker) lies inside the element's
  coding span.

Usage:
    nesting_filter.py --mosaics crossovers.tsv --paf junctions_full/aln.paf \
        --tesorter DIR --markers robust_markers.tsv --fulllength DIR \
        --elements inventory/elements.tsv --out-dir OUT
"""
import argparse
import glob
import os
import re
import statistics
from collections import Counter, defaultdict

import lib_ids as L
from athila_recomb_nt import read_fasta
from junction_scan import read_paf
from tract_scan import load_markers

ap = argparse.ArgumentParser()
for a in ("--mosaics", "--paf", "--tesorter", "--markers", "--fulllength", "--elements", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--end-margin", type=int, default=300)
ap.add_argument("--oversize", type=float, default=1.5)
ap.add_argument("--min-domains", type=int, default=3)
ap.add_argument("--k", type=int, default=31)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
DOM = ["GAG", "PROT", "RT", "RH", "INT"]

mos = {}
for l in open(args.mosaics).readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    mos[f[0]] = [s.split("(")[0] for s in f[6].split(">")]

# element domains from gff3 (element coordinates)
doms = defaultdict(list)
for p in glob.glob(os.path.join(args.tesorter, "*.dom.gff3")):
    acc = L.accession_from_path(p)
    for line in open(p):
        f = line.split("\t")
        if len(f) < 9:
            continue
        ns = L.nsid(acc, f[0])
        if ns not in mos:
            continue
        g = re.search(r"gene=([A-Za-z]+)", f[8])
        if g and g.group(1) in DOM:
            doms[ns].append((g.group(1), int(f[3]), int(f[4]), f[6]))

# family length medians (Athilafinder element length) by k-mer primary family
length = {}
for l in open(args.elements).readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    if f[4] and f[5]:
        length[f[0]] = int(f[5]) - int(f[4])
fam_len = defaultdict(list)
prim = "/".join(args.mosaics.split("/")[:-1]) + "/primary_family.tsv"
for l in open(prim).readlines()[1:]:
    ns, fam = l.rstrip("\n").split("\t")
    if ns in length:
        fam_len[fam.replace("ATHILA", "")].append(length[ns])
fam_med = {f: statistics.median(v) for f, v in fam_len.items() if len(v) >= 20}

by_q = read_paf(args.paf, 80)
marker = load_markers(args.markers)
seqs = {}
for flf in glob.glob(os.path.join(args.fulllength, "*DP_FULLLENGTH_RENAMED.fasta")):
    acc = L.accession_from_path(flf)
    pre = acc + "::"
    if any(m.startswith(pre) for m in mos):
        for el, s in read_fasta(flf).items():
            ns = L.nsid(acc, el)
            if ns in mos:
                seqs[ns] = s.upper()

rows = []
for ns, path in mos.items():
    host, donor = path[0], path[-1]
    d = sorted(doms.get(ns, []), key=lambda x: x[1])
    s = seqs.get(ns, "")
    reasons, notes = [], []
    # internal LTR
    pcs = by_q.get(ns, [])
    for p in pcs:
        if p["part"] != "LTR":
            continue
        left = any(q["part"] == "I" and q["qe"] <= p["qs"] + 50 and q["qe"] - q["qs"] >= 300 for q in pcs)
        right = any(q["part"] == "I" and q["qs"] >= p["qe"] - 50 and q["qe"] - q["qs"] >= 300 for q in pcs)
        if left and right:            # LTR with internal sequence on both sides = nested
            reasons.append("internal_ltr")
            notes.append(f"{p['fam'][6:]} LTR inside at {p['qs']}-{p['qe']}")
            break
    cnt = Counter(x[0] for x in d)
    if any(v > 1 for v in cnt.values()):
        reasons.append("duplicated")
        notes.append("repeated " + ",".join(k for k, v in cnt.items() if v > 1))
    strands = {x[3] for x in d}
    if len(strands) > 1:
        reasons.append("mixed_strand")
    if len(strands) == 1 and cnt and max(cnt.values()) == 1:
        order = [x[0] for x in (d if "+" in strands else d[::-1])]
        if order != [k for k in DOM if k in order]:
            reasons.append("out_of_order")
            notes.append("order " + "-".join(order))
    hf = host
    if hf in fam_med and len(s) > args.oversize * fam_med[hf]:
        notes.append(f"long: {len(s)} bp vs family median {fam_med[hf]:.0f}")
    # switch position from markers
    sw = ""
    if s:
        hits = [(i, marker[s[i:i + args.k]]) for i in range(len(s) - args.k + 1) if s[i:i + args.k] in marker]
        hp = [i for i, f in hits if f == host]
        dp = [i for i, f in hits if f == donor]
        if hp and dp:
            if statistics.median(hp) < statistics.median(dp):
                a, b = max(i for i in hp if i < statistics.median(dp)), min(i for i in dp if i > statistics.median(hp))
            else:
                a, b = max(i for i in dp if i < statistics.median(hp)), min(i for i in hp if i > statistics.median(dp))
            mid = (a + b) // 2
            where = "outside coding span"
            for x in d:
                if x[1] <= mid <= x[2]:
                    where = f"within {x[0]}"
            for x, y in zip(d, d[1:]):
                if x[2] < mid < y[1]:
                    where = f"between {x[0]} and {y[0]}"
            sw = f"{a}-{b} ({where})"
    n_core = len(cnt)
    if reasons:
        cls = "nested_or_merged"
    elif n_core >= args.min_domains and sw and "outside" not in sw:
        cls = "crossover"
    elif n_core == 2 and sw and "outside" not in sw:
        cls = "crossover_partial"
    else:
        cls = "undetermined"
        if n_core < args.min_domains:
            notes.append(f"only {n_core} core domains")
        if not sw:
            notes.append("switch not placed")
    rows.append([ns, ">".join(path), cls, ";".join(reasons), sw,
                 ",".join(x[0] for x in d), " | ".join(notes)])

out = os.path.join(args.out_dir, "crossover_vs_nesting.tsv")
with open(out, "w") as f:
    f.write("ns_id\tmosaic\tclass\tnesting_evidence\tswitch_interval\tdomains\tnotes\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print(f"crossover-like mosaics: {len(rows)} -> " + ", ".join(f"{k} {v}" for k, v in Counter(r[2] for r in rows).most_common()))
print("nesting evidence:", dict(Counter(e for r in rows if r[3] for e in r[3].split(";"))))
pairs = defaultdict(Counter)
for r in rows:
    pairs["-".join(sorted(set(r[1].split(">"))))][r[2]] += 1
print("\nby pair:")
for p, c in sorted(pairs.items(), key=lambda x: -sum(x[1].values())):
    print(f"  {p:>7}: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))
sw = Counter(re.sub(r"^\d+-\d+ \(|\)$", "", r[4]) for r in rows if r[2] == "crossover")
print("\ncrossover switch location:", sw.most_common())
