#!/usr/bin/env python3
"""All ATHILA-derived copies: census, chopping by family/context, truncation hotspots.

Inputs:
  copies.tsv               genome_copies.py (full_length / solo_ltr / fragment_*)
  junctions_frag/          junction_scan.py on fragment sequences (aln.paf, elements.tsv,
                           junctions.tsv)
  junctions_full/aln.paf   full-length elements vs exemplars (for the domain map)
  TEsorter dir             dom.gff3 (domain coordinates on full-length elements)
  CEN178 arrays            centromere cores

Outputs: census tables, per-family / per-context chopping, truncation-breakpoint
positions on each family's internal consensus with projected domain positions.

Usage:
    fragment_analysis.py --copies copies.tsv --frag-dir junctions_frag \
        --full-paf junctions_full/aln.paf --tesorter DIR --arrays cen178_arrays.tsv \
        --out-dir OUT
"""
import argparse
import glob
import os
import re
from collections import Counter, defaultdict

import numpy as np

import lib_ids as L
from centro_class import cores

ap = argparse.ArgumentParser()
for a in ("--copies", "--frag-dir", "--full-paf", "--tesorter", "--arrays", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--peri-bp", type=int, default=2_000_000)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
core = cores(args.arrays, 500_000)


def short(f):
    return f[6:] if f.startswith("ATHILA") else f


def context(acc, chrom, s, e):
    if (acc, chrom) not in core:
        return "unplaced"
    cs, ce, _ = core[(acc, chrom)]
    d = 0 if (e >= cs and s <= ce) else min(abs(s - ce), abs(cs - e))
    return "centromeric" if d == 0 else ("pericentromeric" if d <= args.peri_bp else "arm")


# ---- census ---------------------------------------------------------------- #
rows = [l.rstrip("\n").split("\t") for l in open(args.copies)]
H = rows[0]
C = [dict(zip(H, r)) for r in rows[1:]]
for c in C:
    c["ctx"] = context(c["accession"], c["chrom"], int(c["start"]), int(c["end"]))
# collapse full-length pieces to Athilafinder elements; solo likewise
seen = set()
units = []
for c in C:
    if c["label"] in ("full_length", "solo_ltr"):
        key = (c["accession"], c["athilafinder_id"])
        if key in seen:
            continue
        seen.add(key)
    units.append(c)
lab = Counter(u["label"] for u in units)
accs = sorted({u["accession"] for u in units})
per_acc = {k: np.median([sum(1 for u in units if u["accession"] == a and u["label"] == k)
                         for a in accs]) for k in lab}

out = ["# All ATHILA-derived copies\n",
       f"- genomes: {len(accs)}; copies: {len(units)} "
       f"(full-length pieces collapsed to Athilafinder elements)\n",
       "| copy type | total | median per genome |", "|---|---:|---:|"]
for k in ("full_length", "solo_ltr", "fragment_internal", "fragment_ltr"):
    out.append(f"| {k} | {lab[k]} | {per_acc.get(k, 0):.0f} |")

# chopping by internal family: fragments vs intact
out += ["\n## Internal fragments vs intact elements, by family\n",
        "| family | intact | internal fragments | fragments per intact |", "|---|---:|---:|---:|"]
fam_c = defaultdict(Counter)
for u in units:
    if u["label"] in ("full_length", "fragment_internal") and u["internal_family"]:
        fam_c[short(u["internal_family"])][u["label"]] += 1
for fam, c in sorted(fam_c.items(), key=lambda x: -x[1]["fragment_internal"] / max(1, x[1]["full_length"])):
    if c["full_length"] + c["fragment_internal"] >= 50:
        out.append(f"| {fam} | {c['full_length']} | {c['fragment_internal']} | "
                   f"{c['fragment_internal'] / max(1, c['full_length']):.2f} |")
out += ["\n## By centromeric context\n",
        "| context | intact | solo LTR | internal fragments | LTR fragments | fragments per intact |",
        "|---|---:|---:|---:|---:|---:|"]
ctx_c = defaultdict(Counter)
for u in units:
    ctx_c[u["ctx"]][u["label"]] += 1
for k in ("centromeric", "pericentromeric", "arm"):
    c = ctx_c[k]
    out.append(f"| {k} | {c['full_length']} | {c['solo_ltr']} | {c['fragment_internal']} | "
               f"{c['fragment_ltr']} | {c['fragment_internal'] / max(1, c['full_length']):.2f} |")

with open(os.path.join(args.out_dir, "copies_with_context.tsv"), "w") as fh:
    fh.write("\t".join(H + ["context"]) + "\n")
    for u in units:
        fh.write("\t".join(u[h] for h in H) + f"\t{u['ctx']}\n")

# ---- domain map on each internal exemplar ----------------------------------- #
gff = defaultdict(list)
for p in glob.glob(os.path.join(args.tesorter, "*.dom.gff3")):
    acc = L.accession_from_path(p)
    for line in open(p):
        f = line.split("\t")
        if len(f) < 9:
            continue
        g = re.search(r"gene=([A-Za-z]+)", f[8])
        if g and g.group(1) in L.CORE_DOMAINS:
            gff[L.nsid(acc, f[0])].append((g.group(1), (int(f[3]) + int(f[4])) // 2))
dom_pos = defaultdict(lambda: defaultdict(list))          # ref -> domain -> rel pos
for line in open(args.full_paf):
    f = line.split("\t")
    q, qs, qe, strand, r, rlen, rs, re_ = f[0], int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8])
    if "_LTR" in r or q not in gff:
        continue
    for d, m in gff[q]:
        if qs <= m <= qe:
            rp = rs + (m - qs) if strand == "+" else re_ - (m - qs)
            dom_pos[r][d].append(rp / rlen)
domain_map = {r: {d: float(np.median(v)) for d, v in ds.items() if len(v) >= 10}
              for r, ds in dom_pos.items()}

# ---- truncation breakpoints of internal fragments --------------------------- #
ends = defaultdict(list)                                  # ref -> rel positions of fragment ends
fpaf = os.path.join(args.frag_dir, "aln.paf")
by_q = defaultdict(list)
for line in open(fpaf):
    f = line.split("\t")
    if "_LTR" in f[5] or int(f[9]) < 80:
        continue
    by_q[f[0]].append((int(f[9]), f[5], int(f[6]), int(f[7]), int(f[8])))
for q, hs in by_q.items():
    best = Counter()
    for nm, r, rlen, rs, re_ in hs:
        best[(r, rlen)] += nm
    (r, rlen), _ = best.most_common(1)[0]
    iv = [(rs, re_) for nm, rr, rl, rs, re_ in hs if rr == r]
    lo, hi = min(s for s, e in iv), max(e for s, e in iv)
    if lo > 0.02 * rlen:
        ends[r].append(lo / rlen)
    if hi < 0.98 * rlen:
        ends[r].append(hi / rlen)
with open(os.path.join(args.out_dir, "truncation_breakpoints.tsv"), "w") as fh:
    fh.write("exemplar\tfamily\trel_pos\n")
    for r, v in ends.items():
        for x in v:
            m = re.match(r"(ATHILA\d+[abc]?)", r)
            fh.write(f"{r}\t{short(m.group(1)) if m else r}\t{x:.4f}\n")
with open(os.path.join(args.out_dir, "domain_map.tsv"), "w") as fh:
    fh.write("exemplar\tdomain\trel_pos\n")
    for r, ds in domain_map.items():
        for d, x in ds.items():
            fh.write(f"{r}\t{d}\t{x:.4f}\n")

# where do breakpoints fall relative to domains (pooled, by nearest domain-bounded segment)
seg = Counter()
for r, v in ends.items():
    dm = domain_map.get(r)
    if not dm or len(dm) < 4:
        continue
    order = sorted(dm.items(), key=lambda x: x[1])
    for x in v:
        lab = "5' of " + order[0][0]
        for (d1, p1), (d2, p2) in zip(order, order[1:]):
            if p1 <= x < p2:
                lab = f"{d1}-{d2}"
        if x >= order[-1][1]:
            lab = "3' of " + order[-1][0]
        seg[lab] += 1
out += ["\n## Where internal fragments break (truncation ends on the family consensus)\n",
        f"- fragment ends analysed: {sum(len(v) for v in ends.values())}",
        "- interval between projected domain centres (pooled over exemplars with >=4 mapped domains):\n",
        "| interval | breakpoints |", "|---|---:|"]
for k, n in seg.most_common():
    out.append(f"| {k} | {n} |")

# ---- junctions inside fragments vs intact ----------------------------------- #
def jrate(path):
    r = [l.rstrip("\n").split("\t") for l in open(path)]
    h = r[0]; E = [dict(zip(h, x)) for x in r[1:]]
    keys = [k for k in h if k.startswith("n_")]
    return len(E), {k: np.mean([int(e[k]) for e in E]) * 100 for k in keys}
nf, rf = jrate(os.path.join(args.frag_dir, "elements.tsv"))
out += ["\n## Junctions per 100 copies: internal fragments vs intact elements\n",
        "| junction class | intact | fragments |", "|---|---:|---:|"]
full_el = os.path.join(os.path.dirname(args.full_paf), "elements.tsv")
ni, ri = jrate(full_el)
for k in rf:
    out.append(f"| {k[2:]} | {ri.get(k, 0):.2f} | {rf[k]:.2f} |")
out.append(f"\n(n intact = {ni}, n fragments = {nf})")
open(os.path.join(args.out_dir, "summary.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
