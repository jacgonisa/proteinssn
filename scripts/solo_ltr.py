#!/usr/bin/env python3
"""Solo LTR : intact ratio by LTR family and centromeric context.

A solo LTR is what unequal recombination between an element's two LTRs leaves
behind. Its frequency relative to intact elements of the same family estimates
how much LTR-LTR (homologous) recombination has eliminated that family.

Families are assigned by aligning LTRs to the TAIR12 LTR exemplars. Exemplar LTRs
are shared by some sister families, so LTR families are groups:
  ATHILA6 = 6/6a/6b, ATHILA3 = 0/3, ATHILA7 = 7/7a.

Usage:
    solo_ltr.py --summary-dir DIR --solo-elements junctions_solo/elements.tsv \
        --full-elements junctions_full/elements.tsv --arrays cen178_arrays.tsv --out-dir OUT
"""
import argparse
import glob
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import fisher_exact

import lib_ids as L
from centro_class import cores

ap = argparse.ArgumentParser()
for a in ("--summary-dir", "--solo-elements", "--full-elements", "--arrays", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--peri-bp", type=int, default=2_000_000)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
core = cores(args.arrays, 500_000)


def ltr_fam(path):
    d = {}
    for l in open(path).readlines()[1:]:
        f = l.rstrip("\n").split("\t")
        d[f[0]] = f[4]
    return d


solo_fam, full_fam = ltr_fam(args.solo_elements), ltr_fam(args.full_elements)

rows = []
for p in glob.glob(os.path.join(args.summary_dir, "*SUMMARY_TABLE.txt")):
    acc = L.accession_from_path(p)
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7 or f[0] == "chr" or f[5] not in ("intact", "solo"):
            continue
        ns = L.nsid(acc, f[3])
        fam = (solo_fam if f[5] == "solo" else full_fam).get(ns, "")
        if not fam:
            continue
        s, e = int(f[1]), int(f[2])
        cls = "unplaced"
        if (acc, f[0]) in core:
            cs, ce, _ = core[(acc, f[0])]
            d = 0 if (e >= cs and s <= ce) else min(abs(s - ce), abs(cs - e))
            cls = "centromeric" if d == 0 else ("pericentromeric" if d <= args.peri_bp else "arm")
        ident = [float(x) for x in f[9:13] if x not in ("NA", "")]
        rows.append({"acc": acc, "kind": f[5], "fam": fam, "class": cls,
                     "ltr_ident": np.mean(ident) if ident else np.nan})

with open(os.path.join(args.out_dir, "solo_intact_elements.tsv"), "w") as fh:
    fh.write("accession\tkind\tltr_family\tclass\tltr_identity\n")
    for r in rows:
        fh.write(f"{r['acc']}\t{r['kind']}\t{r['fam']}\t{r['class']}\t{r['ltr_ident']:.2f}\n")

lines = ["# Solo LTR : intact ratio (LTR-LTR recombination)\n",
         f"- intact {sum(r['kind']=='intact' for r in rows)}, solo {sum(r['kind']=='solo' for r in rows)}\n",
         "## By LTR family group\n", "| LTR family | intact | solo | solo fraction |", "|---|---:|---:|---:|"]
by = defaultdict(Counter)
for r in rows:
    by[r["fam"]][r["kind"]] += 1
for fam, c in sorted(by.items(), key=lambda x: -(x[1]["solo"] / max(1, x[1]["solo"] + x[1]["intact"]))):
    n = c["solo"] + c["intact"]
    if n >= 50:
        lines.append(f"| {fam} | {c['intact']} | {c['solo']} | {c['solo']/n:.2f} |")
lines += ["\n## By centromeric context\n", "| class | intact | solo | solo fraction |", "|---|---:|---:|---:|"]
bc = defaultdict(Counter)
for r in rows:
    bc[r["class"]][r["kind"]] += 1
for k in ("centromeric", "pericentromeric", "arm"):
    c = bc[k]; n = c["solo"] + c["intact"]
    lines.append(f"| {k} | {c['intact']} | {c['solo']} | {c['solo']/max(1,n):.2f} |")
a, b = bc["centromeric"]["solo"], bc["centromeric"]["intact"]
c_ = bc["pericentromeric"]["solo"] + bc["arm"]["solo"]
d_ = bc["pericentromeric"]["intact"] + bc["arm"]["intact"]
orr, p = fisher_exact([[a, b], [c_, d_]])
lines.append(f"\nCentromeric vs rest: OR {orr:.2f}, Fisher p = {p:.2g}")
lines += ["\n## By family x context (solo fraction; n>=30 per cell)\n",
          "| LTR family | centromeric | pericentromeric | arm |", "|---|---:|---:|---:|"]
bfc = defaultdict(Counter)
for r in rows:
    bfc[(r["fam"], r["class"])][r["kind"]] += 1
for fam in sorted(by, key=lambda f: -sum(by[f].values())):
    cells = []
    for k in ("centromeric", "pericentromeric", "arm"):
        c = bfc[(fam, k)]; n = c["solo"] + c["intact"]
        cells.append(f"{c['solo']/n:.2f} (n={n})" if n >= 30 else "–")
    if sum(by[fam].values()) >= 50:
        lines.append(f"| {fam} | " + " | ".join(cells) + " |")
open(os.path.join(args.out_dir, "summary.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
