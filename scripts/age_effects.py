#!/usr/bin/env python3
"""Element age (5'-3' LTR identity, Athilafinder) vs context, mosaicism, junctions.

Usage: age_effects.py --summary-dir DIR --centro element_centro_class.tsv \
           --mosaics kmer_recombinants.tsv --junctions junctions_full/elements.tsv --out-dir OUT
"""
import argparse, glob, os
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import lib_ids as L

ap = argparse.ArgumentParser()
for a in ("--summary-dir", "--centro", "--mosaics", "--junctions", "--out-dir"):
    ap.add_argument(a, required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

age = {}
for p in glob.glob(os.path.join(args.summary_dir, "*SUMMARY_TABLE.txt")):
    acc = L.accession_from_path(p)
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) > 12 and f[5] == "intact":
            v = [float(x) for x in f[9:13] if x not in ("NA", "")]
            if v:
                age[L.nsid(acc, f[3])] = float(np.mean(v))
cls, fam = {}, {}
for l in open(args.centro).readlines()[1:]:
    f = l.rstrip("\n").split("\t"); cls[f[0]] = f[5]; fam[f[0]] = f[3]
mos = {l.split("\t")[0] for l in open(args.mosaics).readlines()[1:]}
jrows = [l.rstrip("\n").split("\t") for l in open(args.junctions)]
jh = jrows[0]; J = {r[0]: dict(zip(jh, r)) for r in jrows[1:]}

ids = [i for i in age if i in cls]
A = np.array([age[i] for i in ids])
out = ["# Element age (5'-3' LTR identity; higher = younger)\n",
       f"- intact elements with identity + context: {len(ids)}\n",
       "## Age by centromeric context\n", "| class | n | median LTR identity |", "|---|---:|---:|"]
groups = defaultdict(list)
for i in ids:
    groups[cls[i]].append(age[i])
for k in ("centromeric", "pericentromeric", "arm"):
    out.append(f"| {k} | {len(groups[k])} | {np.median(groups[k]):.1f}% |")
u = mannwhitneyu(groups["centromeric"], groups["pericentromeric"] + groups["arm"])
out.append(f"\nCentromeric vs rest: Mann-Whitney p = {u.pvalue:.2g}")

out += ["\n## Recombination signals by age quintile (intact elements)\n",
        "| LTR identity bin | n | % mosaic (k-mer) | clean deletions per 100 | internal family switches per 100 |",
        "|---|---:|---:|---:|---:|"]
qs = np.quantile(A, [0, .2, .4, .6, .8, 1])
for lo, hi in zip(qs, qs[1:]):
    sel = [i for i in ids if lo <= age[i] <= hi]
    m = np.mean([i in mos for i in sel]) * 100
    dl = np.mean([int(J[i]["n_deletion_clean"]) for i in sel if i in J]) * 100
    sw = np.mean([int(J[i]["n_internal_family_switch"]) for i in sel if i in J]) * 100
    out.append(f"| {lo:.1f}-{hi:.1f}% | {len(sel)} | {m:.1f} | {dl:.2f} | {sw:.2f} |")
r = spearmanr(A, [i in mos for i in ids])
out.append(f"\nSpearman(LTR identity, mosaic) rho = {r.correlation:+.3f}, p = {r.pvalue:.2g}")
open(os.path.join(args.out_dir, "age_summary.md"), "w").write("\n".join(out) + "\n")
with open(os.path.join(args.out_dir, "element_age.tsv"), "w") as fh:
    fh.write("ns_id\tltr_identity\tclass\tfamily\tmosaic\n")
    for i in ids:
        fh.write(f"{i}\t{age[i]:.2f}\t{cls[i]}\t{fam[i]}\t{int(i in mos)}\n")
print("\n".join(out))
