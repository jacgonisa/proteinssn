#!/usr/bin/env python3
"""Crossover- vs conversion-like mosaics counted in INDEPENDENT variant sites.

Overlapping marker k-mers that share one variant are collapsed into an "island"
(consecutive hit positions); each island carries one family label and counts as one
independent site. Paths are built over islands; a run must contain >= --min-sites
islands to count. Null: shuffle family labels among islands within an element.

Usage: island_scan.py --markers M --fulllength DIR --clip-dir DIR --out-dir OUT
"""
import argparse, glob, os, random
from collections import Counter
import lib_ids as L
from athila_recomb_nt import read_fasta
from tract_scan import load_markers, load_spans, classify


def islands(hits):
    isl = []
    for pos, f in hits:
        if isl and pos == isl[-1][3] + 1 and f == isl[-1][0]:
            isl[-1][3] = pos
        else:
            isl.append([f, 1, pos, pos])
    return [(f, a, b) for f, _, a, b in isl]


def path_sites(isl, min_sites):
    runs = []
    for f, a, b in isl:
        if runs and runs[-1][0] == f:
            runs[-1][1] += 1; runs[-1][3] = b
        else:
            runs.append([f, 1, a, b])
    runs = [r for r in runs if r[1] >= min_sites]
    out = []
    for r in runs:
        if out and out[-1][0] == r[0]:
            out[-1][1] += r[1]; out[-1][3] = r[3]
        else:
            out.append(list(r))
    return out


ap = argparse.ArgumentParser()
for a in ("--markers", "--fulllength", "--clip-dir", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--k", type=int, default=31)
ap.add_argument("--min-sites", type=int, nargs="+", default=[2, 3, 5])
ap.add_argument("--null-reps", type=int, default=10)
args = ap.parse_args()
rng = random.Random(0)
marker = load_markers(args.markers); span = load_spans(args.clip_dir)
obs = {m: Counter() for m in args.min_sites}; null = {m: Counter() for m in args.min_sites}
calls = {m: [] for m in args.min_sites}; n = 0
for flf in sorted(glob.glob(os.path.join(args.fulllength, "*DP_FULLLENGTH_RENAMED.fasta"))):
    acc = L.accession_from_path(flf)
    for el, seq in read_fasta(flf).items():
        ns = L.nsid(acc, el)
        if ns not in span:
            continue
        lo, hi = span[ns]; s = seq[lo - 1:hi].upper()
        hits = [(i, marker[s[i:i + args.k]]) for i in range(len(s) - args.k + 1)
                if s[i:i + args.k] in marker]
        isl = islands(hits)
        if len(isl) < 6:
            continue
        n += 1
        labs = [x[0] for x in isl]
        for m in args.min_sites:
            p = path_sites(isl, m); c = classify(p)
            if c:
                obs[m][c] += 1; calls[m].append((ns, c, p))
            for _ in range(args.null_reps):
                sh = labs[:]; rng.shuffle(sh)
                cn = classify(path_sites([(f, a, b) for f, (_, a, b) in zip(sh, isl)], m))
                if cn:
                    null[m][cn] += 1 / args.null_reps
os.makedirs(os.path.join(args.out_dir, "island_scan"), exist_ok=True)
with open(os.path.join(args.out_dir, "island_scan", "calls.tsv"), "w") as f:
    f.write("min_sites\tns_id\tclass\tpath_sites\ttract_family\ttract_sites\ttract_bp\n")
    for m in args.min_sites:
        for ns, c, p in calls[m]:
            pstr = ">".join(f"{x[0]}({x[1]})" for x in p)
            tf = ts = tb = ""
            if c == "conversion":
                tf, ts, tb = p[1][0], p[1][1], p[1][3] - p[1][2] + args.k
            f.write(f"{m}\t{ns}\t{c}\t{pstr}\t{tf}\t{ts}\t{tb}\n")
lines = [f"# Mosaic classes counted in independent variant sites (elements: {n})\n",
         "| min sites per run | class | observed | island-shuffle null | obs/null |",
         "|---:|---|---:|---:|---:|"]
for m in args.min_sites:
    for c in ("crossover", "conversion", "complex"):
        o, nl = obs[m][c], null[m][c]
        lines.append(f"| {m} | {c} | {o} | {nl:.1f} | {o/nl:.1f} |" if nl else f"| {m} | {c} | {o} | 0 | inf |")
open(os.path.join(args.out_dir, "island_scan", "summary.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
