#!/usr/bin/env python3
"""Crossover-like vs conversion-like mosaics from ordered k-mer marker runs.

For each element (internal domain span), family-specific marker hits are ordered by
position and compressed into runs of the same family. Runs shorter than --min-run
markers are treated as noise and dropped, adjacent same-family runs are merged, and
the resulting family path is classified:

  A>B        crossover-like (one switch)
  A>B>A      conversion-like (a tract of B inside A)
  other      complex (>=2 switches through >=3 families, or repeated returns)

A shuffle null permutes family labels across the element's hit positions (keeps
composition, destroys spatial clustering) to estimate chance calls of each class.

Usage:
    tract_scan.py --markers robust_markers.tsv --fulllength DIR --clip-dir DIR \
        --out-dir OUT [--min-run 5 8 12] [--null-reps 5]
"""

from __future__ import annotations

import argparse
import glob
import os
import random
from collections import Counter, defaultdict

import lib_ids as L
from athila_recomb_nt import read_fasta, revcomp


def load_markers(path):
    m = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            km, fam = line.rstrip("\n").split("\t")
            km = km.upper()
            m[km] = fam[6:] if fam.startswith("ATHILA") else fam
            m[revcomp(km)] = m[km]
    return m


def load_spans(clip_dir):
    span = {}
    for gff in glob.glob(os.path.join(clip_dir, "*.dom.gff3")):
        acc = L.accession_from_path(gff)
        for line in open(gff):
            f = line.split("\t")
            if len(f) < 9:
                continue
            ns = L.nsid(acc, f[0])
            a, b = int(f[3]), int(f[4])
            lo, hi = span.get(ns, (a, b))
            span[ns] = (min(lo, a), max(hi, b))
    return span


def path_of(hits, min_run):
    """hits: [(pos, fam)] sorted. Return [(fam, n, start, end)] after denoising."""
    runs = []
    for pos, f in hits:
        if runs and runs[-1][0] == f:
            runs[-1][1] += 1; runs[-1][3] = pos
        else:
            runs.append([f, 1, pos, pos])
    runs = [r for r in runs if r[1] >= min_run]
    merged = []
    for r in runs:
        if merged and merged[-1][0] == r[0]:
            merged[-1][1] += r[1]; merged[-1][3] = r[3]
        else:
            merged.append(list(r))
    return merged


def classify(path):
    fams = [p[0] for p in path]
    if len(fams) < 2:
        return None
    if len(fams) == 2:
        return "crossover"
    if len(fams) == 3 and fams[0] == fams[2]:
        return "conversion"
    return "complex"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", required=True)
    ap.add_argument("--fulllength", required=True)
    ap.add_argument("--clip-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--min-run", type=int, nargs="+", default=[5, 8, 12])
    ap.add_argument("--null-reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "tract_scan")
    os.makedirs(out, exist_ok=True)
    rng = random.Random(args.seed)
    marker = load_markers(args.markers)
    span = load_spans(args.clip_dir)

    obs = {m: Counter() for m in args.min_run}
    null = {m: Counter() for m in args.min_run}
    calls = []                                   # at the middle stringency
    mid = sorted(args.min_run)[len(args.min_run) // 2]
    n_el = 0
    for flf in sorted(glob.glob(os.path.join(args.fulllength,
                                             "*DP_FULLLENGTH_RENAMED.fasta"))):
        acc = L.accession_from_path(flf)
        for el, seq in read_fasta(flf).items():
            ns = L.nsid(acc, el)
            if ns not in span:
                continue
            lo, hi = span[ns]
            s = seq[lo - 1:hi].upper()
            hits = [(i, marker[s[i:i + args.k]]) for i in range(len(s) - args.k + 1)
                    if s[i:i + args.k] in marker]
            if len(hits) < 20:
                continue
            n_el += 1
            labels = [f for _, f in hits]
            for m in args.min_run:
                p = path_of(hits, m)
                c = classify(p)
                if c:
                    obs[m][c] += 1
                    if m == mid:
                        calls.append((ns, c, p, lo))
                for _ in range(args.null_reps):
                    sh = labels[:]
                    rng.shuffle(sh)
                    cn = classify(path_of(list(zip([h[0] for h in hits], sh)), m))
                    if cn:
                        null[m][cn] += 1 / args.null_reps

    with open(os.path.join(out, "tract_calls.tsv"), "w") as f:
        f.write("ns_id\tchrom\tstart\tclass\tpath\ttract_family\ttract_markers\t"
                "tract_bp\tswitch_pos_in_span\n")
        for ns, c, p, lo in calls:
            acc, el = L.split_nsid(ns)
            co = L.element_coords(el) or {}
            pstr = ">".join(f"{x[0]}({x[1]})" for x in p)
            tf = tm = tbp = ""
            if c == "conversion":
                tf, tm, tbp = p[1][0], p[1][1], p[1][3] - p[1][2] + args.k
            sw = ";".join(str(p[i][3]) for i in range(len(p) - 1))
            f.write(f"{ns}\t{co.get('chrom','')}\t{co.get('start','')}\t{c}\t{pstr}\t"
                    f"{tf}\t{tm}\t{tbp}\t{sw}\n")

    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write(f"# Crossover-like vs conversion-like (elements typed: {n_el})\n\n")
        f.write("| min run | class | observed | shuffle null | obs/null |\n"
                "|---:|---|---:|---:|---:|\n")
        for m in args.min_run:
            for c in ["crossover", "conversion", "complex"]:
                o, nl = obs[m][c], null[m][c]
                ratio = f"{o / nl:.1f}" if nl else "inf"
                f.write(f"| {m} | {c} | {o} | {nl:.1f} | {ratio} |\n")
    print(open(os.path.join(out, "summary.md")).read())


if __name__ == "__main__":
    main()
