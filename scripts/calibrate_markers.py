#!/usr/bin/env python3
"""Calibrate k-mer marker settings for crossover detection (k-mers only).

Clean (non-mosaic) elements of each family are split 50/50. Markers are built with
KMC from the TRAIN half (k-mers in >= prevalence of the family, excluded if in >=
--other-prevalence of any other family). On the TEST half:
  false-positive rate = fraction of clean held-out elements called mosaic
  recall              = fraction of synthetic crossovers (A 5' part + B 3' part at a
                        random breakpoint 0.3-0.7) recovered with both families
Calls use independent marker sites (islands) and a minimum number of sites per
segment (--min-sites grid).

Usage: calibrate_markers.py --primary primary_family.tsv --mosaics old_mosaics.tsv \
          --internal DIR --out-dir OUT --kmc-bin DIR
"""
import argparse, glob, math, os, random, subprocess
from collections import defaultdict
import lib_ids as L
from athila_recomb_nt import read_fasta, revcomp
from island_scan_lib import islands, path_sites
from tract_scan import classify

ap = argparse.ArgumentParser()
for a in ("--primary", "--mosaics", "--internal", "--out-dir", "--kmc-bin"):
    ap.add_argument(a, required=True)
ap.add_argument("--prevalence", type=float, nargs="+", default=[0.05, 0.10, 0.25])
ap.add_argument("--other-prevalence", type=float, default=0.05)
ap.add_argument("--min-sites", type=int, nargs="+", default=[2, 3, 5, 8])
ap.add_argument("--min-members", type=int, default=30)
ap.add_argument("--n-test", type=int, default=150)
ap.add_argument("--n-planted", type=int, default=60)
ap.add_argument("--k", type=int, default=31)
args = ap.parse_args()
rng = random.Random(0)
work = os.path.join(args.out_dir, "calib"); os.makedirs(os.path.join(work, "tmp"), exist_ok=True)
kmc, tools = os.path.join(args.kmc_bin, "kmc"), os.path.join(args.kmc_bin, "kmc_tools")

mos = {l.split("\t")[0] for l in open(args.mosaics).readlines()[1:]}
fam_of = {}
for l in open(args.primary).readlines()[1:]:
    ns, f = l.rstrip("\n").split("\t")
    if ns not in mos:
        fam_of[ns] = f.replace("ATHILA", "")
seqs = {}
for p in glob.glob(os.path.join(args.internal, "*.fasta")):
    acc = L.accession_from_path(p)
    for el, s in read_fasta(p).items():
        ns = L.nsid(acc, el)
        if ns in fam_of:
            seqs[ns] = s.upper()
by_fam = defaultdict(list)
for ns, f in fam_of.items():
    if ns in seqs:
        by_fam[f].append(ns)
fams = sorted(f for f, v in by_fam.items() if len(v) >= args.min_members)
train, test = {}, {}
for f in fams:
    v = by_fam[f][:]; rng.shuffle(v); h = len(v) // 2
    train[f], test[f] = v[:h], v[h:]
print("families calibrated:", " ".join(f"{f}({len(train[f])}/{len(test[f])})" for f in fams))


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


for f in fams:
    with open(os.path.join(work, f"{f}.fa"), "w") as o:
        for ns in train[f]:
            o.write(f">{ns}\n{seqs[ns]}\n")


def build(prev):
    for f in fams:
        n = len(train[f])
        run([kmc, f"-k{args.k}", f"-ci{max(2, math.ceil(prev * n))}", "-cs1000000", "-t16", "-fm",
             os.path.join(work, f"{f}.fa"), os.path.join(work, f"{f}.in"), os.path.join(work, "tmp")])
        run([kmc, f"-k{args.k}", f"-ci{max(2, math.ceil(args.other_prevalence * n))}", "-cs1000000", "-t16", "-fm",
             os.path.join(work, f"{f}.fa"), os.path.join(work, f"{f}.pres"), os.path.join(work, "tmp")])
    marker = {}
    for f in fams:
        others = [os.path.join(work, f"{g}.pres") for g in fams if g != f]
        acc_db = others[0]
        for i, o in enumerate(others[1:]):
            nxt = os.path.join(work, f"u{i}")
            run([tools, "simple", acc_db, o, "union", nxt]); acc_db = nxt
        run([tools, "simple", os.path.join(work, f"{f}.in"), acc_db, "kmers_subtract", os.path.join(work, "mk")])
        run([tools, "transform", os.path.join(work, "mk"), "dump", os.path.join(work, "mk.txt")])
        for line in open(os.path.join(work, "mk.txt")):
            km = line.split()[0]
            marker[km] = f; marker[revcomp(km)] = f
    return marker


def call(seq, marker, m):
    hits = [(i, marker[seq[i:i + args.k]]) for i in range(len(seq) - args.k + 1) if seq[i:i + args.k] in marker]
    p = path_sites(islands(hits), m)
    return classify(p), {x[0] for x in p}


test_set = {f: test[f][:args.n_test] for f in fams}
planted = []
pairs = [(a, b) for a in fams for b in fams if a != b]
for _ in range(args.n_planted * 4):
    a, b = rng.choice(pairs)
    A, B = seqs[rng.choice(test[a])], seqs[rng.choice(test[b])]
    x = rng.uniform(0.3, 0.7)
    planted.append((a, b, A[:int(len(A) * x)] + B[int(len(B) * x):]))

out = []
for prev in args.prevalence:
    marker = build(prev)
    nmark = len(marker) // 2
    for m in args.min_sites:
        fp = sum(call(seqs[ns], marker, m)[0] is not None for f in fams for ns in test_set[f])
        ntest = sum(len(v) for v in test_set.values())
        rec = sum(1 for a, b, s in planted if (lambda r: r[0] is not None and {a, b} <= r[1])(call(s, marker, m)))
        out.append((prev, m, nmark, fp / ntest, rec / len(planted)))
        print(f"prevalence {prev:.2f}  min sites {m}: markers {nmark:>7}  false-positive {fp/ntest*100:5.2f}%  "
              f"recall {rec/len(planted)*100:5.1f}%", flush=True)
with open(os.path.join(args.out_dir, "calibration.tsv"), "w") as f:
    f.write("prevalence\tmin_sites\tmarkers\tfalse_positive_rate\trecall\n")
    for r in out:
        f.write("\t".join(map(str, r)) + "\n")


# --------------------------------------------------------------------------- #
# Changepoint likelihood-ratio caller (appended): uses all sites, tolerant of noise
# --------------------------------------------------------------------------- #
import numpy as np


def llr_switch(seq, marker, k=31, min_side=3):
    """Max log-likelihood ratio of a single A->B changepoint vs no change.
    Returns (llr, famA, famB, breakpoint_position) for the two most frequent families."""
    hits = [(i, marker[seq[i:i + k]]) for i in range(len(seq) - k + 1) if seq[i:i + k] in marker]
    isl = islands(hits)
    labs = [f for f, a, b in isl]
    if len(labs) < 2 * min_side:
        return 0.0, None, None, None
    from collections import Counter
    top = [f for f, _ in Counter(labs).most_common(2)]
    if len(top) < 2:
        return 0.0, top[0], None, None
    xs = [(1 if f == top[0] else 0, isl[i][1]) for i, f in enumerate(labs) if f in top]
    x = np.array([v for v, _ in xs], float)
    n = len(x)

    def ll(v):
        if len(v) == 0:
            return 0.0
        p = min(max(v.mean(), 1e-6), 1 - 1e-6)
        return float(np.sum(v * np.log(p) + (1 - v) * np.log(1 - p)))
    base = ll(x)
    best, bt = 0.0, None
    cs = np.cumsum(x)
    for t in range(min_side, n - min_side + 1):
        l, r = x[:t], x[t:]
        p1, p2 = l.mean(), r.mean()
        if not ((p1 > 0.5 and p2 < 0.5) or (p1 < 0.5 and p2 > 0.5)):
            continue
        v = ll(l) + ll(r) - base
        if v > best:
            best, bt = v, t
    if bt is None:
        return 0.0, top[0], top[1], None
    first = top[0] if x[:bt].mean() > 0.5 else top[1]
    second = top[1] if first == top[0] else top[0]
    return best, first, second, xs[bt][1]


if __name__ == "__main__" and os.environ.get("LRT_CALIB"):
    for prev in [float(x) for x in os.environ["LRT_CALIB"].split(",")]:
        marker = build(prev)
        neg = [llr_switch(seqs[ns], marker)[0] for f in fams for ns in test_set[f]]
        pos = [(llr_switch(s, marker), a, b) for a, b, s in planted]
        for q in (0.99, 0.995, 0.999):
            thr = float(np.quantile(neg, q))
            rec = np.mean([r[0] > thr and {r[1], r[2]} == {a, b} for r, a, b in pos])
            print(f"LRT prevalence {prev:.2f}: threshold {thr:6.1f} (FPR {100*(1-q):.1f}%) -> recall {rec*100:5.1f}%", flush=True)
        with open(os.path.join(args.out_dir, f"lrt_null_{prev}.txt"), "w") as fh:
            fh.write("\n".join(f"{v:.3f}" for v in neg))
