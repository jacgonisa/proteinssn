#!/usr/bin/env python3
"""Power simulation: plant B tracts (conversion) or B 3' ends (crossover) into
clean family-A elements and measure recovery by the independent-site classifier."""
import glob, os, random, sys
from collections import defaultdict
import lib_ids as L
from athila_recomb_nt import read_fasta
from tract_scan import load_markers, load_spans, classify
from island_scan_lib import islands, path_sites

R, FL, CLIP = sys.argv[1:4]
K, MIN_SITES, REPS = 31, 3, 60
rng = random.Random(1)
marker = load_markers(f"{R}/kmc_robust/robust_markers.tsv"); span = load_spans(CLIP)
mos = {l.split("\t")[0] for l in open(f"{R}/kmer_recomb_robust_internal/kmer_recombinants.tsv").readlines()[1:]}
pure = defaultdict(list)
for l in open(f"{R}/kmer_recomb_robust_internal/primary_family.tsv").readlines()[1:]:
    ns, f = l.rstrip().split("\t")
    if ns not in mos and ns in span: pure[f[6:]].append(ns)
want = set()
pairs = [("6a", "6b"), ("6", "1"), ("2", "6b"), ("2", "4c"), ("4", "4c"), ("1", "2")]
pick = {}
for a, b in pairs:
    for f in (a, b):
        if f not in pick:
            pick[f] = rng.sample(pure[f], min(40, len(pure[f])))
            want |= set(pick[f])
seqs = {}
for flf in glob.glob(f"{FL}/*DP_FULLLENGTH_RENAMED.fasta"):
    acc = L.accession_from_path(flf)
    for el, s in read_fasta(flf).items():
        ns = L.nsid(acc, el)
        if ns in want:
            lo, hi = span[ns]; seqs[ns] = s[lo - 1:hi].upper()

def call(s):
    hits = [(i, marker[s[i:i+K]]) for i in range(len(s)-K+1) if s[i:i+K] in marker]
    return classify(path_sites(islands(hits), MIN_SITES))

print("pair  host>donor  | crossover recall | conversion recall by tract length (bp)")
print("                 |                  |  100   250   500  1000  2000")
for a, b in pairs:
    base_fp = sum(call(seqs[x]) is not None for x in pick[a]) / len(pick[a])
    xo = 0; conv = {L_: 0 for L_ in (100, 250, 500, 1000, 2000)}
    for _ in range(REPS):
        A = seqs[rng.choice(pick[a])]; B = seqs[rng.choice(pick[b])]
        cut = int(len(A) * rng.uniform(0.35, 0.65))
        xo += call(A[:cut] + B[int(len(B) * cut / len(A)):]) == "crossover"
        for Lt in conv:
            st = int(len(A) * rng.uniform(0.25, 0.6)); sb = int(len(B) * st / len(A))
            conv[Lt] += call(A[:st] + B[sb:sb + Lt] + A[st + Lt:]) == "conversion"
    print(f"{a:>3}>{b:<3} (FP {base_fp:.0%}) | {xo/REPS:>8.0%}         | " +
          " ".join(f"{conv[x]/REPS:>4.0%}" for x in conv))
