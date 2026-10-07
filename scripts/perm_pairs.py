#!/usr/bin/env python3
"""Which ATHILA family pairs recombine more (or less) than expected?

Each mosaic contributes one event per adjacent family switch in its path (A>B>C
gives A-B and B-C). Two units are analysed: every mosaic element, and distinct loci
(same pair, chromosome and --locus-bp window) so that one inherited insertion
copied across many accessions counts once.

Two family-label permutation nulls:
  abundance  each event's two families are redrawn in proportion to how common the
             families are among all typed elements (controls for common families).
  degree     double-edge swaps on the event graph keep how many events each family
             takes part in but shuffle WHO it pairs with (controls for families that
             are simply recombination-prone).

Reports per pair: observed, expected, O/E, one-sided empirical p (enrichment and
depletion) and Benjamini-Hochberg q.

Usage:
    perm_pairs.py --recomb kmer_recombinants.tsv --primary primary_family.tsv \
        --out-dir OUT [--perms 10000] [--locus-bp 500000]
"""

from __future__ import annotations

import argparse
import os
import random
from collections import Counter

import numpy as np


def short(f):
    return f[6:] if f.startswith("ATHILA") else f


def load_events(path, locus_bp):
    elems, loci = [], set()
    for line in open(path).readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        chrom, start = f[2], f[3]
        fams = [s.split("(")[0] for s in f[6].split(">")]
        lb = int(start) // locus_bp if start.isdigit() else None
        for a, b in zip(fams, fams[1:]):
            if a != b:
                p = tuple(sorted((a, b)))
                elems.append(p)
                loci.add((p, chrom, lb))
    return elems, [x[0] for x in loci]


def bh(pvals):
    p = np.asarray(pvals)
    order = np.argsort(p)
    q = np.empty_like(p)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = len(p) - rank
        prev = min(prev, p[i] * len(p) / k)
        q[i] = prev
    return q


def null_abundance(n_events, fams, probs, rng):
    a = rng.choice(len(fams), size=n_events, p=probs)
    b = rng.choice(len(fams), size=n_events, p=probs)
    keep = a != b
    while not keep.all():                                  # redraw self pairs
        m = ~keep
        b[m] = rng.choice(len(fams), size=m.sum(), p=probs)
        keep = a != b
    return Counter(tuple(sorted((fams[i], fams[j]))) for i, j in zip(a, b))


def null_degree(events, rng, swaps_per_edge=10):
    e = [list(x) for x in events]
    n = len(e)
    for _ in range(swaps_per_edge * n):
        i, j = rng.randrange(n), rng.randrange(n)
        if i == j:
            continue
        a, b = e[i]
        c, d = e[j]
        if rng.random() < 0.5:
            c, d = d, c
        if a != d and c != b:
            e[i], e[j] = [a, d], [c, b]
    return Counter(tuple(sorted(x)) for x in e)


def test(events, ab, perms, seed, label):
    rng_np = np.random.default_rng(seed)
    rng = random.Random(seed)
    obs = Counter(events)
    fams = sorted(ab)
    probs = np.array([ab[f] for f in fams], float)
    probs /= probs.sum()
    pairs = sorted({tuple(sorted((a, b))) for i, a in enumerate(fams)
                    for b in fams[i + 1:]})
    res = {}
    for name in ("abundance", "degree"):
        sims = np.zeros((perms, len(pairs)))
        for k in range(perms):
            c = (null_abundance(len(events), fams, probs, rng_np) if name == "abundance"
                 else null_degree(events, rng))
            sims[k] = [c.get(p, 0) for p in pairs]
        o = np.array([obs.get(p, 0) for p in pairs])
        exp = sims.mean(0)
        p_hi = ((sims >= o).sum(0) + 1) / (perms + 1)
        p_lo = ((sims <= o).sum(0) + 1) / (perms + 1)
        res[name] = (o, exp, p_hi, p_lo, bh(p_hi), bh(p_lo))
    rows = []
    for idx, p in enumerate(pairs):
        r = {"unit": label, "pair": f"{p[0]}-{p[1]}", "observed": int(res["abundance"][0][idx])}
        for name in ("abundance", "degree"):
            o, exp, p_hi, p_lo, q_hi, q_lo = res[name]
            r[f"exp_{name}"] = exp[idx]
            r[f"oe_{name}"] = o[idx] / exp[idx] if exp[idx] > 0 else float("nan")
            r[f"p_enrich_{name}"] = p_hi[idx]
            r[f"q_enrich_{name}"] = q_hi[idx]
            r[f"p_deplete_{name}"] = p_lo[idx]
            r[f"q_deplete_{name}"] = q_lo[idx]
        rows.append(r)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recomb", required=True)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--locus-bp", type=int, default=500_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    ab = Counter(short(l.rstrip("\n").split("\t")[1])
                 for l in open(args.primary).readlines()[1:])
    elems, loci = load_events(args.recomb, args.locus_bp)
    for p in set(elems):
        for f in p:
            ab.setdefault(f, 1)
    rows = (test(elems, ab, args.perms, args.seed, "elements")
            + test(loci, ab, args.perms, args.seed + 1, "loci"))

    cols = ["unit", "pair", "observed",
            "exp_abundance", "oe_abundance", "q_enrich_abundance", "q_deplete_abundance",
            "exp_degree", "oe_degree", "q_enrich_degree", "q_deplete_degree",
            "p_enrich_abundance", "p_deplete_abundance", "p_enrich_degree",
            "p_deplete_degree"]
    out = os.path.join(args.out_dir, "pair_permutation.tsv")
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["unit"], -r["observed"])):
            f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")
    print(f"events: {len(elems)} element-level, {len(loci)} distinct loci")
    for unit in ("elements", "loci"):
        print(f"\n== {unit}: pairs significant at q<0.05 ==")
        for r in sorted([r for r in rows if r["unit"] == unit], key=lambda r: -r["observed"]):
            tags = []
            for name in ("abundance", "degree"):
                if r[f"q_enrich_{name}"] < 0.05:
                    tags.append(f"ENRICHED vs {name} (O/E {r[f'oe_{name}']:.1f})")
                if r[f"q_deplete_{name}"] < 0.05:
                    tags.append(f"depleted vs {name} (O/E {r[f'oe_{name}']:.2f})")
            if tags:
                print(f"  {r['pair']:>6}: obs {r['observed']:>4} | " + "; ".join(tags))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
