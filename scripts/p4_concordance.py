#!/usr/bin/env python3
"""Phase 4 — node-level cross-layer concordance (primary analysis).

For the primary full-5 node set, each element i has a domain peptide in all five
layers. We ask whether i's nearest neighbours are the SAME elements in every layer.

Per ordered layer pair (A,B) and node i:
  * jaccard_i(A,B) = |kNN_A(i) ∩ kNN_B(i)| / |kNN_A(i) ∪ kNN_B(i)|   (threshold-free)
  * rho_i(A,B)     = Spearman of similarities over that union

A high jaccard = the family holds across those two domains. Chimera candidates are
elements where ONE domain is discordant with the others while those others agree
(majority vote). "Concordant beyond chance" is calibrated with a degree-preserving
permutation null, and candidate calls are reported at an empirical FDR.

Missing domain = missing data (excluded from that pair, tallied separately). Here,
on the full-5 primary set, all five are present by construction; a layer with no
surviving edge for i simply gives an empty neighbour set.

Usage:
    p4_concordance.py --edges-dir DIR --layers GAG PROT RT RH INT \
        --nodes analysis_nodes.tsv --out-dir OUT \
        [--k 20] [--perms 200] [--null-quantile 0.99] [--seed 0]
"""

from __future__ import annotations

import argparse
import itertools
import os
import random
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr


def read_nodes(path):
    nodes = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            s = line.strip()
            if s:
                nodes.append(s.split("\t")[0])
    return nodes


def read_layer_edges(path, idx):
    """Return symmetric similarity dict {i:{j:ratio}} restricted to known nodes."""
    sim = defaultdict(dict)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b, ratio = f[0], f[1], float(f[4])
            if a in idx and b in idx and a != b:
                if ratio > sim[a].get(b, 0):
                    sim[a][b] = ratio
                    sim[b][a] = ratio
    return sim


def knn_sets(sim, nodes, k):
    """Top-k neighbour set per node (by ratio)."""
    out = {}
    for n in nodes:
        nb = sim.get(n, {})
        if not nb:
            out[n] = set()
            continue
        top = sorted(nb.items(), key=lambda kv: kv[1], reverse=True)[:k]
        out[n] = {x for x, _ in top}
    return out


def jaccard(a, b):
    if not a and not b:
        return np.nan            # no information in either layer
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def pair_rho(simA, simB, i, knnA, knnB):
    union = knnA[i] | knnB[i]
    if len(union) < 3:
        return np.nan
    va = [simA.get(i, {}).get(j, 0.0) for j in union]
    vb = [simB.get(i, {}).get(j, 0.0) for j in union]
    r = spearmanr(va, vb).correlation
    return r


def degree_preserving_perm(knn, nodes, rng):
    """Permute node identities within degree bins (preserves degree sequence)."""
    by_deg = defaultdict(list)
    for n in nodes:
        by_deg[len(knn[n])].append(n)
    perm = {}
    for deg, group in by_deg.items():
        shuffled = group[:]
        rng.shuffle(shuffled)
        for a, b in zip(group, shuffled):
            perm[a] = b
    return perm


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges-dir", required=True)
    ap.add_argument("--layers", nargs="+", default=["GAG", "PROT", "RT", "RH", "INT"])
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--edge-suffix", default=".edges.tsv")
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--null-quantile", type=float, default=0.99)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "phase4")
    os.makedirs(out, exist_ok=True)
    rng = random.Random(args.seed)

    nodes = read_nodes(args.nodes)
    idx = set(nodes)
    print(f"[p4] primary nodes: {len(nodes)}; layers: {args.layers}; k={args.k}")

    sim = {lyr: read_layer_edges(
        os.path.join(args.edges_dir, lyr + args.edge_suffix), idx)
        for lyr in args.layers}
    knn = {lyr: knn_sets(sim[lyr], nodes, args.k) for lyr in args.layers}
    pairs = list(itertools.combinations(args.layers, 2))

    # ---- observed per-node, per-pair jaccard + rho ------------------------ #
    jac = {p: {} for p in pairs}
    rho = {p: {} for p in pairs}
    for (A, B) in pairs:
        for i in nodes:
            jac[(A, B)][i] = jaccard(knn[A][i], knn[B][i])
            rho[(A, B)][i] = pair_rho(sim[A], sim[B], i, knn[A], knn[B])

    scores_path = os.path.join(out, "concordance_scores.tsv")
    with open(scores_path, "w") as f:
        f.write("ns_id\t" + "\t".join(f"jac_{A}_{B}" for A, B in pairs)
                + "\t" + "\t".join(f"rho_{A}_{B}" for A, B in pairs)
                + "\tmean_jac\n")
        for i in nodes:
            jvals = [jac[p][i] for p in pairs]
            rvals = [rho[p][i] for p in pairs]
            mj = np.nanmean(jvals) if not all(np.isnan(jvals)) else np.nan
            f.write(i + "\t"
                    + "\t".join(_fmt(v) for v in jvals) + "\t"
                    + "\t".join(_fmt(v) for v in rvals) + "\t"
                    + _fmt(mj) + "\n")

    # ---- null: degree-preserving permutation per pair -> tau -------------- #
    tau = {}
    for (A, B) in pairs:
        null_vals = []
        for _ in range(args.perms):
            perm = degree_preserving_perm(knn[A], nodes, rng)
            for i in nodes:
                v = jaccard(knn[A][perm[i]], knn[B][i])
                if not np.isnan(v):
                    null_vals.append(v)
        tau[(A, B)] = (float(np.quantile(null_vals, args.null_quantile))
                       if null_vals else 1.0)

    # ---- majority-vote outlier localisation ------------------------------- #
    def call_outliers(jac_table):
        """Return {node: outlier_domain} where 4 domains form a concordant clique
        and exactly one is disconnected from all of them."""
        calls = {}
        for i in nodes:
            concordant = {A: set() for A in args.layers}
            for (A, B) in pairs:
                v = jac_table[(A, B)][i]
                if not np.isnan(v) and v >= tau[(A, B)]:
                    concordant[A].add(B)
                    concordant[B].add(A)
            others = len(args.layers) - 1
            # a domain is "core" if concordant with all other domains that it can be
            core = [d for d in args.layers if len(concordant[d]) == others]
            outliers = [d for d in args.layers if len(concordant[d]) == 0]
            if len(core) == others and len(outliers) == 1:
                calls[i] = outliers[0]
        return calls

    observed = call_outliers(jac)

    # empirical FDR: rerun localisation on permuted layers
    null_calls = 0
    for _ in range(max(1, args.perms // 10)):
        permuted = {}
        for (A, B) in pairs:
            perm = degree_preserving_perm(knn[A], nodes, rng)
            permuted[(A, B)] = {i: jaccard(knn[A][perm[i]], knn[B][i])
                                for i in nodes}
        null_calls += len(call_outliers(permuted))
    exp_null = null_calls / max(1, args.perms // 10)
    fdr = exp_null / len(observed) if observed else 0.0

    cand_path = os.path.join(out, "candidates.tsv")
    with open(cand_path, "w") as f:
        f.write("ns_id\toutlier_domain\tmean_jac\n")
        for i, d in sorted(observed.items()):
            mj = np.nanmean([jac[p][i] for p in pairs])
            f.write(f"{i}\t{d}\t{_fmt(mj)}\n")

    summary = os.path.join(out, "phase4_summary.md")
    with open(summary, "w") as f:
        f.write("# Phase 4 — cross-layer concordance\n\n")
        f.write(f"- primary nodes: {len(nodes)}; k={args.k}; perms={args.perms}\n")
        f.write(f"- candidate single-domain-outlier elements: "
                f"**{len(observed)}**\n")
        f.write(f"- expected under null: {exp_null:.1f}  -> empirical FDR "
                f"~ **{fdr:.2f}**\n\n")
        f.write("## Per-pair concordance threshold (tau, from null)\n\n")
        f.write("| pair | tau | median observed jac |\n|---|---:|---:|\n")
        for (A, B) in pairs:
            med = np.nanmedian([jac[(A, B)][i] for i in nodes])
            f.write(f"| {A}-{B} | {tau[(A,B)]:.3f} | {_fmt(med)} |\n")
        f.write("\n## Outlier-domain tally among candidates\n\n")
        tally = defaultdict(int)
        for d in observed.values():
            tally[d] += 1
        for d in args.layers:
            f.write(f"- {d}: {tally[d]}\n")
        f.write("\n> Discordance != recombination. These are candidates for Phase 6 "
                "artefact filtering and Phase 7 phylogenetic confirmation.\n")

    print(f"[p4] candidates: {len(observed)} (empirical FDR ~ {fdr:.2f})")
    print(f"[p4] wrote {scores_path}, {cand_path}, {summary}")


def _fmt(v):
    return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.4f}"


if __name__ == "__main__":
    main()
