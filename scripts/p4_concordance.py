#!/usr/bin/env python3
"""Phase 4 — node-level cross-layer concordance (primary analysis).

The full-5 SSNs are effectively COMPLETE graphs (every Athila element is alignable
to every other in every domain), so top-k neighbour-set overlap saturates and its
permutation null collapses onto the data. We therefore score concordance on the
FULL similarity profile:

  rho_i(A,B) = Spearman( profile of i to all others in layer A ,
                         profile of i to all others in layer B )

For element i and domain d, its coherence with the rest of the element is
  s_i(d) = mean_{d'!=d} rho_i(d, d')
A single-domain-outlier (candidate chimeric / recombinant) element has one domain
whose s_i(d) sits far below the element's other domains:
  o_i(d) = median_{e!=d} s_i(e) - s_i(d)      (large positive = outlier)

Significance is calibrated with a permutation null that breaks the cross-layer
identity link for one layer at a time (compare i's profile in a layer to a RANDOM
element's profile), giving a null distribution of o and an empirical FDR.

Missing domain = missing data: a layer with no edge for i gives a zero profile and
is handled as such (its rho is defined over the shared columns).

Usage:
    p4_concordance.py --edges-dir DIR --layers GAG PROT RT RH INT \
        --nodes analysis_nodes.tsv --out-dir OUT [--perms 500] [--fdr 0.05] [--seed 0]
"""

from __future__ import annotations

import argparse
import itertools
import os

import numpy as np


def read_nodes(path):
    nodes = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            s = line.strip()
            if s:
                nodes.append(s.split("\t")[0])
    return nodes


def load_matrix(path, idx):
    n = len(idx)
    m = np.zeros((n, n), dtype=np.float64)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b = f[0], f[1]
            if a in idx and b in idx:
                r = float(f[4])
                ia, ib = idx[a], idx[b]
                if r > m[ia, ib]:
                    m[ia, ib] = m[ib, ia] = r
    return m


def rank_rows(m):
    """Row-wise rank transform (for Spearman via Pearson on ranks)."""
    order = m.argsort(axis=1)
    ranks = np.empty_like(order, dtype=np.float64)
    cols = np.arange(m.shape[1])
    for i in range(m.shape[0]):
        ranks[i, order[i]] = cols
    return ranks


def rowwise_corr(X, Y):
    """Per-row Pearson correlation between two equal-shape matrices."""
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    num = (Xc * Yc).sum(axis=1)
    den = np.sqrt((Xc ** 2).sum(axis=1) * (Yc ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges-dir", required=True)
    ap.add_argument("--layers", nargs="+",
                    default=["GAG", "PROT", "RT", "RH", "INT"])
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--edge-suffix", default=".edges.tsv")
    ap.add_argument("--perms", type=int, default=500)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "phase4")
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    nodes = read_nodes(args.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    layers = args.layers
    print(f"[p4] primary nodes: {N}; layers: {layers}; perms={args.perms}")

    ranks = {L: rank_rows(load_matrix(
        os.path.join(args.edges_dir, L + args.edge_suffix), idx)) for L in layers}
    pairs = list(itertools.combinations(layers, 2))

    # ---- observed per-node pairwise profile rho --------------------------- #
    rho = {p: rowwise_corr(ranks[p[0]], ranks[p[1]]) for p in pairs}

    # per-domain coherence s_i(d) and outlier score o_i(d)
    s = {d: np.nanmean(np.vstack([rho[p] for p in pairs if d in p]), axis=0)
         for d in layers}
    S = np.vstack([s[d] for d in layers])                 # (n_layers, N)
    o = np.full_like(S, np.nan)
    for di, d in enumerate(layers):
        others = np.delete(S, di, axis=0)
        o[di] = np.nanmedian(others, axis=0) - S[di]
    best_dom = np.array(layers)[np.nanargmax(o, axis=0)]
    best_o = np.nanmax(o, axis=0)

    # ---- permutation null for o ------------------------------------------- #
    null_o = []
    for _ in range(args.perms):
        for di, d in enumerate(layers):
            perm = rng.permutation(N)
            rp = {}
            for p in pairs:
                if d not in p:
                    rp[p] = rho[p]
                else:
                    a, b = p
                    if a == d:
                        rp[p] = rowwise_corr(ranks[a][perm], ranks[b])
                    else:
                        rp[p] = rowwise_corr(ranks[a], ranks[b][perm])
            s_d = np.nanmean(np.vstack([rp[p] for p in pairs if d in p]), axis=0)
            s_others = np.vstack([np.nanmean(
                np.vstack([rp[p] for p in pairs if e in p]), axis=0)
                for e in layers if e != d])
            null_o.append(np.nanmedian(s_others, axis=0) - s_d)
    null_o = np.concatenate(null_o)
    null_o = null_o[~np.isnan(null_o)]

    # threshold at empirical FDR
    tau = float(np.quantile(null_o, 1 - args.fdr))
    n_perm_blocks = args.perms * len(layers)
    obs_mask = best_o > tau
    n_obs = int(obs_mask.sum())
    exp_null = (null_o > tau).sum() / n_perm_blocks   # per single-layer scan
    fdr = exp_null / n_obs if n_obs else 0.0

    # ---- write ------------------------------------------------------------ #
    with open(os.path.join(out, "concordance_scores.tsv"), "w") as f:
        f.write("ns_id\t" + "\t".join(f"rho_{a}_{b}" for a, b in pairs)
                + "\t" + "\t".join(f"s_{d}" for d in layers)
                + "\toutlier_domain\toutlier_o\n")
        for i, n in enumerate(nodes):
            f.write(n + "\t"
                    + "\t".join(_fmt(rho[p][i]) for p in pairs) + "\t"
                    + "\t".join(_fmt(s[d][i]) for d in layers) + "\t"
                    + f"{best_dom[i]}\t{_fmt(best_o[i])}\n")

    order = np.argsort(-best_o)
    with open(os.path.join(out, "candidates.tsv"), "w") as f:
        f.write("ns_id\toutlier_domain\toutlier_o\ts_outlier\tmedian_s_others\n")
        for i in order:
            if obs_mask[i]:
                di = layers.index(best_dom[i])
                others = np.delete(S[:, i], di)
                f.write(f"{nodes[i]}\t{best_dom[i]}\t{best_o[i]:.4f}\t"
                        f"{S[di, i]:.4f}\t{np.nanmedian(others):.4f}\n")

    from collections import Counter
    tally = Counter(best_dom[i] for i in range(N) if obs_mask[i])
    with open(os.path.join(out, "phase4_summary.md"), "w") as f:
        f.write("# Phase 4 — cross-layer concordance (full-profile)\n\n")
        f.write(f"- primary nodes: {N}; perms: {args.perms}; FDR target: {args.fdr}\n")
        f.write(f"- outlier score threshold tau (null q{1-args.fdr:.2f}): {tau:.3f}\n")
        f.write(f"- candidate single-domain-outlier elements: **{n_obs}** "
                f"(empirical FDR ~ **{fdr:.2f}**)\n\n")
        f.write("## Median per-pair profile concordance (rho)\n\n")
        f.write("| pair | median rho |\n|---|---:|\n")
        for p in pairs:
            f.write(f"| {p[0]}-{p[1]} | {np.nanmedian(rho[p]):.3f} |\n")
        f.write("\n## Per-domain coherence s (median)\n\n")
        f.write("| domain | median s | rank |\n|---|---:|---:|\n")
        meds = sorted(((d, float(np.nanmedian(s[d]))) for d in layers),
                      key=lambda kv: -kv[1])
        for rk, (d, v) in enumerate(meds, 1):
            f.write(f"| {d} | {v:.3f} | {rk} |\n")
        f.write("\n## Outlier-domain tally among candidates\n\n")
        for d in layers:
            f.write(f"- {d}: {tally.get(d, 0)}\n")
        f.write("\n> Discordance != recombination. Candidates feed Phase 6 "
                "artefact filtering and Phase 7 phylogenetic confirmation.\n")

    print(f"[p4] candidates: {n_obs} (tau={tau:.3f}, empirical FDR ~ {fdr:.2f})")
    print(f"[p4] per-domain median coherence: "
          + ", ".join(f"{d}={np.nanmedian(s[d]):.3f}" for d in layers))


def _fmt(v):
    return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.4f}"


if __name__ == "__main__":
    main()
