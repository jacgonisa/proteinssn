#!/usr/bin/env python3
"""Does ATHILA recombination track family homology? (Robin's question)

Uses the robust k-mer mosaic calls (internal domain span) so similar families are
not invisible. For EVERY family pair (zeros included) computes the mosaic count and
several normalisations, then tests the association with family homology:

  raw            mosaic elements
  loci           distinct loci (same chrom, start within --locus-bp) -> collapses the
                 same ancestral insertion seen in many accessions
  per_opp        count / (n_A * n_B)                 abundance opportunity
  obs_exp        count / expected-under-abundance    (same ranking as per_opp)
  per_min        count / min(n_A, n_B)               per element of the rarer family
  per_opp_det    per_opp / min(markers_A, markers_B) detectability-adjusted

Tests: Spearman with homology + permutation p (shuffle homology over pairs), and a
negative-binomial rate model  count ~ homology + log(detectability) +
offset(log n_A n_B)  -> incidence-rate ratio per +10 % homology.

Usage:
    recomb_bias_kmer.py --recomb kmer_recombinants.tsv --primary primary_family.tsv \
        --markers robust_markers.tsv --homology bias_normalizations.tsv --out-dir OUT
"""

from __future__ import annotations

import argparse
import itertools
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import spearmanr


def short(f):
    return f[6:] if f.startswith("ATHILA") else f


def load(args):
    pair_cnt = Counter()
    loci = defaultdict(set)
    for line in open(args.recomb).readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        chrom, start, fams = f[2], f[3], f[7].split(",")
        try:
            lb = int(start) // args.locus_bp
        except ValueError:
            lb = None
        for a, b in itertools.combinations(sorted(fams), 2):
            pair_cnt[(a, b)] += 1
            loci[(a, b)].add((chrom, lb))
    ab = Counter(short(l.rstrip("\n").split("\t")[1])
                 for l in open(args.primary).readlines()[1:])
    mk = Counter(short(l.rstrip("\n").split("\t")[1])
                 for l in open(args.markers).readlines()[1:])
    hom = {}
    for l in open(args.homology).readlines()[1:]:
        c = l.rstrip("\n").split("\t")
        hom[tuple(sorted(c[0].split("-")))] = float(c[2])
    return pair_cnt, loci, ab, mk, hom


def perm_p(x, y, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    obs = spearmanr(x, y).correlation
    y = np.asarray(y)
    null = np.array([spearmanr(x, rng.permutation(y)).correlation for _ in range(n)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recomb", required=True)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--homology", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--locus-bp", type=int, default=500_000)
    ap.add_argument("--perms", type=int, default=10000)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    pair_cnt, loci, ab, mk, hom = load(args)
    tot_opp = sum(ab[a] * ab[b] for a, b in hom if ab[a] and ab[b])
    tot_cnt = sum(pair_cnt[p] for p in hom)

    rows = []
    for (a, b), h in hom.items():
        if not (ab[a] and ab[b]):
            continue
        n = pair_cnt.get((a, b), 0)
        opp = ab[a] * ab[b]
        det = min(mk[a], mk[b]) or 1
        rows.append({
            "pair": f"{a}-{b}", "homology": h, "raw": n,
            "loci": len(loci.get((a, b), ())), "n_a": ab[a], "n_b": ab[b],
            "markers_min": det,
            "per_opp": n / opp * 1e6,
            "obs_exp": n / (tot_cnt * opp / tot_opp) if tot_cnt else 0.0,
            "per_min": n / min(ab[a], ab[b]),
            "per_opp_det": n / opp / det * 1e9,
        })

    metrics = ["raw", "loci", "per_opp", "obs_exp", "per_min", "per_opp_det"]
    H = [r["homology"] for r in rows]
    results = []
    for m in metrics:
        rho, p = perm_p(H, [r[m] for r in rows], args.perms)
        results.append((m, rho, p))

    # negative-binomial rate model (raw counts and distinct loci)
    glm_lines = []
    try:
        import statsmodels.api as sm
        X = np.column_stack([
            np.array(H) / 10.0,
            np.log([r["markers_min"] for r in rows])])
        X = sm.add_constant(X)
        off = np.log([r["n_a"] * r["n_b"] for r in rows])
        for resp in ["raw", "loci"]:
            y = np.array([r[resp] for r in rows])
            for name, fam in [("NegBin", sm.families.NegativeBinomial(alpha=1.0)),
                              ("Poisson", sm.families.Poisson())]:
                fit = sm.GLM(y, X, family=fam, offset=off).fit()
                b, se = fit.params[1], fit.bse[1]
                glm_lines.append(
                    f"{resp:>5} {name:>7}: IRR per +10% homology = {np.exp(b):.2f} "
                    f"(95% CI {np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}), "
                    f"p = {fit.pvalues[1]:.2g}")
    except Exception as e:
        glm_lines.append(f"GLM failed: {e}")

    tsv = os.path.join(args.out_dir, "bias_kmer_pairs.tsv")
    cols = ["pair", "homology", "raw", "loci", "n_a", "n_b", "markers_min"] + metrics[2:]
    with open(tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: -r["raw"]):
            f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")

    md = os.path.join(args.out_dir, "bias_kmer_summary.md")
    with open(md, "w") as f:
        f.write("# Recombination vs homology — robust k-mer mosaics (internal span)\n\n")
        f.write(f"- family pairs tested (both families present): {len(rows)}\n")
        f.write(f"- pairs with >=1 mosaic: {sum(r['raw'] > 0 for r in rows)}\n")
        f.write(f"- mosaic elements: {tot_cnt}\n\n")
        f.write("## Spearman with homology (all pairs, zeros included; permutation p)\n\n")
        f.write("| metric | rho | p |\n|---|---:|---:|\n")
        for m, rho, p in results:
            f.write(f"| {m} | {rho:+.2f} | {p:.3g} |\n")
        f.write("\n## Rate model: count ~ homology + log(detectability) + offset(log n_A n_B)\n\n")
        for l in glm_lines:
            f.write(f"- {l}\n")

    print(open(md).read())
    _plot(rows, results, os.path.join(args.out_dir, "bias_kmer.png"))
    print(f"[bias-kmer] wrote {tsv}, {md}, bias_kmer.png")


def _plot(rows, results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, accent, zero = "#0b0b0b", "#52514e", "#2a78d6", "#b8b7b0"
    panels = [("raw", "mosaic elements"), ("loci", "distinct loci"),
              ("per_opp", "per opportunity (n_A·n_B, ×1e6)"),
              ("per_opp_det", "per opportunity, detectability-adjusted")]
    rho = {m: (r, p) for m, r, p in results}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), facecolor="#fcfcfb")
    for ax, (m, title) in zip(axes.flat, panels):
        ax.set_facecolor("#fcfcfb")
        x = np.array([r["homology"] for r in rows])
        y = np.array([r[m] for r in rows])
        pos = y > 0
        ax.scatter(x[~pos], np.zeros((~pos).sum()), s=28, facecolors="none",
                   edgecolors=zero, linewidths=1.2, label="no mosaic detected")
        ax.scatter(x[pos], y[pos], s=40, color=accent, edgecolors="#fcfcfb",
                   linewidths=1.5, zorder=3, label="mosaics detected")
        top = np.argsort(-y)[:4]
        for i in top:
            if y[i] > 0:
                ax.annotate(rows[i]["pair"], (x[i], y[i]), fontsize=8, color=muted,
                            xytext=(4, 3), textcoords="offset points")
        ax.set_yscale("symlog", linthresh=max(1e-3, np.min(y[pos]) if pos.any() else 1))
        ax.set_xlabel("family homology (% identity, consensus alignment)", color=muted)
        ax.set_title(title, loc="left", fontsize=11, color=ink)
        r, p = rho[m]
        ax.text(0.02, 0.97, f"Spearman ρ = {r:+.2f}  (perm. p = {p:.2g})",
                transform=ax.transAxes, va="top", fontsize=9, color=ink)
        ax.grid(axis="y", color="#e6e5e0", linewidth=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#c9c8c1")
        ax.tick_params(colors=muted, labelsize=8)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 0.9))
    fig.suptitle("Does ATHILA recombination track family homology?  "
                 "(all family pairs, robust k-mer mosaics, internal span)",
                 fontsize=12, color=ink, x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=200, facecolor="#fcfcfb")
    plt.close(fig)


if __name__ == "__main__":
    main()
