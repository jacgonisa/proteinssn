#!/usr/bin/env python3
"""Recombination-bias analysis with multiple normalisations (Robin's question).

For every ATHILA family pair we compute:
  * detected recombinant elements (from the master list),
  * family homology (mean % identity of aligned consensus),
  * DIAGNOSTIC SITES (fixed differences between the two family consensus) =
    a direct measure of DETECTABILITY (few sites => we can't see swaps),
  * family ABUNDANCE (elements per family) and OPPORTUNITY (n_A * n_B).

Then several normalisations of the detected count and their correlation with
homology, so we can separate true bias from detection/abundance artefacts.

Usage:
    recomb_bias.py --events events.tsv --assign domain_family_nt.tsv \
        --exemplars EX.fasta --out-dir OUT --mafft ~/bin/mafft
"""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import tempfile
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import spearmanr

from athila_recomb_nt import read_fasta, family_of


def family_reps(exemplars):
    """Longest internal (_I) exemplar per family."""
    best = {}
    for h, s in read_fasta(exemplars).items():
        if "_I" in h and "_LTR" not in h:
            fam = family_of(h)
            if fam not in best or len(s) > len(best[fam][1]):
                best[fam] = (h, s)
    return {f: s for f, (h, s) in best.items()}


def align(seqs: dict, mafft):
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "in.fa")
        with open(fa, "w") as o:
            for k, s in seqs.items():
                o.write(f">{k}\n{s}\n")
        out = subprocess.run([mafft, "--auto", "--quiet", fa],
                             capture_output=True, text=True, check=True).stdout
    aln, name, buf = {}, None, []
    for line in out.splitlines():
        if line.startswith(">"):
            if name:
                aln[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        aln[name] = "".join(buf)
    return aln


def pair_stats(aln):
    """For each family pair: homology (% id over aligned non-gap) + diagnostic sites."""
    fams = list(aln)
    hom, diag = {}, {}
    for a, b in itertools.combinations(fams, 2):
        sa, sb = aln[a], aln[b]
        same = diff = 0
        for ca, cb in zip(sa, sb):
            if ca in "-" or cb in "-" or ca.upper() == "N" or cb.upper() == "N":
                continue
            if ca.upper() == cb.upper():
                same += 1
            else:
                diff += 1
        tot = same + diff
        hom[(a, b)] = 100.0 * same / tot if tot else float("nan")
        diag[(a, b)] = diff            # fixed differences = diagnostic sites
    return hom, diag


def abundance(assign_path, domain):
    ab = Counter()
    with open(assign_path) as fh:
        next(fh)
        for line in fh:
            ns, dom, f, b, b2, mg, conf = line.rstrip("\n").split("\t")
            if dom == domain and int(conf):
                ab[f] += 1
    return ab


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--mafft", default=os.path.expanduser("~/bin/mafft"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # detected recombination per family pair (intra-Athila, element-weighted)
    det = Counter()
    with open(args.events) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "intra_athila":
                fs = tuple(sorted(f[9].split(",")))
                if len(fs) == 2:
                    det[fs] += 1

    reps = family_reps(args.exemplars)
    aln = align(reps, args.mafft)
    hom, diag = pair_stats(aln)
    ab_rt = abundance(args.assign, "RT")
    ab_gag = abundance(args.assign, "GAG")
    print(f"[bias] families aligned: {len(reps)}; family pairs: {len(hom)}")

    # normalise det so the pair key matches hom/diag keys (both sorted tuples)
    def norm_key(p):
        return tuple(sorted(p))

    rows = []
    for pair in hom:
        a, b = pair
        d = det.get(norm_key(pair), 0)
        opp_rt = ab_rt.get(a, 0) * ab_rt.get(b, 0)
        opp_gag = ab_gag.get(a, 0) * ab_gag.get(b, 0)
        ds = diag[pair]
        rows.append({
            "pair": f"{a[6:]}-{b[6:]}", "a": a, "b": b, "detected": d,
            "homology": hom[pair], "diag_sites": ds,
            "opp_RT": opp_rt, "opp_GAG": opp_gag,
            # normalisations
            "per_opp_RT": d / opp_rt * 1e6 if opp_rt else float("nan"),
            "per_opp_GAG": d / opp_gag * 1e6 if opp_gag else float("nan"),
            "per_diag": d / ds if ds else float("nan"),
            "per_opp_diag": (d / (opp_gag * ds) * 1e9
                             if opp_gag and ds else float("nan")),
        })

    tsv = os.path.join(args.out_dir, "bias_normalizations.tsv")
    cols = ["pair", "detected", "homology", "diag_sites", "opp_GAG", "opp_RT",
            "per_opp_GAG", "per_opp_RT", "per_diag", "per_opp_diag"]
    with open(tsv, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: -r["detected"]):
            out.write("\t".join(
                (f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")

    # correlations with homology
    print("\n[bias] Spearman correlation of each metric with family homology")
    print("       (over ALL pairs, and over pairs with >=1 detection):")
    detected_pairs = [r for r in rows if r["detected"] > 0]
    for metric in ["detected", "per_opp_GAG", "per_opp_RT", "per_diag",
                   "per_opp_diag", "diag_sites"]:
        for label, subset in [("all", rows), ("detected-only", detected_pairs)]:
            xs = [r["homology"] for r in subset if np.isfinite(r[metric])]
            ys = [r[metric] for r in subset if np.isfinite(r[metric])]
            if len(xs) >= 3:
                rho = spearmanr(xs, ys).correlation
                print(f"   {metric:>14} vs homology [{label:>13}]: "
                      f"rho={rho:+.2f} (n={len(xs)})")

    _plots(rows, os.path.join(args.out_dir, "bias_normalizations.png"))
    print(f"\n[bias] wrote {tsv} and bias_normalizations.png")


def _plots(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    metrics = [("detected", "raw detected"),
               ("per_opp_GAG", "per opportunity (abundance, GAG)"),
               ("per_diag", "per diagnostic site (detectability)"),
               ("per_opp_diag", "per opportunity x diagnostic site")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (m, title) in zip(axes.flat, metrics):
        xs = [r["homology"] for r in rows if np.isfinite(r[m])]
        ys = [r[m] for r in rows if np.isfinite(r[m])]
        det = [r["detected"] for r in rows if np.isfinite(r[m])]
        labs = [r["pair"] for r in rows if np.isfinite(r[m])]
        ax.scatter(xs, ys, s=[20 + 6 * d for d in det],
                   c=["#C44E52" if d > 0 else "#CCCCCC" for d in det], zorder=3)
        for x, y, l, d in zip(xs, ys, labs, det):
            if d > 0:
                ax.annotate(l, (x, y), fontsize=7, xytext=(3, 3),
                            textcoords="offset points")
        ax.set_xlabel("family homology (% id)"); ax.set_ylabel(title)
        ax.set_title(title, fontsize=10)
        if len([1 for d in det if d > 0]) >= 3:
            dxs = [x for x, d in zip(xs, det) if d > 0]
            dys = [y for y, d in zip(ys, det) if d > 0]
            rho = spearmanr(dxs, dys).correlation
            ax.text(0.04, 0.95, f"rho(detected)={rho:+.2f}", transform=ax.transAxes,
                    va="top", fontsize=8, bbox=dict(boxstyle="round", fc="#FFF3CD"))
    fig.suptitle("Recombination vs homology under different normalisations\n"
                 "(red = pairs with detected recombination; grey = none)",
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
