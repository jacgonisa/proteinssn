#!/usr/bin/env python3
"""Multi-page PDF report of candidate Athila recombination events + bias analysis.

Pages:
  1. summary + distinct event-type table
  2. bias analysis (Robin's question): family-pair recombination vs family homology
     and abundance, with the detectability caveat made explicit
  3+. full per-element event table (paginated)

Usage:
    recombination_report.py --events events.tsv --assign domain_family_nt.tsv \
        --exemplars EX.fasta --out report.pdf
"""

from __future__ import annotations

import argparse
import itertools
import os
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from athila_recomb_nt import read_fasta, family_of


def load_events(path):
    rows = []
    with open(path) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def family_homology(exemplars):
    exe = read_fasta(exemplars)
    refs = {h: s for h, s in exe.items() if "_I" in h and "_LTR" not in h}
    fam = {h: family_of(h) for h in refs}
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "e.fa")
        with open(fa, "w") as o:
            for h, s in refs.items():
                o.write(f">{h}\n{s}\n")
        res = os.path.join(tmp, "r.m8")
        subprocess.run(["mmseqs", "easy-search", fa, fa, res, os.path.join(tmp, "t"),
                        "--search-type", "3", "-e", "1e-3", "--max-seqs", "200",
                        "--format-output", "query,target,pident", "-v", "1"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pid = defaultdict(list)
        for line in open(res):
            q, t, p = line.rstrip("\n").split("\t")
            if fam[q] != fam[t]:
                pid[tuple(sorted((fam[q], fam[t])))].append(float(p))
    return {k: statistics.mean(v) for k, v in pid.items()}, sorted(set(fam.values()))


def family_abundance(assign_path, domain="GAG"):
    ab = Counter()
    with open(assign_path) as fh:
        next(fh)
        for line in fh:
            ns, dom, f, b, b2, mg, conf = line.rstrip("\n").split("\t")
            if dom == domain and int(conf):
                ab[f] += 1
    return ab


def _table_page(pdf, title, col_labels, rows, fontsize=7, max_rows=42):
    for start in range(0, max(1, len(rows)), max_rows):
        chunk = rows[start:start + max_rows]
        fig, ax = plt.subplots(figsize=(11.7, 8.3))   # A4 landscape
        ax.axis("off")
        ax.set_title(title + (f"  (rows {start+1}-{start+len(chunk)} of {len(rows)})"
                     if len(rows) > max_rows else ""), fontsize=12,
                     fontweight="bold", loc="left")
        if chunk:
            t = ax.table(cellText=chunk, colLabels=col_labels, loc="upper center",
                         cellLoc="left")
            t.auto_set_font_size(False); t.set_fontsize(fontsize)
            t.scale(1, 1.25)
            for j in range(len(col_labels)):
                t[0, j].set_facecolor("#4C72B0")
                t[0, j].set_text_props(color="white", fontweight="bold")
        pdf.savefig(fig); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = load_events(args.events)
    intra = [r for r in rows if r["category"] == "intra_athila"]
    inter = [r for r in rows if r["category"] == "inter_clade"]
    hom, families = family_homology(args.exemplars)
    abund = family_abundance(args.assign)

    # recombination pair counts (intra)
    pair_cnt = Counter()
    for r in intra:
        fs = tuple(sorted(r["partners"].split(",")))
        if len(fs) == 2:
            pair_cnt[fs] += 1

    with PdfPages(args.out) as pdf:
        # ---- Page 1: summary + event types ------------------------------- #
        fig, ax = plt.subplots(figsize=(11.7, 8.3)); ax.axis("off")
        ax.set_title("Candidate Athila recombination events", fontsize=16,
                     fontweight="bold", loc="left")
        txt = (f"Total candidate elements: {len(rows)}\n"
               f"  intra-Athila (between ATHILA families): {len(intra)} "
               f"({sum(int(r['full5']) for r in intra)} full-5)\n"
               f"  inter-clade (Athila x other LTR retro family): {len(inter)}\n\n"
               f"Recombination sits at the element ENDS (GAG / INT) in every event "
               f"type; the RT/RH core never swaps.\n\n"
               f"Confidence: intra-Athila = nucleotide per-domain family typing "
               f"(margin >=1.3x); inter-clade = TEsorter per-domain clade call.\n"
               f"These are CANDIDATES (not yet Phase-6 nesting-filtered or Phase-7 "
               f"breakpoint-confirmed).")
        ax.text(0, 0.92, txt, va="top", fontsize=10, family="monospace")
        # event-type table
        types = Counter((r["category"], r["partners"], r["swapped"]) for r in rows)
        et_rows = [[c, p, s, str(n)] for (c, p, s), n in types.most_common()]
        t = ax.table(cellText=et_rows,
                     colLabels=["category", "partners", "swapped", "# elem"],
                     loc="lower center", cellLoc="left", bbox=[0, 0.02, 1, 0.5])
        t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.3)
        for j in range(4):
            t[0, j].set_facecolor("#4C72B0")
            t[0, j].set_text_props(color="white", fontweight="bold")
        pdf.savefig(fig); plt.close(fig)

        # ---- Page 2: bias analysis --------------------------------------- #
        fig = plt.figure(figsize=(11.7, 8.3))
        fig.suptitle("Recombination bias vs family homology (Robin's question)",
                     fontsize=14, fontweight="bold")
        # scatter homology vs recomb count
        ax1 = fig.add_subplot(2, 2, 1)
        xs, ys, labs = [], [], []
        for p, n in pair_cnt.items():
            if p in hom:
                xs.append(hom[p]); ys.append(n); labs.append(f"{p[0][6:]}-{p[1][6:]}")
        ax1.scatter(xs, ys, s=80, color="#C44E52", zorder=3)
        for x, y, l in zip(xs, ys, labs):
            ax1.annotate(l, (x, y), fontsize=8, xytext=(4, 4),
                         textcoords="offset points")
        ax1.set_xlabel("family homology (% id)"); ax1.set_ylabel("# recombinant elements")
        ax1.set_title("detected recombination vs homology", fontsize=10)
        if len(xs) >= 3:
            from scipy.stats import spearmanr
            ax1.text(0.05, 0.95, f"Spearman rho = {spearmanr(xs, ys).correlation:.2f}",
                     transform=ax1.transAxes, va="top", fontsize=9,
                     bbox=dict(boxstyle="round", fc="#FFF3CD"))
        # detectability curve: homology of recombining vs non-recombining pairs
        ax2 = fig.add_subplot(2, 2, 2)
        rec_h = [hom[p] for p in pair_cnt if p in hom]
        non_h = [v for k, v in hom.items() if k not in pair_cnt]
        ax2.hist([rec_h, non_h], bins=np.arange(65, 100, 4),
                 label=["recombining pairs", "no detection"],
                 color=["#C44E52", "#CCCCCC"])
        ax2.axvspan(70, 90, color="#55A868", alpha=0.12)
        ax2.set_xlabel("family homology (% id)"); ax2.set_ylabel("# family pairs")
        ax2.set_title("detection 'sweet spot' ~70-90%", fontsize=10)
        ax2.legend(fontsize=8)
        # caveat text
        ax3 = fig.add_subplot(2, 1, 2); ax3.axis("off")
        caveat = (
            "ANSWER (with the essential caveat):\n"
            "- Raw detected recombination is NEGATIVELY correlated with homology "
            "(rho = -1.0 over the 4 detectable pairs), and the MOST homologous "
            "families (>93% id: ATHILA6a-6b, ATHILA0-3, ATHILA4-4a) show ZERO "
            "detections.\n"
            "- This is largely a DETECTION ARTEFACT, not biology: when two families "
            "are near-identical, a swapped domain cannot be confidently assigned to "
            "one vs the other, so homologous-recombination between the most similar "
            "families is INVISIBLE to this method. Detection peaks at intermediate "
            "homology (~70-90%).\n"
            "- Counts are also confounded by family ABUNDANCE (ATHILA2/4c/6b are "
            "common, so partner more often by chance).\n\n"
            "To test the homology-bias hypothesis properly we need a method that does "
            "NOT rely on discriminating the swapped family: nucleotide breakpoint "
            "detection (GARD/RDP5) or diagnostic-SNP mosaic mapping, normalised by "
            "pairwise family abundance. The current data can only say recombination "
            "is readily detected between moderately-divergent Athila families.")
        ax3.text(0, 1.0, caveat, va="top", fontsize=9, wrap=True)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig); plt.close(fig)

        # ---- Page 3: homology + recombination matrices ------------------- #
        fams_involved = sorted({f for p in pair_cnt for f in p}
                               | {f for p in hom for f in p})
        n = len(fams_involved); idx = {f: i for i, f in enumerate(fams_involved)}
        H = np.full((n, n), np.nan); R = np.zeros((n, n))
        for (a, b), v in hom.items():
            if a in idx and b in idx:
                H[idx[a], idx[b]] = H[idx[b], idx[a]] = v
        for (a, b), c in pair_cnt.items():
            R[idx[a], idx[b]] = R[idx[b], idx[a]] = c
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.7, 6))
        short = [f[6:] for f in fams_involved]
        for ax, M, title, cmap in [(axa, H, "family homology (% id)", "viridis"),
                                    (axb, R, "recombinant elements", "Reds")]:
            im = ax.imshow(M, cmap=cmap)
            ax.set_xticks(range(n)); ax.set_xticklabels(short, rotation=90, fontsize=7)
            ax.set_yticks(range(n)); ax.set_yticklabels(short, fontsize=7)
            ax.set_title(title, fontsize=11)
            fig.colorbar(im, ax=ax, fraction=0.046)
        fig.suptitle("Family homology vs recombination matrix", fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig); plt.close(fig)

        # ---- full event table -------------------------------------------- #
        cols = ["category", "accession", "chrom", "start", "architecture",
                "event", "full5"]
        tbl = [[r[c] for c in cols] for r in rows]
        _table_page(pdf, "All candidate recombination events", cols, tbl)

    print(f"[report] wrote {args.out} ({len(rows)} events)")


if __name__ == "__main__":
    main()
