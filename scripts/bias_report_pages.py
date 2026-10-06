#!/usr/bin/env python3
"""PDF pages for the k-mer homology-bias analysis: figure, stats, pair table."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

png, pairs_tsv, stats_txt, out = sys.argv[1:5]
ink, muted = "#0b0b0b", "#52514e"
with PdfPages(out) as pdf:
    fig, ax = plt.subplots(figsize=(11.7, 8.3)); ax.axis("off")
    ax.set_title("Robin's question — recombination vs family homology (k-mer method)",
                 loc="left", fontsize=15, fontweight="bold", color=ink)
    ax.text(0, 0.95, open(stats_txt).read(), va="top", fontsize=9.5,
            family="monospace", color=ink, wrap=True)
    pdf.savefig(fig); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.7, 8.3)); ax.axis("off")
    ax.imshow(mpimg.imread(png)); pdf.savefig(fig); plt.close(fig)

    rows = [l.rstrip("\n").split("\t") for l in open(pairs_tsv)]
    hdr, body = rows[0], [r for r in rows[1:] if int(r[2]) > 0]
    keep = ["pair", "homology", "raw", "loci", "n_a", "n_b", "per_opp", "per_opp_det"]
    idx = [hdr.index(k) for k in keep]
    fig, ax = plt.subplots(figsize=(11.7, 8.3)); ax.axis("off")
    ax.set_title("Family pairs with >=1 recombinant mosaic (76 other pairs: zero)",
                 loc="left", fontsize=12, fontweight="bold", color=ink)
    t = ax.table(cellText=[[r[i] for i in idx] for r in body], colLabels=keep,
                 loc="upper center", cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.4)
    for j in range(len(keep)):
        t[0, j].set_facecolor("#2a78d6"); t[0, j].set_text_props(color="white")
    pdf.savefig(fig); plt.close(fig)
