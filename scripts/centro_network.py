#!/usr/bin/env python3
"""ATHILA family recombination network coloured by centrophilic/centrophobic class.

Nodes  = families (area ~ number of elements), coloured red (centrophilic),
         gray (neutral) or blue (centrophobic); arranged by class.
Edges  = family pairs with recombinant mosaics at >=1 distinct locus; width ~ loci.
         Dark = more than expected under the degree-preserving null (q<0.05);
         light = not significant. Dashed = significantly fewer than expected.

Usage: centro_network.py --classes family_centro_class.tsv \
           --perm pair_permutation.tsv --out network.png
"""
import argparse
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COL = {"centrophilic": "#e34948", "neutral": "#9a9890", "centrophobic": "#2a78d6",
       "untested": "#9a9890"}
INK, MUTED, SURF = "#0b0b0b", "#52514e", "#fcfcfb"

ap = argparse.ArgumentParser()
ap.add_argument("--classes", required=True)
ap.add_argument("--perm", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

fam = {}
for l in open(args.classes).readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    fam[f[0]] = {"class": f[1], "n": int(f[2]), "pc": float(f[3])}

edges = []
for l in open(args.perm).readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    if f[0] != "loci":
        continue
    a, b = f[1].split("-")
    if a not in fam or b not in fam:
        continue
    obs, oe_d = int(f[2]), float(f[8])
    q_en_d, q_de_a, q_de_d = float(f[9]), float(f[6]), float(f[10])
    edges.append((a, b, obs, oe_d, q_en_d < 0.05, min(q_de_a, q_de_d) < 0.05))

# layout: classes on arcs of one circle
order = ["centrophilic", "neutral", "centrophobic"]
PREF = ["5", "1", "6b", "6a", "9", "0", "7a", "2", "4c", "6", "3", "4", "4a"]
ranked = sorted(fam, key=lambda x: (order.index(fam[x]["class"]) if fam[x]["class"] in order else 1,
                                    PREF.index(x) if x in PREF else 99, -fam[x]["pc"]))
pos = {}
for i, x in enumerate(ranked):
    t = math.pi / 2 + 2 * math.pi * i / len(ranked)
    pos[x] = (math.cos(t), math.sin(t))

fig, ax = plt.subplots(figsize=(9, 9), facecolor=SURF)
ax.set_facecolor(SURF)
for a, b, obs, oe, enr, dep in sorted(edges, key=lambda e: e[2]):
    (x1, y1), (x2, y2) = pos[a], pos[b]
    if obs > 0:
        ax.plot([x1, x2], [y1, y2], color=INK if enr else "#c9c8c1",
                lw=1 + 2.2 * math.sqrt(obs), solid_capstyle="round", zorder=1,
                ls=(0, (2, 1.2)) if (dep and not enr) else "-")
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        L_ = math.hypot(x2 - x1, y2 - y1) or 1
        nx_, ny_ = -(y2 - y1) / L_, (x2 - x1) / L_
        if nx_ * mx + ny_ * my < 0:          # push label outward, away from centre
            nx_, ny_ = -nx_, -ny_
        mx, my = mx + 0.13 * nx_, my + 0.13 * ny_
        if enr:
            ax.text(mx, my, f"{obs} loci\nO/E {oe:.1f}", fontsize=7.5, color=INK,
                    ha="center", va="center", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none"))
    elif dep:
        ax.plot([x1, x2], [y1, y2], color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
for x in ranked:
    px, py = pos[x]
    r = 0.045 + 0.11 * math.sqrt(fam[x]["n"] / max(v["n"] for v in fam.values()))
    ax.add_patch(plt.Circle((px, py), r, color=COL.get(fam[x]["class"], "#9a9890"),
                            ec=SURF, lw=2, zorder=3))
    ax.text(px, py, x, ha="center", va="center", fontsize=10, color="white",
            fontweight="bold", zorder=5)
    ax.text(px * 1.2, py * 1.2, f"{100*fam[x]['pc']:.0f}% cen", ha="center",
            va="center", fontsize=7.5, color=MUTED)
handles = [Line2D([0], [0], marker="o", ls="", ms=11, mfc=COL[c], mec=SURF, label=c)
           for c in ("centrophilic", "neutral", "centrophobic")]
handles += [Line2D([0], [0], color=INK, lw=3, label="recombination, enriched (q<0.05)"),
            Line2D([0], [0], color="#c9c8c1", lw=3, label="recombination, as expected"),
            Line2D([0], [0], color="#c9c8c1", lw=3, ls=(0, (2, 1.2)),
                   label="recombination, fewer than expected"),
            Line2D([0], [0], color=MUTED, lw=1, ls=(0, (3, 3)),
                   label="no recombination, significantly depleted")]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=2,
          frameon=False, fontsize=8.5)
ax.set_title("ATHILA family recombination network\n"
             "node = family (area ~ copies, % = share inside CEN178 core); "
             "edge = distinct recombinant loci", loc="left", fontsize=11, color=INK)
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.35); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor=SURF, bbox_inches="tight")
print("wrote", args.out)
