#!/usr/bin/env python3
"""Figures for the HTML report (solo LTRs, age, junctions, fragments)."""
import os
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = sys.argv[1]
OUT = os.path.join(R, "report_figs")
os.makedirs(OUT, exist_ok=True)
SURF, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
ACC = "#2a78d6"
CLS = {"centrophilic": "#e34948", "neutral": "#9a9890", "centrophobic": "#2a78d6"}
CTX = {"centromeric": "#e34948", "pericentromeric": "#9a9890", "arm": "#2a78d6"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c1", "axes.labelcolor": MUTED,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 10.5,
                     "axes.titlelocation": "left", "axes.titlecolor": INK})


def axstyle(ax, grid="y"):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.6)
        ax.set_axisbelow(True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=150, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)


fam_class = {l.split("\t")[0]: l.split("\t")[1]
             for l in open(f"{R}/centro/family_centro_class.tsv").readlines()[1:]}

# ---- solo LTRs ------------------------------------------------------------- #
S = [l.rstrip("\n").split("\t") for l in open(f"{R}/solo/solo_intact_elements.tsv").readlines()[1:]]
fam = defaultdict(Counter); ctx = defaultdict(Counter)
for acc, kind, f, c, _ in S:
    fam[f][kind] += 1; ctx[c][kind] += 1
fams = [f for f in fam if sum(fam[f].values()) >= 100]
fams.sort(key=lambda f: fam[f]["solo"] / sum(fam[f].values()))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [2.2, 1]},
                             facecolor=SURF)
y = [fam[f]["solo"] / sum(fam[f].values()) for f in fams]
LTRCLS = {"ATHILA5": "centrophilic", "ATHILA1": "centrophilic", "ATHILA2": "centrophobic",
          "ATHILA4": "centrophobic", "ATHILA4c": "centrophobic", "ATHILA4a": "centrophobic"}
a1.barh([f.replace("ATHILA", "") + ("" if f in LTRCLS else " (group)") for f in fams], y,
        color=[CLS[LTRCLS.get(f, "neutral")] for f in fams], height=0.62)
from matplotlib.patches import Patch
a1.legend(handles=[Patch(color=CLS["centrophilic"], label="centrophilic"),
                   Patch(color=CLS["neutral"], label="mixed / neutral group"),
                   Patch(color=CLS["centrophobic"], label="centrophobic")],
          frameon=False, fontsize=7.5, loc="lower right")
for i, (f, v) in enumerate(zip(fams, y)):
    a1.text(v + 0.01, i, f"{v:.2f}  (n={sum(fam[f].values())})", va="center", fontsize=7.5, color=MUTED)
a1.set_xlim(0, 0.95); a1.set_xlabel("solo LTRs / (solo + intact)")
a1.set_title("Solo-LTR fraction by LTR family group"); axstyle(a1, "x")
ks = ["centromeric", "pericentromeric", "arm"]
v = [ctx[k]["solo"] / max(1, sum(ctx[k].values())) for k in ks]
a2.bar(ks, v, color=[CTX[k] for k in ks], width=0.6)
for i, x in enumerate(v):
    a2.text(i, x + 0.01, f"{x:.2f}", ha="center", fontsize=8, color=INK)
a2.set_ylim(0, 0.6); a2.set_title("by context"); a2.tick_params(axis="x", labelsize=8)
axstyle(a2)
save(fig, "solo_ltr.png")

# ---- age ------------------------------------------------------------------- #
A = [l.rstrip("\n").split("\t") for l in open(f"{R}/age/element_age.tsv").readlines()[1:]]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), facecolor=SURF)
data = [[float(r[1]) for r in A if r[2] == k] for k in ks]
bp = a1.boxplot(data, vert=False, widths=0.55, showfliers=False, patch_artist=True,
                medianprops=dict(color=INK, lw=1.5))
for p, k in zip(bp["boxes"], ks):
    p.set_facecolor(CTX[k]); p.set_alpha(0.75); p.set_edgecolor(SURF)
a1.set_yticks([1, 2, 3], ks); a1.set_xlabel("5′–3′ LTR identity (%)  →  younger")
a1.set_title("Element age by context (intact elements)"); axstyle(a1, "x")
ages = np.array([float(r[1]) for r in A]); mos = np.array([r[4] == "1" for r in A])
qs = np.quantile(ages, [0, .2, .4, .6, .8, 1])
labels, pct = [], []
for lo, hi in zip(qs, qs[1:]):
    sel = (ages >= lo) & (ages <= hi)
    labels.append(f"{lo:.0f}–{hi:.0f}%"); pct.append(mos[sel].mean() * 100)
a2.plot(range(5), pct, color=ACC, lw=2, marker="o", ms=7, mec=SURF, mew=1.5)
for i, x in enumerate(pct):
    a2.text(i, x + 0.25, f"{x:.1f}%", ha="center", fontsize=8, color=INK)
a2.set_xticks(range(5), labels, fontsize=7.5); a2.set_ylim(0, 7)
a2.set_ylabel("% recombinant mosaics"); a2.set_xlabel("LTR identity quintile (old → young)")
a2.set_title("Recombinants by age: two regimes"); axstyle(a2)
a2.annotate("oldest: ATHILA2×4c, 2×6b\n(pericentromeric)", (0, pct[0]), (0.6, 6.1), fontsize=7, color=MUTED)
a2.annotate("youngest: ATHILA6a×6b\n(half centromeric)", (4, pct[4]), (2.7, 6.1), fontsize=7, color=MUTED)
save(fig, "age.png")

# ---- microhomology (unique junctions) -------------------------------------- #
fig, ax0 = plt.subplots(figsize=(6.2, 3.4), facecolor=SURF)
for ax, tag, title in zip([ax0], ("full",), ("Athilafinder intact elements",)):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{R}/junctions_{tag}/junctions.tsv")]
    h = rows[0]; J = [dict(zip(h, r)) for r in rows[1:] if r[h.index("class")] == "deletion_clean"]
    uniq = {}
    for j in J:
        k = (j["fam1"], round(float(j["rel_pos_in_consensus1"] or 0), 2), round(int(j["size"]), -1))
        if j["microhomology_exact"] != "":
            uniq.setdefault(k, int(j["microhomology_exact"]))
    mh = np.array(list(uniq.values()))
    nul = np.array([int(x) for x in open(f"{R}/junctions_{tag}/microhomology_null.txt") if x.strip()])
    xs = np.arange(0, 7)
    o = [np.mean(np.minimum(mh, 6) == x) * 100 for x in xs]
    n = [np.mean(np.minimum(nul, 6) == x) * 100 for x in xs]
    ax.bar(xs - 0.19, o, 0.36, color=ACC, label=f"observed (unique junctions, n={len(mh)})")
    ax.bar(xs + 0.19, n, 0.36, color="#c9c8c1", label="random breakpoints (null)")
    ax.set_xticks(xs, [str(x) for x in xs[:-1]] + ["≥6"]); ax.set_xlabel("microhomology at deletion junction (bp)")
    ax.set_ylabel("% of junctions"); ax.set_title(f"Clean deletions — {title}")
    ax.legend(frameon=False, fontsize=7.5); axstyle(ax)
save(fig, "microhomology.png")

# ---- junction classes per 100 copies --------------------------------------- #
def rates(tag):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{R}/junctions_{tag}/elements.tsv")]
    h = rows[0]; E = [dict(zip(h, r)) for r in rows[1:]]
    return {k[2:]: np.mean([int(e[k]) for e in E]) * 100 for k in h if k.startswith("n_") and k != "n_pieces"}
ri = rates("full")
keys = ["deletion_clean", "deletion_with_insert", "internal_family_switch",
        "internal_family_switch_gapped", "inversion", "duplication"]
names = ["clean deletion", "deletion + insert", "family switch (clean)",
         "family switch (gapped)", "inversion", "duplication"]
fig, ax = plt.subplots(figsize=(10, 3.2), facecolor=SURF)
x = np.arange(len(keys))
ax.bar(x, [ri[k] for k in keys], 0.55, color=ACC)
for i, k in enumerate(keys):
    ax.text(i, ri[k] + 0.12, f"{ri[k]:.2f}", ha="center", fontsize=7.5, color=INK)
ax.set_xticks(x, names, fontsize=8); ax.set_ylabel("junctions per 100 elements")
ax.set_title("Rearrangement junctions in Athilafinder intact elements (n = 20,252)")
axstyle(ax)
save(fig, "junction_rates.png")

# ---- internal completeness of Athilafinder intact elements ------------------ #
rows = [l.rstrip("\n").split("\t") for l in open(f"{R}/junctions_full/elements.tsv")]
h = rows[0]; E = [dict(zip(h, r)) for r in rows[1:] if r[h.index("internal_family")]]
ctxmap = {l.split("\t")[0]: l.split("\t")[5] for l in open(f"{R}/centro/element_centro_class.tsv").readlines()[1:]}
byf = defaultdict(list); byc = defaultdict(list)
for e in E:
    c = float(e["internal_consensus_cov"]); f = e["internal_family"].replace("ATHILA", "")
    byf[f].append(c)
    if e["query"] in ctxmap:
        byc[ctxmap[e["query"]]].append(c)
fl = [f for f in byf if len(byf[f]) >= 100 and f in fam_class]
fl.sort(key=lambda f: np.mean(np.array(byf[f]) < 0.5))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [2.2, 1]},
                             facecolor=SURF)
vals = [np.mean(np.array(byf[f]) < 0.5) * 100 for f in fl]
a1.barh(fl, vals, color=[CLS.get(fam_class[f], "#9a9890") for f in fl], height=0.62)
for i, (f, v) in enumerate(zip(fl, vals)):
    a1.text(v + 0.6, i, f"{v:.0f}%  (n={len(byf[f])})", va="center", fontsize=7.5, color=MUTED)
a1.set_xlim(0, 48); a1.set_xlabel("% of intact elements missing > half of the internal consensus")
a1.set_title("Heavily deleted intact elements, by family")
from matplotlib.patches import Patch
a1.legend(handles=[Patch(color=c, label=k) for k, c in CLS.items()], frameon=False, fontsize=7.5,
          loc="lower right"); axstyle(a1, "x")
v = [np.mean(np.array(byc[k]) < 0.5) * 100 for k in ks]
a2.bar(ks, v, color=[CTX[k] for k in ks], width=0.6)
for i, x in enumerate(v):
    a2.text(i, x + 0.4, f"{x:.0f}%", ha="center", fontsize=8, color=INK)
a2.set_title("by context"); a2.tick_params(axis="x", labelsize=8); axstyle(a2)
save(fig, "completeness.png")

# ---- where internal deletions fall along the consensus ---------------------- #
dm = defaultdict(dict)
for l in open(f"{R}/fragments/domain_map.tsv").readlines()[1:]:
    r, d, x = l.rstrip().split("\t"); dm[r][d] = float(x)
rows = [l.rstrip("\n").split("\t") for l in open(f"{R}/junctions_full/junctions.tsv")]
h = rows[0]; J = [dict(zip(h, r)) for r in rows[1:]]
dels = [j for j in J if j["class"] in ("deletion_clean", "deletion_with_insert") and j["rel_pos_in_consensus1"]]
n_all = len(dels)
uniq = {(j["fam1"], round(float(j["rel_pos_in_consensus1"]), 2), round(int(j["size"]), -2)): float(j["rel_pos_in_consensus1"])
        for j in dels}                     # one entry per distinct deletion (inheritance collapsed)
pos = np.array(list(uniq.values()))
fig, ax = plt.subplots(figsize=(10, 2.9), facecolor=SURF)
hist, edges = np.histogram(pos, bins=40, range=(0, 1))
ax.bar(edges[:-1], hist / hist.mean(), width=1 / 40, align="edge", color=ACC, edgecolor=SURF, lw=0.6)
ax.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)))
dpos = defaultdict(list)
for r, d in dm.items():
    if len(d) >= 5:
        for k, x in d.items():
            dpos[k].append(x)
for k in ("GAG", "PROT", "RT", "RH", "INT"):
    if dpos[k]:
        m = float(np.median(dpos[k]))
        ax.axvline(m, color=INK, lw=0.8)
        ax.text(m, ax.get_ylim()[1] * 0.98, k, ha="center", va="top", fontsize=8, color=INK,
                bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none"))
ax.set_xlim(0, 1); ax.set_xlabel("position along family internal consensus (5′ → 3′)")
ax.set_ylabel("deletion starts / uniform")
ax.set_title(f"Where internal deletions start in intact elements ({len(pos):,} distinct deletions from {n_all:,} junctions; dashed = uniform)")
axstyle(ax)
save(fig, "deletion_positions.png")
print("figures in", OUT, sorted(os.listdir(OUT)))
