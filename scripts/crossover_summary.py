#!/usr/bin/env python3
"""Summarise crossover calls: loci clustered across accessions, pair table, figures.

A locus = crossover calls with the same family pair, chromosome and breakpoint
interval whose starts lie within --link-bp of each other (single linkage); the
number of accessions carrying a locus is its replication.

Usage: crossover_summary.py RESULTS_DIR
"""
import os
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = sys.argv[1]
LINK = 1_000_000
OUT = os.path.join(R, "crossovers")
SURF, INK, MUTED, GRID, ACC = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0", "#2a78d6"

rows = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/crossover_calls.tsv")]
h = rows[0]
calls = [dict(zip(h, r)) for r in rows[1:]]
cross = [c for c in calls if c["class"] == "crossover"]

groups = defaultdict(list)
for c in cross:
    pair = "-".join(sorted((c["family_5p"], c["family_3p"])))
    groups[(pair, c["chrom"], c["breakpoint"])].append(c)
loci = []
for (pair, chrom, bp), cs in groups.items():
    cs.sort(key=lambda c: int(c["start"]))
    cur = [cs[0]]
    for c in cs[1:]:
        if int(c["start"]) - int(cur[-1]["start"]) <= LINK:
            cur.append(c)
        else:
            loci.append((pair, chrom, bp, cur)); cur = [c]
    loci.append((pair, chrom, bp, cur))

hom = {}
for l in open(f"{R}/bias/bias_normalizations.tsv").readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    hom["-".join(sorted(f[0].split("-")))] = float(f[2])
cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{R}/centro/family_centro_class.tsv").readlines()[1:]}

with open(f"{OUT}/crossover_loci.tsv", "w") as f:
    f.write("pair\tchrom\tstart_min_Mb\tstart_max_Mb\tbreakpoint\tn_elements\tn_accessions\tdirection\taccessions\n")
    for pair, chrom, bp, cs in sorted(loci, key=lambda x: -len({c["accession"] for c in x[3]})):
        accs = sorted({c["accession"] for c in cs})
        dirs = Counter(f"{c['family_5p']}>{c['family_3p']}" for c in cs).most_common(1)[0][0]
        st = [int(c["start"]) / 1e6 for c in cs]
        f.write(f"{pair}\t{chrom}\t{min(st):.2f}\t{max(st):.2f}\t{bp}\t{len(cs)}\t{len(accs)}\t{dirs}\t{','.join(accs)}\n")

pairs = defaultdict(lambda: Counter())
for pair, chrom, bp, cs in loci:
    n_acc = len({c["accession"] for c in cs})
    pairs[pair]["elements"] += len(cs)
    pairs[pair]["loci"] += 1
    pairs[pair]["replicated"] += n_acc >= 2
with open(f"{OUT}/crossover_pairs.tsv", "w") as f:
    f.write("pair\telements\tloci\treplicated_loci\thomology\tclass_a\tclass_b\n")
    for p, c in sorted(pairs.items(), key=lambda x: (-x[1]["replicated"], -x[1]["elements"])):
        a, b = p.split("-")
        f.write(f"{p}\t{c['elements']}\t{c['loci']}\t{c['replicated']}\t{hom.get(p, float('nan')):.0f}\t"
                f"{cls.get(a, 'NA')}\t{cls.get(b, 'NA')}\n")

# figure: pairs (loci, replicated)
ps = sorted(pairs, key=lambda p: (pairs[p]["replicated"], pairs[p]["loci"]))
fig, ax = plt.subplots(figsize=(8.5, 0.32 * len(ps) + 1.2), facecolor=SURF)
y = np.arange(len(ps))
ax.barh(y, [pairs[p]["loci"] for p in ps], color="#c9c8c1", height=0.6, label="all loci")
ax.barh(y, [pairs[p]["replicated"] for p in ps], color=ACC, height=0.6, label="replicated in ≥2 accessions")
for i, p in enumerate(ps):
    ax.text(pairs[p]["loci"] + 0.15, i, f"{pairs[p]['elements']} elements", va="center", fontsize=7.5, color=MUTED)
ax.set_yticks(y, [f"ATHILA{p.replace('-', ' × ')}" for p in ps], fontsize=8.5)
ax.set_xlabel("crossover loci"); ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("Crossover loci by family pair", loc="left", fontsize=11, color=INK)
ax.set_facecolor(SURF); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_pairs.png", dpi=150, facecolor=SURF); plt.close(fig)

# figure: breakpoint location, replicated loci counted once
order = ["between GAG and PROT", "between PROT and RT", "within RT", "within RH", "within INT",
         "outside annotated domains"]
bc_all = Counter(bp for pair, chrom, bp, cs in loci)
bc_rep = Counter(bp for pair, chrom, bp, cs in loci if len({c["accession"] for c in cs}) >= 2)
keys = [k for k in order if bc_all[k]] + [k for k in bc_all if k not in order]
fig, ax = plt.subplots(figsize=(8.5, 2.9), facecolor=SURF)
x = np.arange(len(keys))
ax.bar(x, [bc_all[k] for k in keys], 0.6, color="#c9c8c1", label="all loci")
ax.bar(x, [bc_rep[k] for k in keys], 0.6, color=ACC, label="replicated loci")
ax.set_xticks(x, [k.replace(" annotated", "\nannotated").replace("between ", "between\n") for k in keys], fontsize=8)
ax.set_ylabel("loci"); ax.legend(frameon=False, fontsize=8)
ax.set_title("Where the switch falls in the element", loc="left", fontsize=11, color=INK)
ax.set_facecolor(SURF); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_breakpoints.png", dpi=150, facecolor=SURF); plt.close(fig)

rep = [l for l in loci if len({c["accession"] for c in l[3]}) >= 2]
print(f"crossover calls {len(cross)}; loci {len(loci)}; replicated (>=2 accessions) {len(rep)}; "
      f">=3: {sum(1 for l in loci if len({c['accession'] for c in l[3]}) >= 3)}")
print(open(f"{OUT}/crossover_pairs.tsv").read())
