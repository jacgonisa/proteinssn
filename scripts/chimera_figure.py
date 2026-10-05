#!/usr/bin/env python3
"""Make a chimera figure for one element:

  Panel A  domain architecture along the element (boxes coloured by TEsorter clade,
           arrows = strand) — shows the Athila vs foreign modules.
  Panel B  INT sequence-similarity network (Athila sample + all non-Athila INTs),
           coloured by clade, with the target element starred — it lands in the
           foreign cluster.
  Panel C  GAG network, same idea — the same element stays in the Athila cluster.

Usage:
    chimera_figure.py --data-dir DIR --accession ACC --element ELEM --out FIG.png
        [--diamond ~/bin/diamond] [--sample-athila 150]
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import networkx as nx

import lib_ids as L

CLADE_COLORS = {
    "Athila": "#4C72B0", "Retand": "#C44E52", "Tekay": "#DD8452",
    "Tork": "#55A868", "Ivana": "#8172B3", "CRM": "#17becf",
    "Reina": "#937860", "Ale": "#DA8BC3",
}
OTHER = "#BBBBBB"


def clade_color(c):
    return CLADE_COLORS.get(c, OTHER)


def find_file(data_dir, acc, ext):
    hits = glob.glob(os.path.join(data_dir, f"{acc}*{ext}"))
    if not hits:
        sys.exit(f"no {ext} for accession {acc}")
    return hits[0]


# --------------------------------------------------------------------------- #
# Panel A: architecture from dom.gff3
# --------------------------------------------------------------------------- #
def parse_arch(gff3, element):
    doms = []
    with open(gff3) as fh:
        for line in fh:
            if not line.startswith(element + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            doms.append({"start": int(f[3]), "end": int(f[4]), "strand": f[6],
                         "gene": attrs.get("gene", "?"),
                         "clade": attrs.get("clade", "?"),
                         "evalue": attrs.get("evalue", "")})
    return sorted(doms, key=lambda d: d["start"])


def draw_arch(ax, element, doms):
    span = max(d["end"] for d in doms)
    ax.plot([0, span], [0, 0], color="#333333", lw=1.5, zorder=1)
    for d in doms:
        col = clade_color(d["clade"])
        w = d["end"] - d["start"]
        dx = w if d["strand"] == "+" else -w
        x0 = d["start"] if d["strand"] == "+" else d["end"]
        ax.add_patch(FancyArrow(
            x0, 0, dx, 0, width=0.22, head_width=0.32,
            head_length=min(abs(dx) * 0.35, span * 0.02),
            length_includes_head=True, color=col, zorder=3,
            edgecolor="white", linewidth=0.6))
        ax.text((d["start"] + d["end"]) / 2, 0.42,
                f"{d['gene']}\n{d['clade']}", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=col)
    ax.set_xlim(-span * 0.03, span * 1.03)
    ax.set_ylim(-0.8, 1.1)
    ax.set_yticks([])
    ax.set_xlabel("position in element (bp)")
    ax.set_title(f"A  Domain architecture — {element}", loc="left",
                 fontsize=12, fontweight="bold")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


# --------------------------------------------------------------------------- #
# Panels B/C: per-domain SSN with cross-clade references
# --------------------------------------------------------------------------- #
def collect_domain(data_dir, domain, target_ns, sample_athila, seed=0):
    """Longest peptide per element for `domain`, keep all non-Athila + a sample of
    Athila + the target. Returns (fasta_records dict ns->seq, clade map)."""
    import random
    rng = random.Random(seed)
    seqs, clade = {}, {}
    for ff in sorted(glob.glob(os.path.join(data_dir, "*.dom.faa"))):
        acc = L.accession_from_path(ff)
        for element, dom, attrs, seq in L.parse_dom_faa(ff):
            if dom != domain:
                continue
            ns = L.nsid(acc, element)
            if ns not in seqs or len(seq) > len(seqs[ns]):
                seqs[ns] = seq
                clade[ns] = attrs.get("clade", "?")
    non_ath = [n for n, c in clade.items() if c != "Athila"]
    ath = [n for n, c in clade.items() if c == "Athila"]
    rng.shuffle(ath)
    keep = set(non_ath) | set(ath[:sample_athila]) | {target_ns}
    keep &= set(seqs)
    return {n: seqs[n] for n in keep}, {n: clade[n] for n in keep}


def build_network(records, diamond, workdir):
    faa = os.path.join(workdir, "d.faa")
    with open(faa, "w") as o:
        for n, s in records.items():
            o.write(f">{n}\n{s}\n")
    db = os.path.join(workdir, "d")
    raw = os.path.join(workdir, "d.tsv")
    subprocess.run([diamond, "makedb", "--in", faa, "--db", db, "--quiet"],
                   check=True)
    subprocess.run([diamond, "blastp", "-q", faa, "--db", db, "-o", raw,
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore",
                    "--very-sensitive", "--max-target-seqs", "5000",
                    "--evalue", "1e-3", "--quiet"], check=True)
    selfbit, pair = {}, {}
    with open(raw) as fh:
        for line in fh:
            q, s, b = line.split("\t"); b = float(b)
            if q == s:
                selfbit[q] = max(selfbit.get(q, 0), b)
            else:
                k = (q, s) if q < s else (s, q)
                pair[k] = max(pair.get(k, 0), b)
    g = nx.Graph()
    g.add_nodes_from(records)
    for (a, b), bit in pair.items():
        denom = min(selfbit.get(a, bit), selfbit.get(b, bit)) or bit
        ratio = min(1.0, bit / denom)
        if ratio >= 0.20:           # keep structure legible
            g.add_edge(a, b, weight=ratio)
    return g


def draw_ssn(ax, g, clade, target_ns, title, seed=1):
    pos = nx.spring_layout(g, seed=seed, weight="weight", k=0.3)
    nx.draw_networkx_edges(g, pos, ax=ax, alpha=0.15, width=0.4,
                           edge_color="#999999")
    for c in sorted(set(clade.values())):
        ns = [n for n in g if clade[n] == c and n != target_ns]
        if not ns:
            continue
        ax.scatter([pos[n][0] for n in ns], [pos[n][1] for n in ns],
                   s=45, c=clade_color(c), edgecolors="white", linewidths=0.4,
                   label=f"{c} ({sum(1 for x in clade.values() if x==c)})", zorder=2)
    if target_ns in pos:
        ax.scatter([pos[target_ns][0]], [pos[target_ns][1]], s=520, marker="*",
                   c="#111111", edgecolors="yellow", linewidths=1.6, zorder=5,
                   label="chimeric element")
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.axis("off")
    ax.legend(loc="best", fontsize=8, framealpha=0.9)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--diamond", default=os.path.expanduser("~/bin/diamond"))
    ap.add_argument("--sample-athila", type=int, default=150)
    args = ap.parse_args()

    target_ns = L.nsid(args.accession, args.element)
    gff3 = find_file(args.data_dir, args.accession, ".dom.gff3")
    arch = parse_arch(gff3, args.element)
    if not arch:
        sys.exit(f"element {args.element} not found in {gff3}")

    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.4], hspace=0.28, wspace=0.12)
    draw_arch(fig.add_subplot(gs[0, :]), args.element, arch)

    with tempfile.TemporaryDirectory() as tmp:
        for col, (dom, title) in enumerate([
                ("INT", "B  INT network (integrase) — foreign module"),
                ("GAG", "C  GAG network (capsid) — Athila module")]):
            recs, clade = collect_domain(args.data_dir, dom, target_ns,
                                         args.sample_athila)
            g = build_network(recs, args.diamond, tmp)
            present = "yes" if target_ns in recs else "NO"
            print(f"[fig] {dom}: {len(recs)} seqs, target present={present}, "
                  f"target clade-call={clade.get(target_ns,'NA')}")
            draw_ssn(fig.add_subplot(gs[1, col]), g, clade, target_ns, title,
                     seed=col + 1)

    fig.suptitle(f"Chimeric Athila element  {args.accession} : {args.element}\n"
                 f"Athila capsid/protease (5′) + Retand integrase (3′)",
                 fontsize=14, fontweight="bold")
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"[fig] wrote {args.out}")


if __name__ == "__main__":
    main()
