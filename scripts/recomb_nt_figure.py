#!/usr/bin/env python3
"""Figure for a within-Athila recombinant (nucleotide per-domain SSN).

  Panel A  domain architecture, boxes coloured by the element's per-domain ATHILA
           family (nucleotide assignment), strand-aware.
  Panels B+  one nucleotide SSN per chosen domain: a family-stratified sample of
           all elements' domain-NT, coloured by nucleotide family, with the target
           starred — showing it sits in different family clusters in different domains.

Usage:
    recomb_nt_figure.py --assign domain_family_nt.tsv --data-dir DIR \
        --fulllength DIR --accession ACC --element ELEM --domains GAG RT \
        --out FIG.png [--per-family 70]
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import tempfile
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
from matplotlib.lines import Line2D
import networkx as nx

import lib_ids as L
from athila_recomb_nt import read_fasta, revcomp

FAM_COLORS = {
    "ATHILA2": "#4C72B0", "ATHILA4c": "#C44E52", "ATHILA1": "#55A868",
    "ATHILA6b": "#DD8452", "ATHILA6a": "#8172B3", "ATHILA7": "#17becf",
    "ATHILA7a": "#DA8BC3", "ATHILA0": "#937860", "ATHILA5": "#CCB974",
    "ATHILA3": "#64B5CD", "ATHILA4": "#aa8800", "ATHILA8a": "#999999",
    "ATHILA8b": "#bbbbbb", "ATHILA9": "#e377c2",
}
GREY = "#CCCCCC"


def fam_color(f):
    return FAM_COLORS.get(f, GREY)


def load_assign(path):
    a = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            ns, dom, fam, b, b2, mg, conf = line.rstrip("\n").split("\t")
            a[(ns, dom)] = (fam, int(conf))
    return a


def extract_domain_nt(data_dir, fl_dir, wanted):
    """wanted: set of (ns, domain). Returns {(ns,domain): seq}."""
    by_acc = defaultdict(set)
    for ns, dom in wanted:
        acc, _ = L.split_nsid(ns)
        by_acc[acc].add((ns, dom))
    out = {}
    for acc, items in by_acc.items():
        flf = glob.glob(os.path.join(fl_dir, f"{acc}*DP_FULLLENGTH_RENAMED.fasta"))
        gff = glob.glob(os.path.join(data_dir, f"{acc}*.dom.gff3"))
        if not flf or not gff:
            continue
        seqs = read_fasta(flf[0])
        want_doms = {(L.split_nsid(ns)[1], dom) for ns, dom in items}
        for line in open(gff[0]):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            dom = attrs.get("gene")
            el = f[0]
            if (el, dom) not in want_doms or el not in seqs:
                continue
            s, e, strand = int(f[3]), int(f[4]), f[6]
            sub = seqs[el][s - 1:e]
            if strand == "-":
                sub = revcomp(sub)
            out[(L.nsid(acc, el), dom)] = sub
    return out


def nt_network(records, workdir, threads=6, min_ratio=0.30):
    """records: {id: seq}. mmseqs all-vs-all nt -> graph with bitscore-ratio edges."""
    fa = os.path.join(workdir, "n.fa")
    with open(fa, "w") as o:
        for k, s in records.items():
            o.write(f">{k}\n{s}\n")
    res = os.path.join(workdir, "n.m8")
    subprocess.run(["mmseqs", "easy-search", fa, fa, res, os.path.join(workdir, "t"),
                    "--search-type", "3", "-e", "1e-3", "--max-seqs", "2000",
                    "--threads", str(threads),
                    "--format-output", "query,target,bits", "-v", "1"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    selfbit, pair = {}, {}
    with open(res) as fh:
        for line in fh:
            q, t, b = line.rstrip("\n").split("\t"); b = float(b)
            if q == t:
                selfbit[q] = max(selfbit.get(q, 0), b)
            else:
                k = (q, t) if q < t else (t, q)
                pair[k] = max(pair.get(k, 0), b)
    g = nx.Graph(); g.add_nodes_from(records)
    for (a, b), bit in pair.items():
        denom = min(selfbit.get(a, bit), selfbit.get(b, bit)) or bit
        r = min(1.0, bit / denom)
        if r >= min_ratio:
            g.add_edge(a, b, weight=r)
    return g


def draw_arch(ax, element, data_dir, acc, assign):
    gff = glob.glob(os.path.join(data_dir, f"{acc}*.dom.gff3"))[0]
    doms = []
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if not f[0] == element or len(f) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        fam = assign.get((L.nsid(acc, element), attrs.get("gene")), (None,))[0]
        doms.append((int(f[3]), int(f[4]), f[6], attrs.get("gene"), fam,
                     attrs.get("clade")))
    doms.sort()
    span = max(d[1] for d in doms)
    ax.plot([0, span], [0, 0], color="#333", lw=1.5, zorder=1)
    for s, e, strand, gene, fam, clade in doms:
        col = fam_color(fam) if fam else GREY
        w = e - s; dx = w if strand == "+" else -w
        x0 = s if strand == "+" else e
        ax.add_patch(FancyArrow(x0, 0, dx, 0, width=0.22, head_width=0.32,
                     head_length=min(abs(dx) * 0.35, span * 0.02),
                     length_includes_head=True, color=col, edgecolor="white",
                     linewidth=0.6, zorder=3))
        ax.text((s + e) / 2, 0.42, f"{gene}\n{fam or clade}", ha="center",
                va="bottom", fontsize=9, fontweight="bold", color=col)
    ax.set_xlim(-span * 0.03, span * 1.03); ax.set_ylim(-0.8, 1.1)
    ax.set_yticks([]); ax.set_xlabel("position in element (bp)")
    ax.set_title(f"A  Domain architecture (nt family) — {element}", loc="left",
                 fontsize=12, fontweight="bold")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)


def layout_positions(g, seed=1):
    """Scalable layout: igraph FR (small) / DrL (large); fall back to networkx."""
    nodes = list(g.nodes())
    if len(nodes) > 1200:
        try:
            import igraph as ig
            idx = {n: i for i, n in enumerate(nodes)}
            edges = [(idx[u], idx[v]) for u, v in g.edges()]
            w = [g[u][v]["weight"] for u, v in g.edges()]
            ig_g = ig.Graph(n=len(nodes), edges=edges)
            lay = (ig_g.layout_drl() if len(nodes) > 5000
                   else ig_g.layout_fr(weights=w))
            return {nodes[i]: tuple(lay[i]) for i in range(len(nodes))}
        except Exception as e:
            print(f"  (igraph layout failed: {e}; using networkx)")
    return nx.spring_layout(g, seed=seed, weight="weight", k=0.35)


def draw_ssn(ax, g, fam_of, target, title, letter, seed=1):
    pos = layout_positions(g, seed=seed)
    nx.draw_networkx_edges(g, pos, ax=ax, alpha=0.12, width=0.4, edge_color="#999")
    fams = sorted(set(fam_of.values()))
    for fam in fams:
        ns = [n for n in g if fam_of[n] == fam and n != target]
        if ns:
            ax.scatter([pos[n][0] for n in ns], [pos[n][1] for n in ns], s=38,
                       c=fam_color(fam), edgecolors="white", linewidths=0.3,
                       label=f"{fam} ({sum(1 for v in fam_of.values() if v==fam)})",
                       zorder=2)
    if target in pos:
        ax.scatter([pos[target][0]], [pos[target][1]], s=560, marker="*",
                   c="#111", edgecolors="yellow", linewidths=1.8, zorder=6,
                   label="recombinant")
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=12, fontweight="bold")
    ax.axis("off")
    ax.legend(loc="best", fontsize=7.5, framealpha=0.9)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--fulllength", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--domains", nargs="+", default=["GAG", "RT"])
    ap.add_argument("--families", nargs="*", default=None,
                    help="restrict SSN nodes to these families for clarity")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-family", type=int, default=70)
    ap.add_argument("--min-ratio", type=float, default=0.30,
                    help="edge bitscore-ratio threshold (higher = cleaner clusters)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    fam_filter = set(args.families) if args.families else None

    assign = load_assign(args.assign)
    target = L.nsid(args.accession, args.element)

    import random
    rng = random.Random(args.seed)
    panels = {}
    wanted = set()
    for dom in args.domains:
        byfam = defaultdict(list)
        for (ns, d), (fam, conf) in assign.items():
            if d == dom and conf and (fam_filter is None or fam in fam_filter):
                byfam[fam].append(ns)
        chosen = []
        for fam, members in byfam.items():
            rng.shuffle(members)
            chosen += [(ns, dom) for ns in members[:args.per_family]]
        chosen.append((target, dom))
        panels[dom] = chosen
        wanted |= set(chosen)

    seqs = extract_domain_nt(args.data_dir, args.fulllength, wanted)
    print(f"[fig] extracted {len(seqs)} domain-NT sequences")

    ncol = len(args.domains)
    fig = plt.figure(figsize=(7 * ncol, 11))
    gs = fig.add_gridspec(2, ncol, height_ratios=[1, 2.4], hspace=0.28, wspace=0.12)
    draw_arch(fig.add_subplot(gs[0, :]), args.element, args.data_dir,
              args.accession, assign)

    with tempfile.TemporaryDirectory() as tmp:
        for col, dom in enumerate(args.domains):
            recs = {f"{ns}": seqs[(ns, dom)] for ns, d in panels[dom]
                    if (ns, dom) in seqs}
            fam_of = {ns: assign[(ns, dom)][0] for ns in recs}
            g = nt_network(recs, tmp, min_ratio=args.min_ratio)
            tfam = assign.get((target, dom), ("?",))[0]
            print(f"[fig] {dom}: {len(recs)} seqs, target family={tfam}")
            draw_ssn(fig.add_subplot(gs[1, col]), g, fam_of, target,
                     f"{dom} nucleotide SSN", chr(66 + col), seed=col + 1)

    arch = ";".join(f"{d}:{assign.get((target,d),('?',))[0]}"
                    for d in ["GAG", "PROT", "RT", "RH", "INT"]
                    if (target, d) in assign)
    fig.suptitle(f"Within-Athila recombinant  {args.accession} : {args.element}\n{arch}",
                 fontsize=14, fontweight="bold")
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"[fig] wrote {args.out}")


if __name__ == "__main__":
    main()
