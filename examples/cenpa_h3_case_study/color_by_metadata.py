#!/usr/bin/env python3
"""Colour a proteinssn network by an external metadata table (e.g. taxonomic clade),
which the built-in `ssn plot` (colours by cluster) does not do.

Usage:
    python color_by_metadata.py EDGES.tsv META.tsv OUT.png \
        [--color-col clade] [--shape-col group] [--layout organic]

META.tsv: tab-separated, first column = sequence name (matching the FASTA / edge list),
plus a colour column and (optionally) a shape column.

Layouts:  organic = ForceAtlas2 (Gephi/Cytoscape-style, the usual SSN look); spring = fast fallback.
"""
import argparse, itertools
import pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PALETTE = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf","#999999"]
# sensible fixed colours for common taxonomic clades (fall back to PALETTE otherwise)
CLADE_COLORS = {"Vertebrates":"#e41a1c","Invertebrates":"#377eb8","Viridiplantae":"#4daf4a",
                "Fungi":"#984ea3","Protist":"#ff7f00"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("edges"); ap.add_argument("meta"); ap.add_argument("out")
    ap.add_argument("--color-col", default=None, help="metadata column for node colour (default: 2nd col)")
    ap.add_argument("--shape-col", default=None, help="optional metadata column for node shape")
    ap.add_argument("--layout", choices=["organic","spring"], default="organic")
    ap.add_argument("--giant-only", action="store_true",
                    help="plot only the largest connected component (drops singletons / tiny clusters "
                         "that ForceAtlas2 flings to the periphery)")
    ap.add_argument("--title", default="Sequence similarity network")
    a = ap.parse_args()

    ed = pd.read_csv(a.edges, sep="\t")
    src, tgt = ed.columns[0], ed.columns[1]
    meta = pd.read_csv(a.meta, sep="\t").set_index(pd.read_csv(a.meta, sep="\t").columns[0])
    color_col = a.color_col or meta.columns[0]
    G = nx.from_pandas_edgelist(ed, src, tgt)
    if a.giant_only:
        giant = max(nx.connected_components(G), key=len)
        G = G.subgraph(giant).copy()
    else:
        G.add_nodes_from(meta.index)                   # include singletons

    cats = sorted(meta[color_col].dropna().unique())
    pal = {c: CLADE_COLORS.get(c, PALETTE[i % len(PALETTE)]) for i, c in enumerate(cats)}
    col = meta[color_col].to_dict()

    print(f"layout: {a.layout} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)…", flush=True)
    if a.layout == "organic":
        pos = nx.forceatlas2_layout(G, max_iter=300, scaling_ratio=2.0, seed=1)
    else:
        pos = nx.spring_layout(G, k=0.15, iterations=50, seed=1)

    fig, ax = plt.subplots(figsize=(13, 13))
    nx.draw_networkx_edges(G, pos, alpha=0.05, width=0.3, edge_color="#888", ax=ax)
    if a.shape_col:
        shp = meta[a.shape_col].to_dict()
        markers = dict(zip(sorted(set(shp.values())), itertools.cycle(["o","^","s","D","v"])))
        for sval, mk in markers.items():
            ns = [n for n in G if shp.get(n) == sval]
            nx.draw_networkx_nodes(G, pos, nodelist=ns, node_shape=mk,
                node_size=26 if mk != "o" else 10,
                node_color=[pal.get(col.get(n), "#999") for n in ns],
                linewidths=0.2, edgecolors="#333", ax=ax)
    else:
        nx.draw_networkx_nodes(G, pos, node_size=14,
            node_color=[pal.get(col.get(n), "#999") for n in G], linewidths=0.2, edgecolors="#333", ax=ax)

    ax.set_axis_off(); ax.set_title(a.title, fontsize=13)
    legend = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=10, label=k)
              for k, c in pal.items()]
    if a.shape_col:
        legend += [Line2D([0],[0], marker=mk, color='w', markerfacecolor='#555', markersize=10, label=s)
                   for s, mk in markers.items()]
    ax.legend(handles=legend, loc="upper left", fontsize=10, frameon=True)
    plt.tight_layout()
    fig.savefig(a.out, dpi=200, bbox_inches="tight", facecolor="white")
    print("saved", a.out)

if __name__ == "__main__":
    main()
