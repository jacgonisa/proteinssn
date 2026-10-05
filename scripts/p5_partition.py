#!/usr/bin/env python3
"""Phase 5 — partition-level view (supporting figure, not the primary statistic).

Leiden community detection per layer on the shared primary node set, then:
  * NMI / ARI for every layer pair (one-number "how well families hold");
  * a cluster-flow table GAG->RT->INT (Sankey input) + heatmaps;
  * multiplex Leiden consensus and a list of community-switching nodes.

Cluster boundaries are threshold artefacts, so candidate calls (Phase 4) must not
depend on this; it is a figure.

Usage:
    p5_partition.py --edges-dir DIR --layers GAG PROT RT RH INT \
        --nodes analysis_nodes.tsv --out-dir OUT [--resolution 1.0]
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
from collections import Counter, defaultdict

import numpy as np


def read_nodes(path):
    nodes = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            s = line.strip()
            if s:
                nodes.append(s.split("\t")[0])
    return nodes


def build_graph(edges_path, node_index):
    import igraph as ig
    edges, weights = [], []
    with open(edges_path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b, ratio = f[0], f[1], float(f[4])
            if a in node_index and b in node_index and a != b:
                edges.append((node_index[a], node_index[b]))
                weights.append(ratio)
    g = ig.Graph(n=len(node_index), edges=edges)
    g.es["weight"] = weights
    return g


def leiden(g, resolution):
    import leidenalg as la
    part = la.find_partition(g, la.RBConfigurationVertexPartition,
                             weights="weight", resolution_parameter=resolution,
                             seed=0)
    return np.array(part.membership)


def contingency(a, b):
    labels_a = sorted(set(a)); labels_b = sorted(set(b))
    ia = {l: i for i, l in enumerate(labels_a)}
    ib = {l: i for i, l in enumerate(labels_b)}
    m = np.zeros((len(labels_a), len(labels_b)), dtype=float)
    for x, y in zip(a, b):
        m[ia[x], ib[y]] += 1
    return m


def nmi(a, b):
    m = contingency(a, b)
    n = m.sum()
    if n == 0:
        return float("nan")
    pa = m.sum(1) / n; pb = m.sum(0) / n
    mi = 0.0
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if m[i, j] > 0:
                mi += (m[i, j] / n) * math.log((m[i, j] / n) / (pa[i] * pb[j]))
    ha = -sum(p * math.log(p) for p in pa if p > 0)
    hb = -sum(p * math.log(p) for p in pb if p > 0)
    denom = math.sqrt(ha * hb)
    return mi / denom if denom > 0 else 1.0


def ari(a, b):
    m = contingency(a, b)
    n = m.sum()
    comb2 = lambda x: x * (x - 1) / 2
    sum_ij = sum(comb2(v) for v in m.flatten())
    sum_a = sum(comb2(v) for v in m.sum(1))
    sum_b = sum(comb2(v) for v in m.sum(0))
    expected = sum_a * sum_b / comb2(n) if n > 1 else 0
    max_index = 0.5 * (sum_a + sum_b)
    denom = max_index - expected
    return (sum_ij - expected) / denom if denom != 0 else 1.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges-dir", required=True)
    ap.add_argument("--layers", nargs="+",
                    default=["GAG", "PROT", "RT", "RH", "INT"])
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--edge-suffix", default=".edges.tsv")
    ap.add_argument("--resolution", type=float, default=1.0)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "phase5")
    os.makedirs(out, exist_ok=True)

    nodes = read_nodes(args.nodes)
    node_index = {n: i for i, n in enumerate(nodes)}
    print(f"[p5] nodes {len(nodes)}; layers {args.layers}")

    membership = {}
    for lyr in args.layers:
        g = build_graph(os.path.join(args.edges_dir, lyr + args.edge_suffix),
                        node_index)
        membership[lyr] = leiden(g, args.resolution)
        print(f"[p5] {lyr}: {len(set(membership[lyr]))} communities")

    # partitions table
    with open(os.path.join(out, "partitions.tsv"), "w") as f:
        f.write("ns_id\t" + "\t".join(args.layers) + "\n")
        for n in nodes:
            i = node_index[n]
            f.write(n + "\t" + "\t".join(str(membership[l][i])
                    for l in args.layers) + "\n")

    # NMI / ARI
    pairs = list(itertools.combinations(args.layers, 2))
    with open(os.path.join(out, "nmi_ari.tsv"), "w") as f:
        f.write("layer_a\tlayer_b\tNMI\tARI\n")
        nmi_mat = np.ones((len(args.layers), len(args.layers)))
        for (A, B) in pairs:
            v_nmi = nmi(membership[A], membership[B])
            v_ari = ari(membership[A], membership[B])
            f.write(f"{A}\t{B}\t{v_nmi:.4f}\t{v_ari:.4f}\n")
            ia, ib = args.layers.index(A), args.layers.index(B)
            nmi_mat[ia, ib] = nmi_mat[ib, ia] = v_nmi
    _heatmap(nmi_mat, args.layers, os.path.join(out, "nmi_heatmap.png"),
             "NMI between layer partitions")

    # cluster flow (Sankey input) for consecutive pol-informative layers
    flow_chain = [l for l in ["GAG", "RT", "INT"] if l in args.layers]
    for A, B in zip(flow_chain, flow_chain[1:]):
        flow = Counter()
        for n in nodes:
            i = node_index[n]
            flow[(membership[A][i], membership[B][i])] += 1
        with open(os.path.join(out, f"flow_{A}_{B}.tsv"), "w") as f:
            f.write(f"{A}_comm\t{B}_comm\tn\n")
            for (ca, cb), c in flow.most_common():
                f.write(f"{ca}\t{cb}\t{c}\n")

    # multiplex consensus + switchers
    switchers = _multiplex(args, node_index, nodes, membership, out)

    with open(os.path.join(out, "phase5_summary.md"), "w") as f:
        f.write("# Phase 5 — partition view\n\n")
        f.write(f"- nodes: {len(nodes)}; resolution: {args.resolution}\n")
        for lyr in args.layers:
            f.write(f"- {lyr}: {len(set(membership[lyr]))} communities\n")
        f.write("\n## Mean NMI / ARI (family coherence)\n\n")
        f.write("| pair | NMI | ARI |\n|---|---:|---:|\n")
        for (A, B) in pairs:
            f.write(f"| {A}-{B} | {nmi(membership[A], membership[B]):.3f} | "
                    f"{ari(membership[A], membership[B]):.3f} |\n")
        f.write(f"\n- community-switching nodes (multiplex): **{switchers}**\n")
        f.write("\n> Supporting figure; the statistics are in Phase 4.\n")

    print(f"[p5] wrote partitions, nmi_ari, flows, multiplex to {out}")


def _multiplex(args, node_index, nodes, membership, out):
    try:
        import igraph as ig
        import leidenalg as la
    except Exception:
        return 0
    graphs = []
    for lyr in args.layers:
        graphs.append(build_graph(
            os.path.join(args.edges_dir, lyr + args.edge_suffix), node_index))
    parts = [la.RBConfigurationVertexPartition(g, weights="weight")
             for g in graphs]
    optimiser = la.Optimiser()
    optimiser.optimise_partition_multiplex(parts, n_iterations=-1)
    consensus = np.array(parts[0].membership)

    # a node switches if, in some layer, its layer-community maps (by majority)
    # to a consensus community different from its own consensus membership.
    n_switch = 0
    with open(os.path.join(out, "multiplex.tsv"), "w") as f:
        f.write("ns_id\tconsensus\tis_switcher\n")
        layer_maps = []
        for lyr in args.layers:
            comm2cons = defaultdict(Counter)
            for n in nodes:
                i = node_index[n]
                comm2cons[membership[lyr][i]][consensus[i]] += 1
            layer_maps.append({c: cc.most_common(1)[0][0]
                               for c, cc in comm2cons.items()})
        for n in nodes:
            i = node_index[n]
            mapped = {layer_maps[k][membership[lyr][i]]
                      for k, lyr in enumerate(args.layers)}
            switch = int(len(mapped) > 1 or
                         any(m != consensus[i] for m in mapped))
            n_switch += switch
            f.write(f"{n}\t{consensus[i]}\t{switch}\n")
    return n_switch


def _heatmap(mat, labels, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                    color="white" if mat[i, j] < 0.6 else "black", fontsize=9)
    ax.set_title(title)
    fig.colorbar(im, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
