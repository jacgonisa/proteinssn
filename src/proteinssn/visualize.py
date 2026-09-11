"""Static (matplotlib) and interactive (pyvis) renderings of an SSN."""

from __future__ import annotations

import math
import sys

from .network import NetworkResult, cluster_sizes

# A modern, reasonably colourblind-aware palette. Clusters beyond it are grey.
PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#EEC643", "#DA8BC3", "#00A6A6", "#B07AA1", "#5FA2CE",
    "#F28E2B", "#59A14F", "#E15759", "#9C755F", "#76B7B2",
]
GREY = "#C9CBD3"
BG = "#FBFBFD"          # near-white canvas
EDGE_RGB = (0.36, 0.40, 0.52)  # slate, alpha applied per-edge


def _cluster_colors(result: NetworkResult, top_n: int) -> dict[str, str]:
    """Colour the `top_n` biggest real clusters; singletons/overflow stay grey."""
    colors, i = {}, 0
    for label, size in cluster_sizes(result):
        if size > 1 and i < top_n:
            colors[label] = PALETTE[i % len(PALETTE)]
            i += 1
        else:
            colors[label] = GREY
    return colors


def _packed_layout(graph, seed: int) -> dict:
    """Lay out each connected component, then shelf-pack them into a tidy block.

    Plain ``spring_layout`` lets disconnected clusters drift far apart, leaving
    most of the canvas empty. Here each component is laid out on its own and
    scaled by ``sqrt(size)``, then packed left-to-right / top-to-bottom so the
    figure stays compact and every cluster reads clearly.
    """
    import networkx as nx

    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    gap = 1.6
    boxes = []  # (positions dict, radius)
    for comp in components:
        if len(comp) == 1:
            node = next(iter(comp))
            boxes.append(({node: (0.0, 0.0)}, 0.5))
            continue
        sub = graph.subgraph(comp)
        pos = nx.spring_layout(sub, seed=seed, k=1.4 / math.sqrt(len(comp)),
                               iterations=80)
        cx = sum(p[0] for p in pos.values()) / len(pos)
        cy = sum(p[1] for p in pos.values()) / len(pos)
        maxr = max(math.hypot(x - cx, y - cy) for x, y in pos.values()) or 1.0
        scale = math.sqrt(len(comp))
        pos = {n: ((x - cx) / maxr * scale, (y - cy) / maxr * scale)
               for n, (x, y) in pos.items()}
        boxes.append((pos, scale))

    widths = [2 * r + gap for _, r in boxes]
    target = max(max(widths), math.sqrt(sum(w * w for w in widths)) * 1.5)

    positions: dict = {}
    x = y = row_h = 0.0
    for (pos, r), w in zip(boxes, widths):
        if x > 1e-9 and x + w > target:
            x, y, row_h = 0.0, y - row_h, 0.0
        cx, cy = x + w / 2, y - w / 2
        for n, (px, py) in pos.items():
            positions[n] = (cx + px, cy + py)
        x += w
        row_h = max(row_h, w)
    return positions


def static_plot(
    result: NetworkResult,
    path: str,
    top_n: int = 14,
    max_nodes: int = 3000,
    seed: int = 42,
    dpi: int = 220,
) -> bool:
    """Publication-style PNG/SVG coloured by cluster. Returns False if skipped."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import networkx as nx

    graph = result.graph
    if graph.number_of_nodes() == 0:
        print("  (nothing to plot: network is empty)", file=sys.stderr)
        return False
    if graph.number_of_nodes() > max_nodes:
        print(
            f"  skipping static plot: {graph.number_of_nodes()} nodes "
            f"> --max-nodes {max_nodes} (use the interactive HTML instead).",
            file=sys.stderr,
        )
        return False

    colors = _cluster_colors(result, top_n)
    pos = _packed_layout(graph, seed)

    # Node size grows (gently) with connectivity.
    deg = dict(graph.degree())
    node_size = [70 + 26 * math.sqrt(deg[n]) for n in graph.nodes()]
    node_colors = [colors[result.clusters[n]] for n in graph.nodes()]

    # Edge shade/width scale with % identity.
    idents = [d.get("pident", 100.0) for *_e, d in graph.edges(data=True)]
    lo, hi = (min(idents), max(idents)) if idents else (0.0, 100.0)
    span = (hi - lo) or 1.0
    edge_colors, edge_widths = [], []
    for *_e, d in graph.edges(data=True):
        t = (d.get("pident", 100.0) - lo) / span
        edge_colors.append((*EDGE_RGB, 0.12 + 0.45 * t))
        edge_widths.append(0.4 + 1.6 * t)

    fig, ax = plt.subplots(figsize=(13, 13))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color=edge_colors,
                           width=edge_widths)
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors,
                           node_size=node_size, linewidths=0.6,
                           edgecolors="white")

    # Legend: biggest real clusters, then singletons / overflow grouped in grey.
    multi = [(label, size) for label, size in cluster_sizes(result) if size > 1]
    shown = multi[:top_n]

    def _handle(color, text):
        return Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color,
                      markeredgecolor="white", markersize=11, label=text)

    handles = [_handle(colors[label], f"{label}  ({size})")
               for label, size in shown]
    extra = len(multi) - len(shown)
    if extra > 0:
        handles.append(_handle(GREY, f"+{extra} more clusters"))
    if result.singletons:
        handles.append(_handle(GREY, f"singletons  ({len(result.singletons)})"))
    legend = ax.legend(handles=handles, loc="upper left", frameon=True,
                       fontsize=10, title="clusters", title_fontsize=11,
                       borderpad=0.8, labelspacing=0.6)
    legend.get_frame().set_edgecolor("#DfE1E8")
    legend.get_frame().set_facecolor("white")

    ax.set_title(
        f"Sequence Similarity Network\n"
        f"{graph.number_of_nodes()} proteins · {graph.number_of_edges()} edges · "
        f"{result.n_clusters} clusters",
        fontsize=16, fontweight="bold", color="#22252B", pad=18,
    )
    ax.margins(0.04)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"  static network   -> {path}", file=sys.stderr)
    return True


_INTERACTIVE_OPTIONS = """
{
  "nodes": {
    "borderWidth": 2,
    "borderWidthSelected": 4,
    "color": { "border": "#ffffff", "highlight": { "border": "#22252B" } },
    "shadow": { "enabled": true, "size": 12, "x": 0, "y": 2,
                "color": "rgba(0,0,0,0.20)" },
    "font": { "color": "#e8e8ef", "size": 12, "face": "Helvetica" }
  },
  "edges": {
    "color": { "color": "#7c8296", "opacity": 0.5,
               "highlight": "#ffd166", "hover": "#ffd166" },
    "smooth": { "type": "continuous", "roundness": 0.15 },
    "scaling": { "min": 0.4, "max": 5 }
  },
  "interaction": { "hover": true, "tooltipDelay": 80,
                   "hideEdgesOnDrag": true, "navigationButtons": true },
  "physics": {
    "solver": "forceAtlas2Based",
    "forceAtlas2Based": { "gravitationalConstant": -60,
                          "centralGravity": 0.008,
                          "springLength": 110, "springConstant": 0.08,
                          "damping": 0.5, "avoidOverlap": 0.6 },
    "stabilization": { "iterations": 220 },
    "minVelocity": 0.6
  }
}
"""


def interactive_plot(
    result: NetworkResult,
    path: str,
    top_n: int = 14,
    max_nodes: int = 20000,
) -> bool:
    """Self-contained interactive HTML (pan/zoom/hover). Returns False if skipped."""
    try:
        from pyvis.network import Network
    except ImportError:
        print(
            "  skipping interactive plot: pyvis not installed "
            "(pip install pyvis).",
            file=sys.stderr,
        )
        return False

    graph = result.graph
    if graph.number_of_nodes() == 0:
        print("  (nothing to plot: network is empty)", file=sys.stderr)
        return False
    if graph.number_of_nodes() > max_nodes:
        print(
            f"  skipping interactive plot: {graph.number_of_nodes()} nodes "
            f"> --max-nodes {max_nodes}.",
            file=sys.stderr,
        )
        return False

    colors = _cluster_colors(result, top_n)
    deg = dict(graph.degree())
    net = Network(height="100vh", width="100%", bgcolor="#12141C",
                  font_color="#e8e8ef", notebook=False, cdn_resources="in_line")
    net.set_options(_INTERACTIVE_OPTIONS)

    for node in graph.nodes():
        cluster = result.clusters[node]
        net.add_node(
            node, label=node,
            title=f"{node}\ncluster: {cluster}\nconnections: {deg[node]}",
            color=colors[cluster],
            size=10 + 3.0 * math.sqrt(deg[node]),
        )
    for u, v, data in graph.edges(data=True):
        title = f"{data['pident']:.1f}% identity"
        if "evalue" in data:
            title += f"\ne-value: {data['evalue']:.1e}"
        net.add_edge(u, v, value=float(data["pident"]), title=title)

    net.save_graph(path)
    print(f"  interactive html  -> {path}", file=sys.stderr)
    return True


def sweep_plot(rows: list[dict], path: str, dpi: int = 220) -> None:
    """Plot threshold-sweep results: #clusters and largest cluster vs identity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ident = [r["identity"] for r in rows]
    fig, ax1 = plt.subplots(figsize=(8.5, 5))
    fig.patch.set_facecolor(BG)
    ax1.set_facecolor(BG)

    c1, c2 = "#4C72B0", "#C44E52"
    ax1.plot(ident, [r["n_clusters"] for r in rows], "-o", color=c1,
             lw=2.2, ms=7, label="clusters")
    ax1.set_xlabel("identity threshold (%)", fontsize=11)
    ax1.set_ylabel("number of clusters", color=c1, fontsize=11)
    ax1.tick_params(axis="y", labelcolor=c1)
    ax1.grid(True, alpha=0.25, linestyle="--")

    ax2 = ax1.twinx()
    ax2.plot(ident, [r["largest_cluster"] for r in rows], "-s", color=c2,
             lw=2.2, ms=7, label="largest cluster")
    ax2.set_ylabel("size of largest cluster", color=c2, fontsize=11)
    ax2.tick_params(axis="y", labelcolor=c2)

    for spine in ("top",):
        ax1.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)

    ax1.set_title("Threshold sweep", fontsize=14, fontweight="bold",
                  color="#22252B", pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print(f"  sweep plot -> {path}", file=sys.stderr)
