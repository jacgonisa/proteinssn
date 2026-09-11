"""Static (matplotlib) and interactive (pyvis) renderings of an SSN."""

from __future__ import annotations

import sys

from .network import NetworkResult, cluster_sizes

# A colourblind-friendly palette; clusters beyond it fall back to grey.
PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
]
GREY = "#D0D0D0"


def _cluster_colors(result: NetworkResult, top_n: int) -> dict[str, str]:
    """Assign a colour to the `top_n` biggest clusters; rest are grey."""
    ranked = [label for label, _ in cluster_sizes(result)]
    colors = {}
    for i, label in enumerate(ranked):
        colors[label] = PALETTE[i % len(PALETTE)] if i < top_n else GREY
    return colors


def static_plot(
    result: NetworkResult,
    path: str,
    top_n: int = 14,
    max_nodes: int = 3000,
    seed: int = 42,
    dpi: int = 200,
) -> bool:
    """Spring-layout PNG/SVG coloured by cluster. Returns False if skipped."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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
    node_colors = [colors[result.clusters[n]] for n in graph.nodes()]

    pos = nx.spring_layout(graph, seed=seed, k=None)
    fig, ax = plt.subplots(figsize=(12, 12))
    nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.25, width=0.5,
                           edge_color="#999999")
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors,
                           node_size=40, linewidths=0.3, edgecolors="white")
    ax.set_title(
        f"SSN: {graph.number_of_nodes()} proteins, "
        f"{result.n_clusters} clusters",
        fontsize=14,
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  static network   -> {path}", file=sys.stderr)
    return True


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
    net = Network(height="800px", width="100%", bgcolor="#ffffff",
                  font_color="#222222", notebook=False, cdn_resources="in_line")
    net.barnes_hut(gravity=-8000, spring_length=120)

    for node in graph.nodes():
        cluster = result.clusters[node]
        net.add_node(
            node, label=node, title=f"{node}\n{cluster}",
            color=colors[cluster], size=12,
        )
    for u, v, data in graph.edges(data=True):
        net.add_edge(u, v, value=float(data["pident"]),
                     title=f"{data['pident']:.1f}% id")

    net.save_graph(path)
    print(f"  interactive html  -> {path}", file=sys.stderr)
    return True


def sweep_plot(rows: list[dict], path: str, dpi: int = 200) -> None:
    """Plot threshold-sweep results: #clusters and largest cluster vs identity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ident = [r["identity"] for r in rows]
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(ident, [r["n_clusters"] for r in rows], "-o", color="#4C72B0",
             label="clusters")
    ax1.set_xlabel("identity threshold (%)")
    ax1.set_ylabel("number of clusters", color="#4C72B0")
    ax1.tick_params(axis="y", labelcolor="#4C72B0")

    ax2 = ax1.twinx()
    ax2.plot(ident, [r["largest_cluster"] for r in rows], "-s", color="#C44E52",
             label="largest cluster")
    ax2.set_ylabel("size of largest cluster", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")

    ax1.set_title("Threshold sweep")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  sweep plot -> {path}", file=sys.stderr)
