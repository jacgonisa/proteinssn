"""Turn an all-vs-all DIAMOND table into a filtered network and clusters.

DIAMOND is run with output format 6 and these 11 columns:

    qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen

An edge between two proteins is kept when the alignment passes the identity,
(optional) coverage and (optional) e-value thresholds. Reciprocal hits (A-B and
B-A) and self-hits are collapsed. Clusters are the connected components of the
resulting graph.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx

# Column order that `ssn build` asks DIAMOND to emit.
DIAMOND_COLUMNS = [
    "qseqid", "sseqid", "evalue", "pident", "bitscore",
    "qstart", "qend", "qlen", "sstart", "send", "slen",
]


@dataclass
class NetworkResult:
    graph: nx.Graph
    clusters: dict[str, str]  # protein id -> "cluster01"
    n_edges: int = 0
    n_self_hits: int = 0
    n_reciprocal: int = 0
    n_below_threshold: int = 0
    singletons: list[str] = field(default_factory=list)

    @property
    def n_nodes(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def n_clusters(self) -> int:
        return len({c for c in self.clusters.values()})


def _coverage(start: str, end: str, length: str) -> float:
    """Percent of a sequence spanned by the alignment."""
    length = int(length)
    if length == 0:
        return 0.0
    return 100.0 * (int(end) - int(start)) / length


def build_graph(
    diamond_table: str,
    identity: float,
    coverage: float | None = None,
    evalue: float | None = None,
    extra_nodes: list[str] | None = None,
) -> NetworkResult:
    """Read a raw DIAMOND table and return the filtered graph + clusters.

    Streams the file line by line so it stays cheap on large all-vs-all tables.
    Pass ``extra_nodes`` (e.g. every protein id in the FASTA) to make sure
    proteins with no surviving edges still show up as singleton clusters.
    """
    graph = nx.Graph()
    seen: set[tuple[str, str]] = set()
    n_self = n_recip = n_below = 0

    with open(diamond_table) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            q, s, ev, pident = parts[0], parts[1], parts[2], parts[3]

            if q == s:
                n_self += 1
                continue

            key = (q, s) if q < s else (s, q)
            if key in seen:
                n_recip += 1
                continue
            seen.add(key)

            if float(pident) < identity:
                n_below += 1
                continue

            if coverage is not None:
                cov_q = _coverage(parts[5], parts[6], parts[7])
                cov_s = _coverage(parts[8], parts[9], parts[10])
                if cov_q < coverage or cov_s < coverage:
                    n_below += 1
                    continue

            if evalue is not None and float(ev) > evalue:
                n_below += 1
                continue

            graph.add_edge(q, s, pident=float(pident), evalue=float(ev))

    if extra_nodes:
        graph.add_nodes_from(extra_nodes)

    clusters, singletons = _label_clusters(graph)
    return NetworkResult(
        graph=graph,
        clusters=clusters,
        n_edges=graph.number_of_edges(),
        n_self_hits=n_self,
        n_reciprocal=n_recip,
        n_below_threshold=n_below,
        singletons=singletons,
    )


def _label_clusters(graph: nx.Graph) -> tuple[dict[str, str], list[str]]:
    """Label connected components cluster01, cluster02, ... largest first."""
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    width = max(2, len(str(len(components))))
    mapping: dict[str, str] = {}
    singletons: list[str] = []
    for i, component in enumerate(components, start=1):
        label = f"cluster{str(i).zfill(width)}"
        for node in component:
            mapping[node] = label
        if len(component) == 1:
            singletons.extend(component)
    return mapping, singletons


def write_edges(result: NetworkResult, path: str) -> None:
    """Write a tab-separated edge list: source, target, pident."""
    with open(path, "w") as out:
        out.write("source\ttarget\tpident\n")
        for u, v, data in result.graph.edges(data=True):
            out.write(f"{u}\t{v}\t{data['pident']:.1f}\n")


def write_clusters(result: NetworkResult, path: str) -> None:
    """Write a node->cluster table (Cytoscape-friendly attribute file)."""
    # Order by cluster label then id so the file is stable and easy to scan.
    rows = sorted(result.clusters.items(), key=lambda kv: (kv[1], kv[0]))
    with open(path, "w") as out:
        out.write("name\tcluster\n")
        for name, cluster in rows:
            out.write(f"{name}\t{cluster}\n")


def cluster_sizes(result: NetworkResult) -> list[tuple[str, int]]:
    """Cluster label -> member count, largest first."""
    counts: dict[str, int] = {}
    for cluster in result.clusters.values():
        counts[cluster] = counts.get(cluster, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
