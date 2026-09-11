"""Command line interface: `ssn build | plot | stats | sweep`."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__, diamond, network, visualize


def _print_summary(result: network.NetworkResult) -> None:
    print("\nNetwork summary", file=sys.stderr)
    print(f"  proteins (nodes) : {result.n_nodes}", file=sys.stderr)
    print(f"  edges            : {result.n_edges}", file=sys.stderr)
    print(f"  clusters         : {result.n_clusters}", file=sys.stderr)
    print(f"  singletons       : {len(result.singletons)}", file=sys.stderr)
    print(f"  self-hits dropped: {result.n_self_hits}", file=sys.stderr)
    print(f"  reciprocal dups  : {result.n_reciprocal}", file=sys.stderr)
    print(f"  below threshold  : {result.n_below_threshold}", file=sys.stderr)
    top = network.cluster_sizes(result)[:5]
    if top:
        shown = ", ".join(f"{label} ({size})" for label, size in top)
        print(f"  biggest clusters : {shown}", file=sys.stderr)


def _visualise(result, outdir: Path, stem: str, args) -> None:
    if args.no_plot:
        return
    want_static = args.viz in ("both", "static")
    want_html = args.viz in ("both", "interactive")
    if want_static:
        visualize.static_plot(result, str(outdir / f"{stem}.network.png"),
                              top_n=args.top_clusters, max_nodes=args.max_nodes)
    if want_html:
        visualize.interactive_plot(result, str(outdir / f"{stem}.network.html"),
                                   top_n=args.top_clusters)


# --------------------------------------------------------------------------- #
# ssn build
# --------------------------------------------------------------------------- #
def cmd_build(args: argparse.Namespace) -> None:
    fasta = args.fasta
    if not Path(fasta).is_file():
        sys.exit(f"error: input FASTA not found: {fasta}")

    stem = Path(fasta).stem
    outdir = Path(args.outdir or f"{stem}_ssn")
    outdir.mkdir(parents=True, exist_ok=True)

    n_seqs = diamond.count_sequences(fasta)
    print(f"proteinssn {__version__}", file=sys.stderr)
    print(f"input   : {fasta} ({n_seqs} proteins)", file=sys.stderr)
    print(f"outdir  : {outdir}", file=sys.stderr)
    print(f"filters : identity>={args.identity}%"
          + (f", coverage>={args.coverage}%" if args.coverage else "")
          + (f", evalue<={args.evalue}" if args.evalue else ""),
          file=sys.stderr)

    if args.table:
        table = Path(args.table)
        print(f"reusing existing DIAMOND table: {table}", file=sys.stderr)
    else:
        table = diamond.all_vs_all(
            fasta, outdir, threads=args.threads, sensitivity=args.sensitivity,
            max_target_seqs=args.max_target_seqs, evalue=args.diamond_evalue,
            iterate=not args.no_iterate, quiet=args.quiet,
        )

    result = network.build_graph(
        str(table), identity=args.identity, coverage=args.coverage,
        evalue=args.evalue, extra_nodes=diamond.fasta_ids(fasta),
    )

    edges_path = outdir / f"{stem}.edges.tsv"
    clusters_path = outdir / f"{stem}.clusters.tsv"
    network.write_edges(result, str(edges_path))
    network.write_clusters(result, str(clusters_path))
    print(f"\n  edge list     -> {edges_path}", file=sys.stderr)
    print(f"  cluster table -> {clusters_path}", file=sys.stderr)

    _visualise(result, outdir, stem, args)
    _print_summary(result)
    print("\nDone. Import the edge list + cluster table into Cytoscape, "
          "or open the .html in a browser.", file=sys.stderr)


# --------------------------------------------------------------------------- #
# ssn plot  (re-plot an existing edge list)
# --------------------------------------------------------------------------- #
def cmd_plot(args: argparse.Namespace) -> None:
    edges = Path(args.edges)
    if not edges.is_file():
        sys.exit(f"error: edge list not found: {edges}")

    import networkx as nx
    graph = nx.Graph()
    with open(edges) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        # Tolerate files with or without a header row.
        if header and header[0].lower() not in ("source", "node1", "#name"):
            handle.seek(0)
            reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            pident = float(row[2]) if len(row) > 2 else 100.0
            graph.add_edge(row[0], row[1], pident=pident)

    clusters, singletons = network._label_clusters(graph)
    result = network.NetworkResult(
        graph=graph, clusters=clusters,
        n_edges=graph.number_of_edges(), singletons=singletons,
    )

    stem = edges.stem.replace(".edges", "")
    outdir = Path(args.outdir or edges.parent)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.viz in ("both", "static"):
        visualize.static_plot(result, str(outdir / f"{stem}.network.png"),
                              top_n=args.top_clusters, max_nodes=args.max_nodes)
    if args.viz in ("both", "interactive"):
        visualize.interactive_plot(result, str(outdir / f"{stem}.network.html"),
                                   top_n=args.top_clusters)
    _print_summary(result)


# --------------------------------------------------------------------------- #
# ssn stats  (cluster size distribution from a cluster table)
# --------------------------------------------------------------------------- #
def cmd_stats(args: argparse.Namespace) -> None:
    path = Path(args.clusters)
    if not path.is_file():
        sys.exit(f"error: cluster table not found: {path}")

    counts: dict[str, int] = {}
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header and header[0].lower() not in ("name", "#name", "node"):
            handle.seek(0)
            reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            counts[row[1]] = counts.get(row[1], 0) + 1

    if not counts:
        sys.exit("error: no clusters found in table.")

    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    sizes = [c for _, c in ranked]
    total = sum(sizes)
    singletons = sum(1 for c in sizes if c == 1)

    print(f"proteins : {total}")
    print(f"clusters : {len(ranked)}")
    print(f"singletons : {singletons}")
    print(f"largest  : {sizes[0]}  ({ranked[0][0]})")
    print(f"median cluster size : {sizes[len(sizes)//2]}")
    print("\ncluster\tsize")
    for label, size in ranked[: args.top]:
        print(f"{label}\t{size}")
    if len(ranked) > args.top:
        print(f"... ({len(ranked) - args.top} more; use --top to show more)")


# --------------------------------------------------------------------------- #
# ssn sweep  (scan identity thresholds to help pick one)
# --------------------------------------------------------------------------- #
def cmd_sweep(args: argparse.Namespace) -> None:
    if args.table:
        table = Path(args.table)
        if not table.is_file():
            sys.exit(f"error: DIAMOND table not found: {table}")
        fasta_ids = None
        stem = table.stem.replace(".diamond", "")
        outdir = Path(args.outdir or table.parent)
        outdir.mkdir(parents=True, exist_ok=True)
    else:
        fasta = args.fasta
        if not fasta or not Path(fasta).is_file():
            sys.exit("error: provide a FASTA (to run DIAMOND) or --table.")
        stem = Path(fasta).stem
        outdir = Path(args.outdir or f"{stem}_ssn")
        outdir.mkdir(parents=True, exist_ok=True)
        table = diamond.all_vs_all(
            fasta, outdir, threads=args.threads, sensitivity=args.sensitivity,
            quiet=args.quiet,
        )
        fasta_ids = diamond.fasta_ids(fasta)

    thresholds = _frange(args.start, args.stop, args.step)
    rows = []
    print("identity\tclusters\tedges\tlargest\tsingletons", file=sys.stderr)
    for ident in thresholds:
        result = network.build_graph(
            str(table), identity=ident, coverage=args.coverage,
            evalue=args.evalue, extra_nodes=fasta_ids,
        )
        sizes = network.cluster_sizes(result)
        largest = sizes[0][1] if sizes else 0
        rows.append({
            "identity": ident, "n_clusters": result.n_clusters,
            "n_edges": result.n_edges, "largest_cluster": largest,
            "singletons": len(result.singletons),
        })
        print(f"{ident:g}\t{result.n_clusters}\t{result.n_edges}\t"
              f"{largest}\t{len(result.singletons)}", file=sys.stderr)

    tsv_path = outdir / f"{stem}.sweep.tsv"
    with open(tsv_path, "w") as out:
        out.write("identity\tn_clusters\tn_edges\tlargest_cluster\tsingletons\n")
        for r in rows:
            out.write(f"{r['identity']:g}\t{r['n_clusters']}\t{r['n_edges']}\t"
                      f"{r['largest_cluster']}\t{r['singletons']}\n")
    print(f"\n  sweep table -> {tsv_path}", file=sys.stderr)
    if not args.no_plot:
        visualize.sweep_plot(rows, str(outdir / f"{stem}.sweep.png"))


def _frange(start: float, stop: float, step: float) -> list[float]:
    vals, v = [], start
    while v <= stop + 1e-9:
        vals.append(round(v, 6))
        v += step
    return vals


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #
def _add_viz_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--viz", choices=["both", "static", "interactive"],
                   default="both", help="which plots to make (default: both)")
    p.add_argument("--top-clusters", type=int, default=14, metavar="N",
                   help="colour the N biggest clusters; rest grey (default: 14)")
    p.add_argument("--max-nodes", type=int, default=3000, metavar="N",
                   help="skip the static plot above this many nodes "
                        "(default: 3000)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ssn",
        description="Build and visualise protein Sequence Similarity Networks.",
    )
    parser.add_argument("--version", action="version",
                        version=f"proteinssn {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # build
    b = sub.add_parser("build", help="FASTA -> network + clusters + plots")
    b.add_argument("fasta", help="protein FASTA file")
    b.add_argument("-i", "--identity", type=float, default=30.0,
                   help="min %% identity for an edge (default: 30)")
    b.add_argument("-c", "--coverage", type=float, default=None,
                   help="min %% coverage of BOTH sequences (default: off)")
    b.add_argument("-e", "--evalue", type=float, default=None,
                   help="max e-value for an edge (default: off)")
    b.add_argument("-o", "--outdir", default=None,
                   help="output directory (default: <name>_ssn)")
    b.add_argument("-t", "--threads", type=int, default=4,
                   help="DIAMOND threads (default: 4)")
    b.add_argument("--sensitivity", default="sensitive",
                   choices=list(diamond.SENSITIVITY_FLAGS),
                   help="DIAMOND sensitivity (default: sensitive)")
    b.add_argument("--max-target-seqs", type=int, default=3000,
                   help="DIAMOND --max-target-seqs (default: 3000)")
    b.add_argument("--diamond-evalue", type=float, default=1e-5,
                   help="e-value passed to DIAMOND search (default: 1e-5)")
    b.add_argument("--no-iterate", action="store_true",
                   help="disable DIAMOND --iterate")
    b.add_argument("--table", default=None,
                   help="reuse an existing DIAMOND table, skip the search")
    b.add_argument("--no-plot", action="store_true", help="skip plotting")
    b.add_argument("--quiet", action="store_true", help="less DIAMOND output")
    _add_viz_opts(b)
    b.set_defaults(func=cmd_build)

    # plot
    p = sub.add_parser("plot", help="re-plot an existing edge list")
    p.add_argument("edges", help="edge list (source<TAB>target[<TAB>pident])")
    p.add_argument("-o", "--outdir", default=None)
    _add_viz_opts(p)
    p.set_defaults(func=cmd_plot)

    # stats
    s = sub.add_parser("stats", help="cluster size distribution from a table")
    s.add_argument("clusters", help="cluster table (name<TAB>cluster)")
    s.add_argument("--top", type=int, default=20,
                   help="rows to show (default: 20)")
    s.set_defaults(func=cmd_stats)

    # sweep
    w = sub.add_parser("sweep", help="scan identity thresholds to pick one")
    w.add_argument("fasta", nargs="?", default=None,
                   help="protein FASTA (omit if using --table)")
    w.add_argument("--table", default=None,
                   help="reuse an existing DIAMOND table instead of a FASTA")
    w.add_argument("--start", type=float, default=20.0,
                   help="first identity threshold (default: 20)")
    w.add_argument("--stop", type=float, default=90.0,
                   help="last identity threshold (default: 90)")
    w.add_argument("--step", type=float, default=10.0,
                   help="step between thresholds (default: 10)")
    w.add_argument("-c", "--coverage", type=float, default=None)
    w.add_argument("-e", "--evalue", type=float, default=None)
    w.add_argument("-o", "--outdir", default=None)
    w.add_argument("-t", "--threads", type=int, default=4)
    w.add_argument("--sensitivity", default="sensitive",
                   choices=list(diamond.SENSITIVITY_FLAGS))
    w.add_argument("--no-plot", action="store_true")
    w.add_argument("--quiet", action="store_true")
    w.set_defaults(func=cmd_sweep)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
