#!/usr/bin/env python3
"""Phase 3 — all-vs-all similarity for one layer, normalised to a bitscore ratio.

For one domain layer: DIAMOND all-vs-all, then every edge is normalised as

    ratio(i,j) = bitscore(i,j) / min(selfbit(i), selfbit(j))

so layers of different conservation (RT >> GAG) are comparable. Raw bitscores or a
fixed E-value are never compared across layers.

Reuses `proteinssn.diamond` (find_diamond + sensitivity flags) where available.

Usage:
    p3_allvsall.py --faa LAYER.faa --out-prefix OUT/GAG [--nodes analysis_nodes.tsv]
        [--sensitivity very-sensitive] [--max-target-seqs 5000] [--evalue 1e-3]
        [--threads 16]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

OUTFMT = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
          "qlen", "slen"]


def find_diamond() -> str:
    try:
        from proteinssn.diamond import find_diamond as _f
        return _f()
    except Exception:
        exe = shutil.which("diamond")
        if not exe:
            sys.exit("error: 'diamond' not found (conda install -c bioconda diamond)")
        return exe


def sensitivity_flags(name: str) -> list[str]:
    try:
        from proteinssn.diamond import SENSITIVITY_FLAGS
        return SENSITIVITY_FLAGS.get(name, ["--very-sensitive"])
    except Exception:
        return {["fast"][0]: []}.get(name, [f"--{name}"])


def read_nodes(path: str | None) -> set[str] | None:
    if not path:
        return None
    nodes = set()
    with open(path) as fh:
        header = next(fh, None)
        for line in fh:
            s = line.strip()
            if s:
                nodes.add(s.split("\t")[0])
    return nodes


def subset_fasta(faa: str, nodes: set[str] | None, out: str) -> int:
    n = 0
    keep = True
    with open(faa) as fh, open(out, "w") as o:
        for line in fh:
            if line.startswith(">"):
                sid = line[1:].strip().split()[0]
                keep = nodes is None or sid in nodes
                if keep:
                    n += 1
            if keep:
                o.write(line)
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--nodes", default=None)
    ap.add_argument("--sensitivity", default="very-sensitive")
    ap.add_argument("--max-target-seqs", type=int, default=5000)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    nodes = read_nodes(args.nodes)

    faa = args.faa
    if nodes is not None:
        faa = args.out_prefix + ".input.faa"
        n = subset_fasta(args.faa, nodes, faa)
        print(f"[p3] {os.path.basename(args.out_prefix)}: {n} sequences in node set")

    diamond = find_diamond()
    db = args.out_prefix + ".dmnd"
    raw = args.out_prefix + ".diamond.tsv"
    tmp = args.out_prefix + "_tmp"
    os.makedirs(tmp, exist_ok=True)

    subprocess.run([diamond, "makedb", "--in", faa, "--db", db, "--quiet"],
                   check=True)
    cmd = [diamond, "blastp", "-q", faa, "--db", db, "-o", raw,
           "--outfmt", "6", *OUTFMT,
           "--max-target-seqs", str(args.max_target_seqs),
           "--evalue", str(args.evalue), "--threads", str(args.threads),
           "--tmpdir", tmp, "--quiet", *sensitivity_flags(args.sensitivity)]
    print("  $ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

    # ---- self bitscores + best pair bitscore ------------------------------ #
    selfbit: dict[str, float] = {}
    pair_bit: dict[tuple[str, str], float] = {}
    pair_ev: dict[tuple[str, str], float] = {}
    pair_pid: dict[tuple[str, str], float] = {}
    with open(raw) as fh:
        for line in fh:
            q, s, pid, _ln, ev, bit, _ql, _sl = line.rstrip("\n").split("\t")
            bit = float(bit)
            if q == s:
                if bit > selfbit.get(q, 0):
                    selfbit[q] = bit
                continue
            key = (q, s) if q < s else (s, q)
            if bit > pair_bit.get(key, 0):
                pair_bit[key] = bit
                pair_ev[key] = float(ev)
                pair_pid[key] = float(pid)

    # fallback self-bitscore for nodes lacking a self-hit
    for key, bit in pair_bit.items():
        for node in key:
            if node not in selfbit or bit > selfbit[node]:
                # only as a floor; real self-hit (if present) already set above
                selfbit.setdefault(node, bit)

    edges = args.out_prefix + ".edges.tsv"
    n_edges = 0
    with open(edges, "w") as out:
        out.write("source\ttarget\tpident\tbitscore\tratio\tevalue\n")
        for (a, b), bit in pair_bit.items():
            denom = min(selfbit.get(a, bit), selfbit.get(b, bit)) or bit
            ratio = min(1.0, bit / denom)
            out.write(f"{a}\t{b}\t{pair_pid[(a,b)]:.1f}\t{bit:.1f}\t"
                      f"{ratio:.4f}\t{pair_ev[(a,b)]:.2e}\n")
            n_edges += 1

    shutil.rmtree(tmp, ignore_errors=True)
    print(f"[p3] {os.path.basename(args.out_prefix)}: {n_edges} edges -> {edges}")


if __name__ == "__main__":
    main()
