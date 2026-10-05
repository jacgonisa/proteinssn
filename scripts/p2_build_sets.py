#!/usr/bin/env python3
"""Phase 2 — per-layer sequence sets.

Steps (see plan):
  1. Namespace every id by accession (mandatory; element ids collide across genomes).
  2. Keep LTR/Gypsy/Athila elements only.
  3. Per (element, domain) keep the LONGEST peptide (dom.faa can list a domain twice).
  4. Length filter: drop a domain peptide < `length_frac_min` x that domain's median.
  5. Redundancy collapse (Track A): mmseqs easy-cluster at 95% id / 80% cov on the
     CONCATENATED per-element peptide (not per domain), to pick representatives.
  6. Emit node_table.tsv, per-layer FASTAs (rep + all), and analysis_nodes.tsv
     (the primary node set = full-5 representatives).

Usage:
    p2_build_sets.py --input-dir DIR --out-dir OUT \
        [--length-frac-min 0.60] [--min-id 0.95] [--min-cov 0.80] [--threads 16] \
        [--mmseqs mmseqs]
"""

from __future__ import annotations

import argparse
import glob
import os
import statistics
import subprocess
import sys
from collections import defaultdict

import lib_ids as L


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--length-frac-min", type=float, default=0.60)
    ap.add_argument("--min-id", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--cov-mode", type=int, default=1)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--mmseqs", default="mmseqs")
    args = ap.parse_args()

    p2 = os.path.join(args.out_dir, "phase2")
    layers_dir = os.path.join(p2, "layers")
    os.makedirs(layers_dir, exist_ok=True)

    # ---- 1-2. node records from cls.tsv (Athila only) --------------------- #
    nodes: dict[str, dict] = {}   # ns_id -> record
    for cf in sorted(glob.glob(os.path.join(args.input_dir, "*.cls.tsv"))):
        acc = L.accession_from_path(cf)
        for rec in L.parse_cls_tsv(cf):
            if not L.is_athila(rec):
                continue
            ns = L.nsid(acc, rec["element"])
            present = set(rec["domains"])
            nodes[ns] = {
                "accession": acc, "element": rec["element"],
                "complete": rec["complete"], "strand": rec["strand"],
                "domains": present,
                "is_full5": set(L.CORE_DOMAINS) <= present,
            }
    print(f"[p2] Athila elements: {len(nodes)}")

    # ---- 3. longest peptide per (ns_id, domain) --------------------------- #
    # seqs[domain][ns_id] = (length, sequence)
    seqs: dict[str, dict[str, tuple[int, str]]] = {d: {} for d in L.CORE_DOMAINS}
    for ff in sorted(glob.glob(os.path.join(args.input_dir, "*.dom.faa"))):
        acc = L.accession_from_path(ff)
        for element, domain, _attrs, seq in L.parse_dom_faa(ff):
            if domain not in L.CORE_DOMAINS:
                continue
            ns = L.nsid(acc, element)
            if ns not in nodes:         # non-Athila element
                continue
            cur = seqs[domain].get(ns)
            if cur is None or len(seq) > cur[0]:
                seqs[domain][ns] = (len(seq), seq)

    # ---- 4. length filter per domain -------------------------------------- #
    kept: dict[str, dict[str, str]] = {d: {} for d in L.CORE_DOMAINS}
    for d in L.CORE_DOMAINS:
        lengths = [l for l, _ in seqs[d].values()]
        if not lengths:
            continue
        med = statistics.median(lengths)
        cut = args.length_frac_min * med
        for ns, (l, s) in seqs[d].items():
            if l >= cut:
                kept[d][ns] = s
        print(f"[p2] {d}: {len(seqs[d])} peptides, median len {med:.0f}, "
              f">= {cut:.0f} kept {len(kept[d])}")

    # ---- 5. concatenated peptide + mmseqs clustering ---------------------- #
    concat = os.path.join(p2, "concat.faa")
    n_concat = 0
    with open(concat, "w") as out:
        for ns in nodes:
            parts = [kept[d][ns] for d in L.CORE_DOMAINS if ns in kept[d]]
            if not parts:
                continue
            out.write(f">{ns}\n{''.join(parts)}\n")
            n_concat += 1
    print(f"[p2] concatenated peptides written: {n_concat}")

    rep_of = _mmseqs_cluster(concat, p2, args)
    reps = set(rep_of.values())
    print(f"[p2] clusters/representatives: {len(reps)} "
          f"(from {len(rep_of)} clustered elements)")

    # ---- 6a. per-layer FASTAs (all copies + representatives) --------------- #
    for d in L.CORE_DOMAINS:
        with open(os.path.join(layers_dir, f"{d}.all.faa"), "w") as fa_all, \
             open(os.path.join(layers_dir, f"{d}.rep.faa"), "w") as fa_rep:
            for ns, s in kept[d].items():
                fa_all.write(f">{ns}\n{s}\n")
                if rep_of.get(ns) == ns:    # this element is its cluster rep
                    fa_rep.write(f">{ns}\n{s}\n")

    # ---- 6b. node_table + analysis_nodes ---------------------------------- #
    node_table = os.path.join(p2, "node_table.tsv")
    with open(node_table, "w") as out:
        out.write("ns_id\taccession\telement\tcomplete\tstrand\t"
                  "domains_kept\tis_full5\tcluster_rep\tis_rep\n")
        for ns, r in nodes.items():
            dk = [d for d in L.CORE_DOMAINS if ns in kept[d]]
            rep = rep_of.get(ns, "")
            out.write("\t".join(str(x) for x in [
                ns, r["accession"], r["element"], r["complete"], r["strand"],
                ";".join(dk), int(r["is_full5"]), rep, int(rep == ns),
            ]) + "\n")

    # primary node set: full-5 AND representative AND all 5 domains survived filter
    analysis = os.path.join(p2, "analysis_nodes.tsv")
    primary: list[str] = []
    for ns, r in nodes.items():
        if (r["is_full5"] and rep_of.get(ns) == ns
                and all(ns in kept[d] for d in L.CORE_DOMAINS)):
            primary.append(ns)
    with open(analysis, "w") as out:
        out.write("ns_id\n")
        for ns in primary:
            out.write(ns + "\n")
    n_primary = len(primary)

    # per-layer FASTAs restricted to the primary node set (same nodes in all 5
    # layers -> the multiplex Phase 4 consumes directly).
    primary_set = set(primary)
    for d in L.CORE_DOMAINS:
        with open(os.path.join(layers_dir, f"{d}.primary.faa"), "w") as fa:
            for ns in primary:
                if ns in kept[d]:
                    fa.write(f">{ns}\n{kept[d][ns]}\n")

    print(f"[p2] node_table: {len(nodes)} rows -> {node_table}")
    print(f"[p2] PRIMARY node set (full-5 representatives): {n_primary} -> {analysis}")


def _mmseqs_cluster(concat_faa: str, outdir: str, args) -> dict[str, str]:
    """Run mmseqs easy-cluster; return element -> representative map."""
    tmp = os.path.join(outdir, "mmseqs_tmp")
    prefix = os.path.join(outdir, "clu")
    os.makedirs(tmp, exist_ok=True)
    cmd = [
        args.mmseqs, "easy-cluster", concat_faa, prefix, tmp,
        "--min-seq-id", str(args.min_id), "-c", str(args.min_cov),
        "--cov-mode", str(args.cov_mode), "--threads", str(args.threads),
        "-v", "1",
    ]
    print("  $ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    rep_of: dict[str, str] = {}
    with open(prefix + "_cluster.tsv") as fh:   # cols: representative  member
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            rep_of[member] = rep
    return rep_of


if __name__ == "__main__":
    main()
