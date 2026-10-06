#!/usr/bin/env python3
"""Build ROBUST family-specific k-mer markers with KMC.

For each ATHILA family:
  * "in"  DB  = k-mers present in >= --prevalence of the family's PURE (non-mosaic)
               elements (KMC -ci = prevalence * n_members). Families with fewer than
               --min-members pure elements fall back to their exemplar k-mers.
  * "out" DB  = union of every OTHER family's k-mers seen in >= 2 pure elements,
               plus all other families' exemplar k-mers.
  * markers   = in - out   (kmc_tools kmers_subtract).

Writes markers.tsv (canonical k-mer, family) for kmer_recomb.py --markers.

Usage:
    kmc_markers.py --primary primary_family.tsv --mosaics kmer_recombinants.tsv \
        --fulllength DIR --exemplars EX.fasta --out-dir OUT --kmc-bin DIR
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import subprocess
from collections import defaultdict

import lib_ids as L
from athila_recomb_nt import read_fasta, family_of


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def kmc_count(kmc, fasta, prefix, tmp, k, ci, threads):
    run([kmc, f"-k{k}", f"-ci{ci}", "-cs1000000", f"-t{threads}", "-fm",
         fasta, prefix, tmp])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary", required=True)
    ap.add_argument("--mosaics", required=True)
    ap.add_argument("--fulllength", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--kmc-bin", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--prevalence", type=float, default=0.25)
    ap.add_argument("--min-members", type=int, default=30)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    kmc = os.path.join(args.kmc_bin, "kmc")
    tools = os.path.join(args.kmc_bin, "kmc_tools")
    work = os.path.join(args.out_dir, "kmc_work")
    tmp = os.path.join(work, "tmp")
    os.makedirs(tmp, exist_ok=True)

    mosaic = {l.split("\t")[0] for l in open(args.mosaics).readlines()[1:]}
    members = defaultdict(set)
    for l in open(args.primary).readlines()[1:]:
        ns, fam = l.rstrip("\n").split("\t")
        if ns not in mosaic:
            members[fam].add(ns)

    # write per-family FASTAs of pure elements
    fam_fa = {f: open(os.path.join(work, f"{f}.pure.fa"), "w") for f in members}
    for flf in glob.glob(os.path.join(args.fulllength, "*DP_FULLLENGTH_RENAMED.fasta")):
        acc = L.accession_from_path(flf)
        for el, s in read_fasta(flf).items():
            ns = L.nsid(acc, el)
            for f, m in members.items():
                if ns in m:
                    fam_fa[f].write(f">{ns}\n{s}\n")
                    break
    for fh in fam_fa.values():
        fh.close()

    # exemplar FASTAs per family
    exe = defaultdict(list)
    for h, s in read_fasta(args.exemplars).items():
        if "_I" in h and "_LTR" not in h:
            exe[family_of(h)].append((h, s))
    families = sorted(set(exe) | set(members))
    for f in families:
        with open(os.path.join(work, f"{f}.exe.fa"), "w") as o:
            for h, s in exe.get(f, []):
                o.write(f">{h}\n{s}\n")

    # "in" and "presence" DBs
    source = {}
    for f in families:
        n = len(members.get(f, ()))
        if n >= args.min_members:
            ci = max(2, math.ceil(args.prevalence * n))
            kmc_count(kmc, os.path.join(work, f"{f}.pure.fa"),
                      os.path.join(work, f"{f}.in"), tmp, args.k, ci, args.threads)
            kmc_count(kmc, os.path.join(work, f"{f}.pure.fa"),
                      os.path.join(work, f"{f}.pres"), tmp, args.k, 2, args.threads)
            source[f] = f"pure n={n}, ci={ci}"
        else:
            kmc_count(kmc, os.path.join(work, f"{f}.exe.fa"),
                      os.path.join(work, f"{f}.in"), tmp, args.k, 1, args.threads)
            source[f] = f"exemplar (pure n={n})"
        kmc_count(kmc, os.path.join(work, f"{f}.exe.fa"),
                  os.path.join(work, f"{f}.exeall"), tmp, args.k, 1, args.threads)

    # markers = in(F) - union(other families' presence + exemplars)
    out_tsv = os.path.join(args.out_dir, "robust_markers.tsv")
    with open(out_tsv, "w") as out:
        out.write("kmer\tfamily\n")
        for f in families:
            others = []
            for g in families:
                if g == f:
                    continue
                others.append(os.path.join(work, f"{g}.exeall"))
                if os.path.exists(os.path.join(work, f"{g}.pres.kmc_pre")):
                    others.append(os.path.join(work, f"{g}.pres"))
            acc_db = others[0]
            for i, o in enumerate(others[1:]):
                nxt = os.path.join(work, f"_u_{f}_{i}")
                run([tools, "simple", acc_db, o, "union", nxt])
                acc_db = nxt
            mk = os.path.join(work, f"{f}.markers")
            run([tools, "simple", os.path.join(work, f"{f}.in"), acc_db,
                 "kmers_subtract", mk])
            dump = mk + ".txt"
            run([tools, "transform", mk, "dump", dump])
            n = 0
            for line in open(dump):
                out.write(f"{line.split()[0]}\t{f}\n")
                n += 1
            print(f"[kmc] {f:>10}: {n:>6} robust markers  ({source[f]})")
            for p in glob.glob(os.path.join(work, f"_u_{f}_*")):
                os.remove(p)
    print(f"[kmc] wrote {out_tsv}")


if __name__ == "__main__":
    main()
