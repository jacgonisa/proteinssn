#!/usr/bin/env python3
"""Find every ATHILA-derived copy in each genome, including truncated fragments.

Per genome: minimap2 (-x asm20) maps each chromosome (query) against the TAIR12
ATHILA exemplars (target); every collinear hit is an ATHILA-derived segment. Hits
are merged into copies (gap <= --merge-gap), annotated with family / part (LTR vs
internal) by aligned bp, and labelled against Athilafinder's own calls:

  full_length   overlaps an Athilafinder intact element
  solo_ltr      overlaps an Athilafinder solo LTR
  fragment      everything else (truncated internals, LTR fragments, remnants)

Fragment sequences (+/- --flank bp of genomic context) are written for the junction
scan; the flank lets truncation breakpoints be inspected later.

Usage:
    genome_copies.py --genomes DIR --summary-dir DIR --exemplars EX.fasta \
        --out-dir OUT [--jobs 6] [--threads 4]
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import lib_ids as L


def ref_meta(name):
    base = name.split("#")[0]
    m = re.match(r"(ATHILA\d+[abc]?)", base)
    return (m.group(1) if m else base), ("LTR" if "_LTR" in base else "I")


def athilafinder_calls(summary):
    calls = defaultdict(list)
    for line in open(summary):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7 or f[0] == "chr":
            continue
        try:
            calls[f[0]].append((int(f[1]), int(f[2]), f[5], f[3]))
        except ValueError:
            pass
    return calls


def read_genome(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def run_one(job):
    gpath, summary, exemplars, out, threads, merge_gap, min_len, min_ident, flank, mm2 = job
    acc = L.accession_from_path(gpath)
    paf = os.path.join(out, "paf", f"{acc}.paf")
    os.makedirs(os.path.dirname(paf), exist_ok=True)
    if not os.path.exists(paf):
        with open(paf + ".tmp", "w") as fh:
            subprocess.run([mm2, "-x", "asm20", "-c", "--secondary=no", "-t", str(threads),
                            exemplars, gpath], stdout=fh, stderr=subprocess.DEVNULL,
                           check=True)
        os.rename(paf + ".tmp", paf)
    hits = defaultdict(list)
    for line in open(paf):
        f = line.split("\t")
        qs, qe, nm, al = int(f[2]), int(f[3]), int(f[9]), int(f[10])
        if qe - qs < min_len or not al or nm / al < min_ident:
            continue
        fam, part = ref_meta(f[5])
        hits[f[0]].append((qs, qe, fam, part, nm))
    calls = athilafinder_calls(summary)
    genome = None
    copies = []
    for chrom, hs in hits.items():
        hs.sort()
        cur = None
        for qs, qe, fam, part, nm in hs:
            if cur and qs - cur["end"] <= merge_gap:
                cur["end"] = max(cur["end"], qe)
                cur["bp"][(fam, part)] = cur["bp"].get((fam, part), 0) + (qe - qs)
                cur["n"] += 1
            else:
                if cur:
                    copies.append(cur)
                cur = {"chrom": chrom, "start": qs, "end": qe, "n": 1,
                       "bp": {(fam, part): qe - qs}}
        if cur:
            copies.append(cur)
    rows, frag_fa = [], []
    for c in copies:
        label, af_id = "fragment", ""
        for s, e, kind, tid in calls.get(c["chrom"], []):
            ov = min(e, c["end"]) - max(s, c["start"])
            if ov > 0.5 * min(e - s, c["end"] - c["start"]):
                label = "full_length" if kind == "intact" else "solo_ltr"
                af_id = tid
                break
        ibp = sum(v for (fm, p), v in c["bp"].items() if p == "I")
        lbp = sum(v for (fm, p), v in c["bp"].items() if p == "LTR")
        fam_i = max(((fm, v) for (fm, p), v in c["bp"].items() if p == "I"),
                    key=lambda x: x[1], default=("", 0))[0]
        fam_l = max(((fm, v) for (fm, p), v in c["bp"].items() if p == "LTR"),
                    key=lambda x: x[1], default=("", 0))[0]
        if label == "fragment":
            label = ("fragment_internal" if ibp >= 300 else "fragment_ltr")
        cid = f"{acc}::{c['chrom']}:{c['start']}-{c['end']}"
        rows.append([cid, acc, c["chrom"], c["start"], c["end"], c["end"] - c["start"],
                     label, af_id, fam_i, fam_l, ibp, lbp, c["n"]])
        if label == "fragment_internal":
            frag_fa.append((cid, c["chrom"], c["start"], c["end"]))
    if frag_fa:
        genome = read_genome(gpath)
        with open(os.path.join(out, "fragments", f"{acc}.fa"), "w") as fh:
            for cid, chrom, s, e in frag_fa:
                seq = genome.get(chrom, "")
                a, b = max(0, s - flank), min(len(seq), e + flank)
                fh.write(f">{cid}|flank5={s - a}|flank3={b - e}\n{seq[a:b]}\n")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--minimap2", default=os.path.expanduser("~/minimap2/minimap2"))
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--merge-gap", type=int, default=500)
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--min-ident", type=float, default=0.70)
    ap.add_argument("--flank", type=int, default=50)
    args = ap.parse_args()
    out = os.path.join(args.out_dir, "genome_copies")
    os.makedirs(os.path.join(out, "fragments"), exist_ok=True)

    jobs = []
    for g in sorted(glob.glob(os.path.join(args.genomes, "*.fasta"))):
        acc = L.accession_from_path(g)
        summ = glob.glob(os.path.join(args.summary_dir, f"{acc}.fasta.DP_FULLLENGTH_AND_SOLO_SUMMARY_TABLE.txt"))
        if not summ:
            continue
        jobs.append((g, summ[0], args.exemplars, out, args.threads, args.merge_gap,
                     args.min_len, args.min_ident, args.flank, args.minimap2))
    print(f"[copies] genomes: {len(jobs)}", flush=True)
    rows = []
    with ProcessPoolExecutor(args.jobs) as ex:
        for i, r in enumerate(ex.map(run_one, jobs), 1):
            rows.extend(r)
            if i % 10 == 0:
                print(f"[copies] {i}/{len(jobs)} genomes done", flush=True)
    cols = ["copy_id", "accession", "chrom", "start", "end", "length", "label",
            "athilafinder_id", "internal_family", "ltr_family", "internal_bp", "ltr_bp",
            "n_hits"]
    with open(os.path.join(out, "copies.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    from collections import Counter
    print(f"[copies] {len(rows)} ATHILA-derived copies: "
          + ", ".join(f"{k} {v}" for k, v in Counter(r[6] for r in rows).most_common()))


if __name__ == "__main__":
    main()
