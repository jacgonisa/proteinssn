#!/usr/bin/env python3
"""Phase 0 — inventory. Walk every genome's cls.tsv and emit a tidy element table
plus a human-readable summary. Read-only on the inputs.

Usage:
    p0_inventory.py --input-dir DIR --out-dir OUT [--counts Athila_counts.tsv]
"""

from __future__ import annotations

import argparse
import glob
import os
from collections import Counter, defaultdict

import lib_ids as L


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--counts", default=None,
                    help="optional Athila_counts.tsv for cross-reference")
    args = ap.parse_args()

    inv_dir = os.path.join(args.out_dir, "inventory")
    os.makedirs(inv_dir, exist_ok=True)

    cls_files = sorted(glob.glob(os.path.join(args.input_dir, "*.cls.tsv")))
    if not cls_files:
        raise SystemExit(f"no *.cls.tsv under {args.input_dir}")

    clade_counts = Counter()
    complete_counts = Counter()
    domain_counts = Counter()          # elements carrying each domain
    ncore_hist = Counter()             # how many of the 5 core domains present
    per_acc = Counter()
    n_full5 = n_has_pol = n_total = 0

    elements_tsv = os.path.join(inv_dir, "elements.tsv")
    with open(elements_tsv, "w") as out:
        out.write("ns_id\taccession\telement\tchrom\tstart\tend\ttag\t"
                  "clade\tcomplete\tstrand\tn_core_domains\tdomains\t"
                  "is_full5\thas_pol\tis_athila\n")
        for cf in cls_files:
            acc = L.accession_from_path(cf)
            for rec in L.parse_cls_tsv(cf):
                n_total += 1
                per_acc[acc] += 1
                clade_counts[rec["clade"]] += 1
                complete_counts[rec["complete"]] += 1

                present = set(rec["domains"])
                for d in present:
                    domain_counts[d] += 1
                core_present = [d for d in L.CORE_DOMAINS if d in present]
                ncore_hist[len(core_present)] += 1

                is_full5 = set(L.CORE_DOMAINS) <= present
                has_pol = any(d in present for d in L.POL_DOMAINS)
                is_ath = L.is_athila(rec)
                n_full5 += is_full5
                n_has_pol += has_pol

                coords = L.element_coords(rec["element"]) or {}
                out.write("\t".join(str(x) for x in [
                    L.nsid(acc, rec["element"]), acc, rec["element"],
                    coords.get("chrom", ""), coords.get("start", ""),
                    coords.get("end", ""), coords.get("tag", ""),
                    rec["clade"], rec["complete"], rec["strand"],
                    len(core_present),
                    ";".join(d for d in L.CORE_DOMAINS if d in present)
                    + ("|" + ";".join(sorted(present - set(L.CORE_DOMAINS)))
                       if present - set(L.CORE_DOMAINS) else ""),
                    int(is_full5), int(has_pol), int(is_ath),
                ]) + "\n")

    # optional cross-reference
    xref = ""
    if args.counts and os.path.exists(args.counts):
        counts_accs = set()
        with open(args.counts) as fh:
            next(fh, None)
            for line in fh:
                counts_accs.add(line.split("\t")[0].strip())
        file_accs = set(per_acc)
        xref = (f"- Athila_counts.tsv accessions: {len(counts_accs)} "
                f"(naming differs from filenames; filenames are authoritative)\n"
                f"  - in files not in counts: {len(file_accs - counts_accs)}\n"
                f"  - in counts not in files: {len(counts_accs - file_accs)}\n")

    summary = os.path.join(inv_dir, "summary.md")
    with open(summary, "w") as md:
        md.write("# Phase 0 inventory\n\n")
        md.write(f"- genomes (cls.tsv files): **{len(cls_files)}**\n")
        md.write(f"- total elements: **{n_total}**\n")
        md.write(f"- full-5 (GAG+PROT+RT+RH+INT): **{n_full5}**\n")
        md.write(f"- has >=1 pol domain (RT/RH/INT): **{n_has_pol}**\n\n")
        md.write(_table("Clade", clade_counts))
        md.write(_table("Complete", complete_counts))
        md.write(_table("Per-domain element counts", domain_counts,
                        order=L.CORE_DOMAINS))
        md.write("## Core-domain count per element\n\n")
        md.write("| #core domains | elements |\n|---:|---:|\n")
        for k in sorted(ncore_hist):
            md.write(f"| {k} | {ncore_hist[k]} |\n")
        md.write("\n")
        md.write("## Elements per accession (min / median / max)\n\n")
        vals = sorted(per_acc.values())
        med = vals[len(vals) // 2] if vals else 0
        md.write(f"- min {vals[0]}, median {med}, max {vals[-1]} "
                 f"across {len(per_acc)} accessions\n\n")
        if xref:
            md.write("## Cross-reference\n\n" + xref + "\n")

    print(f"[p0] {n_total} elements, {n_full5} full-5, "
          f"{len(cls_files)} genomes")
    print(f"[p0] wrote {elements_tsv}")
    print(f"[p0] wrote {summary}")


def _table(title: str, counter: Counter, order=None) -> str:
    keys = list(order) if order else []
    keys += [k for k, _ in counter.most_common() if k not in keys]
    rows = "".join(f"| {k} | {counter.get(k, 0)} |\n" for k in keys
                   if counter.get(k, 0))
    return f"## {title}\n\n| value | count |\n|---|---:|\n{rows}\n"


if __name__ == "__main__":
    main()
