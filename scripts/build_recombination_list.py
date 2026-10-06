#!/usr/bin/env python3
"""Consolidate all candidate recombination events into one master table.

Two categories:
  * intra_athila   — element whose domains map to >=2 different ATHILA families
                     (nucleotide per-domain typing; athila_recomb_nt).
  * inter_clade    — Athila element carrying a domain from a DIFFERENT LTR
                     retrotransposon clade (Retand/Tork/Ivana/CRM/...), i.e.
                     retro x retro recombination (TEsorter per-domain clade call).

DNA-transposon domains inside Athila (e.g. EnSpm_CACTA) are reported separately as
nested insertions, NOT recombination.

Usage:
    build_recombination_list.py --data-dir DIR --nt-recomb recombinants_nt.tsv \
        --assign domain_family_nt.tsv --out-dir OUT
"""

from __future__ import annotations

import argparse
import glob
import os
from collections import Counter, defaultdict

import lib_ids as L

LTR_CLADES = {"Retand", "Tork", "Ivana", "CRM", "Tekay", "Reina", "Ale",
              "Galadriel", "Ogre", "Tcn1", "TatI", "TatII", "TatIII"}


def full5_typed(assign_path):
    """ns_ids with all four long domains confidently nt-typed (for full5 flag)."""
    fam = defaultdict(dict)
    with open(assign_path) as fh:
        next(fh)
        for line in fh:
            ns, dom, f, b, b2, mg, conf = line.rstrip("\n").split("\t")
            if int(conf):
                fam[ns][dom] = f
    return {ns for ns, dm in fam.items()
            if all(d in dm for d in ("GAG", "RT", "RH", "INT"))}, fam


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--nt-recomb", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    full5_set, fam = full5_typed(args.assign)

    # full-5 element membership (all 5 core domains present per cls)
    full5_present = set()
    cls_clade = {}
    for cf in sorted(glob.glob(os.path.join(args.data_dir, "*.cls.tsv"))):
        acc = L.accession_from_path(cf)
        for rec in L.parse_cls_tsv(cf):
            ns = L.nsid(acc, rec["element"])
            cls_clade[ns] = rec
            if L.is_athila(rec) and set(L.CORE_DOMAINS) <= set(rec["domains"]):
                full5_present.add(ns)

    rows = []

    # ---- intra-Athila (nucleotide) ---------------------------------------- #
    with open(args.nt_recomb) as fh:
        next(fh)
        for line in fh:
            ns, acc, chrom, start, end, nld, arch, fams, odd = \
                line.rstrip("\n").split("\t")
            f5 = ns in full5_present and ns in full5_set
            rows.append({
                "category": "intra_athila", "ns_id": ns, "accession": acc,
                "chrom": chrom, "start": start, "end": end,
                "architecture": arch, "event": f"{odd}->{fams}",
                "swapped": odd, "partners": fams, "confidence": "nt margin>=1.3",
                "full5": int(f5),
            })

    # ---- inter-clade (TEsorter per-domain clade) -------------------------- #
    for cf in sorted(glob.glob(os.path.join(args.data_dir, "*.cls.tsv"))):
        acc = L.accession_from_path(cf)
        for rec in L.parse_cls_tsv(cf):
            if rec["clade"] != "Athila":
                continue
            foreign = {d: c for d, c in rec["domains"].items() if c in LTR_CLADES}
            athila_core = [d for d in L.CORE_DOMAINS if rec["domains"].get(d) == "Athila"]
            if foreign and len(athila_core) >= 2:
                ns = L.nsid(acc, rec["element"])
                c = L.element_coords(rec["element"]) or {}
                arch = " ".join(f"{d}:{rec['domains'][d]}" for d in L.CORE_DOMAINS
                                if d in rec["domains"])
                fd = ";".join(f"{d}={c2}" for d, c2 in foreign.items())
                rows.append({
                    "category": "inter_clade", "ns_id": ns, "accession": acc,
                    "chrom": c.get("chrom", ""), "start": c.get("start", ""),
                    "end": c.get("end", ""), "architecture": arch,
                    "event": fd, "swapped": ";".join(foreign),
                    "partners": "Athila+" + ";".join(sorted(set(foreign.values()))),
                    "confidence": "TEsorter clade", "full5": int(ns in full5_present),
                })

    master = os.path.join(args.out_dir, "recombination_events.tsv")
    cols = ["category", "ns_id", "accession", "chrom", "start", "end",
            "architecture", "event", "swapped", "partners", "confidence", "full5"]
    with open(master, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["category"], r["chrom"],
                                             str(r["start"]))):
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- event-type collapse (recurrent copies -> one type) --------------- #
    def event_type(r):
        return (r["category"], r["partners"], r["swapped"])

    types = Counter(event_type(r) for r in rows)

    summary = os.path.join(args.out_dir, "recombination_summary.md")
    intra = [r for r in rows if r["category"] == "intra_athila"]
    inter = [r for r in rows if r["category"] == "inter_clade"]
    with open(summary, "w") as f:
        f.write("# Candidate recombination events — master list\n\n")
        f.write(f"- **intra-Athila (between ATHILA families):** {len(intra)} elements "
                f"({sum(r['full5'] for r in intra)} full-5)\n")
        f.write(f"- **inter-clade (Athila x other LTR retro family):** {len(inter)} "
                f"elements ({sum(r['full5'] for r in inter)} full-5)\n")
        f.write(f"- **total candidate elements:** {len(rows)}\n\n")
        f.write("## Distinct event types (recurrent copies collapsed)\n\n")
        f.write("| category | partners | swapped domain(s) | # elements |\n")
        f.write("|---|---|---|---:|\n")
        for (cat, partners, swapped), n in types.most_common():
            f.write(f"| {cat} | {partners} | {swapped} | {n} |\n")

    print(f"[list] intra-Athila: {len(intra)} ({sum(r['full5'] for r in intra)} full-5)")
    print(f"[list] inter-clade:  {len(inter)} ({sum(r['full5'] for r in inter)} full-5)")
    print(f"[list] distinct event types: {len(types)}")
    print(f"[list] wrote {master}, {summary}")


if __name__ == "__main__":
    main()
