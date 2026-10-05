#!/usr/bin/env python3
"""Phase 1 — the free screen: per-domain clade disagreement.

TEsorter assigns a clade to every domain hit (the ``Domains`` column of cls.tsv).
Before building any network we ask: within one element, do the domains agree on a
clade? This is a cheap baseline expectation for Phase 4 — NOT the answer, because
REXdb sub-clade resolution within Athila is coarse.

Usage:
    p1_clade_disagree.py --input-dir DIR --out-dir OUT
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
    args = ap.parse_args()

    inv = os.path.join(args.out_dir, "inventory")
    os.makedirs(inv, exist_ok=True)

    domain_clade = defaultdict(Counter)   # domain -> Counter(clade)
    pair_disagree = Counter()             # (domA,domB) -> # elements disagreeing
    pair_total = Counter()                # (domA,domB) -> # elements co-occurring
    n_elements = n_disagree = 0
    disagree_rows = []

    for cf in sorted(glob.glob(os.path.join(args.input_dir, "*.cls.tsv"))):
        acc = L.accession_from_path(cf)
        for rec in L.parse_cls_tsv(cf):
            dom2clade = rec["domains"]
            if not dom2clade:
                continue
            n_elements += 1
            for d, c in dom2clade.items():
                domain_clade[d][c] += 1
            doms = sorted(dom2clade, key=lambda d: (L.CORE_DOMAINS.index(d)
                          if d in L.CORE_DOMAINS else 99))
            # pairwise co-occurrence / disagreement
            for i in range(len(doms)):
                for j in range(i + 1, len(doms)):
                    a, b = doms[i], doms[j]
                    pair_total[(a, b)] += 1
                    if dom2clade[a] != dom2clade[b]:
                        pair_disagree[(a, b)] += 1
            distinct = set(dom2clade.values())
            if len(distinct) > 1:
                n_disagree += 1
                disagree_rows.append((L.nsid(acc, rec["element"]),
                                      ";".join(f"{d}:{dom2clade[d]}"
                                               for d in doms)))

    # write disagreeing elements
    dpath = os.path.join(inv, "disagreeing_elements.tsv")
    with open(dpath, "w") as out:
        out.write("ns_id\tdomain_clades\n")
        for ns, s in disagree_rows:
            out.write(f"{ns}\t{s}\n")

    # domain x clade matrix
    all_clades = sorted({c for cc in domain_clade.values() for c in cc})
    mpath = os.path.join(inv, "clade_by_domain.tsv")
    with open(mpath, "w") as out:
        out.write("domain\t" + "\t".join(all_clades) + "\n")
        for d in L.CORE_DOMAINS + [d for d in domain_clade
                                   if d not in L.CORE_DOMAINS]:
            if d not in domain_clade:
                continue
            out.write(d + "\t" + "\t".join(str(domain_clade[d].get(c, 0))
                      for c in all_clades) + "\n")

    md = os.path.join(inv, "clade_disagreement.md")
    with open(md, "w") as out:
        out.write("# Phase 1 — per-domain clade disagreement\n\n")
        out.write(f"- elements with >=1 domain clade call: **{n_elements}**\n")
        pct = 100.0 * n_disagree / n_elements if n_elements else 0
        out.write(f"- elements whose domains disagree on clade: "
                  f"**{n_disagree}** ({pct:.1f}%)\n\n")
        out.write("> Baseline only. A much larger Phase-4 discordance would be a "
                  "sequence-level signal, not a clade-call artefact.\n\n")
        out.write("## Top disagreeing domain pairs\n\n")
        out.write("| domain pair | disagree / co-occur | rate |\n|---|---:|---:|\n")
        for pair, tot in pair_total.most_common():
            dis = pair_disagree.get(pair, 0)
            if dis:
                out.write(f"| {pair[0]}-{pair[1]} | {dis} / {tot} | "
                          f"{100.0*dis/tot:.1f}% |\n")
        out.write(f"\nFull matrix: `clade_by_domain.tsv`; "
                  f"disagreeing elements: `disagreeing_elements.tsv`\n")

    print(f"[p1] {n_disagree}/{n_elements} elements with domain clade disagreement")
    print(f"[p1] wrote {md}")


if __name__ == "__main__":
    main()
