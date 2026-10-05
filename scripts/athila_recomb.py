#!/usr/bin/env python3
"""Find recombination BETWEEN pipeline-defined Athila families (ATHILA0..9).

The Athilafinder pipeline defines Athila sub-families by consensus internal
sequences (athila_cns_internal: ATHILA0, 1, 2, 3, 4, 4c, 5, 6a/b, 7/7a, 8a/b, 9),
with per-family per-domain peptides. We assign EACH domain of EACH Athila element
to its best-matching ATHILA family (DIAMOND, best bitscore), then flag elements
whose domains map to DIFFERENT families = inter-family (within-Athila) recombinants.

No downsampling, no full-5 restriction: every Athila element with >=2 confidently
assigned domains is tested.

Usage:
    athila_recomb.py --data-dir DIR --consensus CONS.dom.faa --out-dir OUT \
        [--diamond ~/bin/diamond] [--min-bits 40] [--margin 1.15] [--threads 8]
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from collections import Counter, defaultdict

import lib_ids as L

Q_SEP = "@@"


def family_of(header_seqid):
    """'ATHILA4c.v1_I|...:Ty3-GAG' -> 'ATHILA4c.v1'."""
    fam = header_seqid.split("|", 1)[0]
    return fam[:-2] if fam.endswith("_I") else fam


def build_consensus_db(consensus_faa, workdir, diamond):
    """Write consensus domain peptides as FAMILY__DOMAIN and make a DIAMOND db."""
    out = os.path.join(workdir, "cons.faa")
    n = 0
    with open(out, "w") as o:
        for element, domain, attrs, seq in L.parse_dom_faa(consensus_faa):
            fam = family_of(element) if "|" not in element else family_of(element)
            # parse_dom_faa gives element = token before first '|'; that's FAMILY_I
            fam = element[:-2] if element.endswith("_I") else element
            dom = domain
            if dom:
                o.write(f">{fam}__{dom}\n{seq}\n")
                n += 1
    db = os.path.join(workdir, "cons")
    subprocess.run([diamond, "makedb", "--in", out, "--db", db, "--quiet"],
                   check=True)
    print(f"[recomb] consensus db: {n} family-domain peptides")
    return db


def build_query(data_dir, workdir):
    """All Athila elements' core-domain peptides (longest per element,domain)."""
    athila = set()
    for cf in sorted(glob.glob(os.path.join(data_dir, "*.cls.tsv"))):
        acc = L.accession_from_path(cf)
        for rec in L.parse_cls_tsv(cf):
            if L.is_athila(rec):
                athila.add(L.nsid(acc, rec["element"]))
    best = {}   # (ns, domain) -> (len, seq)
    for ff in sorted(glob.glob(os.path.join(data_dir, "*.dom.faa"))):
        acc = L.accession_from_path(ff)
        for element, domain, attrs, seq in L.parse_dom_faa(ff):
            if domain not in L.CORE_DOMAINS:
                continue
            ns = L.nsid(acc, element)
            if ns not in athila:
                continue
            k = (ns, domain)
            if k not in best or len(seq) > best[k][0]:
                best[k] = (len(seq), seq)
    q = os.path.join(workdir, "query.faa")
    with open(q, "w") as o:
        for (ns, domain), (_l, seq) in best.items():
            o.write(f">{ns}{Q_SEP}{domain}\n{seq}\n")
    print(f"[recomb] Athila elements: {len(athila)}; query domain peptides: {len(best)}")
    return q, athila


def assign_families(db, query, workdir, diamond, threads, min_bits, margin):
    """Best + second-best family per (element,domain) with a margin filter."""
    raw = os.path.join(workdir, "hits.tsv")
    subprocess.run([diamond, "blastp", "-q", query, "--db", db, "-o", raw,
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore",
                    "--very-sensitive", "--max-target-seqs", "25",
                    "--evalue", "1e-3", "--threads", str(threads), "--quiet"],
                   check=True)
    # best family bitscore per query (collapse family variants, e.g. ATHILA4c.v1)
    perq = defaultdict(dict)      # q -> {family: best_bits}
    with open(raw) as fh:
        for line in fh:
            q, s, b = line.rstrip("\n").split("\t")
            fam = s.split("__")[0]
            fam_base = fam.split(".")[0].rstrip("abcv")  # ATHILA4c.v1->ATHILA4, 6a->6? keep variant
            fam_use = fam.split(".")[0]                   # ATHILA4c.v1 -> ATHILA4c
            b = float(b)
            if b > perq[q].get(fam_use, 0):
                perq[q][fam_use] = b
    assign = {}
    for q, fams in perq.items():
        ranked = sorted(fams.items(), key=lambda kv: -kv[1])
        best_fam, best_b = ranked[0]
        second_b = ranked[1][1] if len(ranked) > 1 else 0.0
        confident = best_b >= min_bits and (second_b == 0 or best_b / max(second_b, 1e-9) >= margin)
        ns, domain = q.split(Q_SEP)
        assign[(ns, domain)] = (best_fam, best_b, second_b, confident)
    return assign


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--diamond", default=os.path.expanduser("~/bin/diamond"))
    ap.add_argument("--min-bits", type=float, default=40.0)
    ap.add_argument("--margin", type=float, default=1.15,
                    help="best/second-best bitscore ratio to call a confident family")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "athila_recomb")
    os.makedirs(out, exist_ok=True)
    tmp = os.path.join(out, "tmp")
    os.makedirs(tmp, exist_ok=True)

    db = build_consensus_db(args.consensus, tmp, args.diamond)
    query, athila = build_query(args.data_dir, tmp)
    assign = assign_families(db, query, tmp, args.diamond, args.threads,
                             args.min_bits, args.margin)

    # per element: domain -> family (confident only)
    by_elem = defaultdict(dict)
    for (ns, domain), (fam, b, b2, conf) in assign.items():
        if conf:
            by_elem[ns][domain] = fam

    # write per-domain family assignment
    with open(os.path.join(out, "domain_family_assignment.tsv"), "w") as f:
        f.write("ns_id\tdomain\tfamily\tbitscore\tsecond_bits\tconfident\n")
        for (ns, domain), (fam, b, b2, conf) in sorted(assign.items()):
            f.write(f"{ns}\t{domain}\t{fam}\t{b:.0f}\t{b2:.0f}\t{int(conf)}\n")

    # recombinants: >=2 confident domains, >=2 distinct families
    recs = []
    for ns, dmap in by_elem.items():
        fams = set(dmap.values())
        if len(dmap) >= 2 and len(fams) >= 2:
            recs.append((ns, dmap))
    print(f"[recomb] Athila elements with >=2 confident domains: "
          f"{sum(1 for d in by_elem.values() if len(d)>=2)}")
    print(f"[recomb] INTER-FAMILY recombinants (domains hit >=2 ATHILA families): "
          f"{len(recs)}")

    # which domain is the 'odd one out' + which family-pairs recombine
    odd_domain = Counter()
    pair_counter = Counter()
    with open(os.path.join(out, "recombinants.tsv"), "w") as f:
        f.write("ns_id\taccession\tchrom\tstart\tend\tn_domains\t"
                "architecture\tfamilies\todd_domain\n")
        for ns, dmap in sorted(recs):
            acc, element = L.split_nsid(ns)
            c = L.element_coords(element) or {}
            arch = ";".join(f"{d}:{dmap[d]}" for d in L.CORE_DOMAINS if d in dmap)
            famc = Counter(dmap.values())
            majfam, _ = famc.most_common(1)[0]
            odd = [d for d in L.CORE_DOMAINS if d in dmap and dmap[d] != majfam]
            for d in odd:
                odd_domain[d] += 1
            for fa in sorted(set(dmap.values())):
                for fb in sorted(set(dmap.values())):
                    if fa < fb:
                        pair_counter[(fa, fb)] += 1
            f.write(f"{ns}\t{acc}\t{c.get('chrom','')}\t{c.get('start','')}\t"
                    f"{c.get('end','')}\t{len(dmap)}\t{arch}\t"
                    f"{','.join(sorted(set(dmap.values())))}\t{','.join(odd) or 'NA'}\n")

    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("# Within-Athila inter-family recombination\n\n")
        f.write(f"- Athila elements tested (>=2 confident domains): "
                f"{sum(1 for d in by_elem.values() if len(d)>=2)}\n")
        f.write(f"- **inter-family recombinants: {len(recs)}**\n\n")
        f.write("## Odd-domain tally (domain that differs from the element majority)\n\n")
        for d in L.CORE_DOMAINS:
            f.write(f"- {d}: {odd_domain.get(d,0)}\n")
        f.write("\n## Top recombining family pairs\n\n")
        f.write("| family A | family B | elements |\n|---|---|---:|\n")
        for (fa, fb), n in pair_counter.most_common(15):
            f.write(f"| {fa} | {fb} | {n} |\n")

    print(f"[recomb] odd-domain tally: {dict(odd_domain)}")
    print(f"[recomb] wrote {out}/recombinants.tsv, summary.md")


if __name__ == "__main__":
    main()
