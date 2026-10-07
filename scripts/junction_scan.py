#!/usr/bin/env python3
"""Recombination junctions from split alignments of ATHILA copies to family consensus.

Every query sequence (full-length element, truncated element, or solo LTR) is aligned
with minimap2 (-x asm20, base-level) to the TAIR12 ATHILA exemplars (internal _I and
_LTR, labelled by family). A non-overlapping chain of alignments is picked greedily
by matching bases, ordered along the query, and every pair of consecutive pieces is
a junction:

  internal_deletion  same family+part, same strand, consensus jumps forward > --min-gap
  duplication        same family+part, same strand, consensus jumps backward > --min-gap
  inversion          same family, strand flips
  family_switch      pieces from different families (homologous-recombination-like
                     if long junction homology, illegitimate-like if microhomology)
  ltr_internal       expected LTR <-> internal boundary of the same family (not counted)

Junction homology = overlap of the two pieces on the query (microhomology when
small, long homology when large); a negative overlap is an untemplated insertion.

Usage:
    junction_scan.py --query Q.fasta[ Q2.fasta ...] --exemplars EX.fasta --out-dir OUT \
        [--minimap2 PATH] [--threads 16] [--tag full]
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from collections import Counter, defaultdict

import numpy as np


def ref_meta(name):
    """'ATHILA6a_I_nonauto#LTR/Ty3' -> ('ATHILA6a', 'I')."""
    base = name.split("#")[0]
    m = re.match(r"(ATHILA\d+[abc]?)", base)
    fam = m.group(1) if m else base
    part = "LTR" if "_LTR" in base else "I"
    return fam, part


def read_paf(path, min_match):
    by_q = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        q, qlen, qs, qe, strand, r, rlen, rs, re_, nmatch, alen = (
            f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]),
            int(f[7]), int(f[8]), int(f[9]), int(f[10]))
        if nmatch < min_match:
            continue
        fam, part = ref_meta(r)
        by_q[q].append({"q": q, "qlen": qlen, "qs": qs, "qe": qe, "strand": strand,
                        "ref": r, "fam": fam, "part": part, "rlen": rlen,
                        "rs": rs, "re": re_, "match": nmatch, "alen": alen,
                        "ident": nmatch / alen if alen else 0})
    return by_q


def chain(alns, max_overlap):
    chosen = []
    for a in sorted(alns, key=lambda x: -x["match"]):
        ok = True
        for c in chosen:
            ov = min(a["qe"], c["qe"]) - max(a["qs"], c["qs"])
            if ov > max_overlap and ov > 0.3 * min(a["qe"] - a["qs"], c["qe"] - c["qs"]):
                ok = False
                break
        if ok:
            chosen.append(a)
    return sorted(chosen, key=lambda x: x["qs"])


def classify(s1, s2, min_gap, clean_gap=20, min_piece=300):
    hom = s1["qe"] - s2["qs"]                     # >0 junction homology, <0 insertion
    if s1["fam"] != s2["fam"]:
        if s1["part"] != s2["part"]:
            return "ltr_internal_mismatch", hom, None
        if s1["part"] == "LTR":
            return "ltr_family_switch", hom, None
        clean = (-hom <= 50 and min(s1["qe"] - s1["qs"], s2["qe"] - s2["qs"]) >= min_piece)
        return ("internal_family_switch" if clean else "internal_family_switch_gapped"), hom, None
    if s1["strand"] != s2["strand"]:
        return "inversion", hom, None
    if s1["part"] != s2["part"]:
        return "ltr_internal", hom, None
    if s1["ref"] != s2["ref"]:
        return "same_family_other_exemplar", hom, None
    if s1["strand"] == "+":
        rgap = s2["rs"] - s1["re"]
    else:
        rgap = s1["rs"] - s2["re"]
    qgap = max(0, -hom)
    if rgap - qgap > min_gap:
        return ("deletion_clean" if qgap <= clean_gap else "deletion_with_insert"), hom, rgap - qgap
    if rgap < -min_gap:
        return "duplication", hom, -rgap
    return None, hom, None


def read_fa(path):
    d, n, b = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if n:
                d[n] = "".join(b).upper()
            n, b = line[1:].split()[0], []
        else:
            b.append(line.strip())
    if n:
        d[n] = "".join(b).upper()
    return d


def microhomology(seq, a, b, cap=30):
    """Identical bases shared by the two breakpoints a<b of a deletion in seq."""
    k = 0
    while k < cap and a - k - 1 >= 0 and seq[a - k - 1] == seq[b - k - 1]:
        k += 1
    r = 0
    while k + r < cap and b + r < len(seq) and seq[a + r] == seq[b + r]:
        r += 1
    return k + r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", nargs="+", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tag", default="full")
    ap.add_argument("--minimap2", default=os.path.expanduser("~/minimap2/minimap2"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-match", type=int, default=80)
    ap.add_argument("--max-overlap", type=int, default=300)
    ap.add_argument("--min-gap", type=int, default=100)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, f"junctions_{args.tag}")
    os.makedirs(out, exist_ok=True)
    paf = os.path.join(out, "aln.paf")
    if not os.path.exists(paf):
        with open(paf, "w") as fh:
            subprocess.run([args.minimap2, "-x", "asm20", "-c", "--secondary=no",
                            "-t", str(args.threads), args.exemplars, *args.query],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
    by_q = read_paf(paf, args.min_match)
    refseq = read_fa(args.exemplars)
    import random
    rng = random.Random(0)
    mh_null = []

    junctions, elements = [], []
    for q, alns in by_q.items():
        ch = chain(alns, args.max_overlap)
        qlen = ch[0]["qlen"]
        # family/part composition by aligned query bases
        cov = Counter()
        for a in ch:
            cov[(a["fam"], a["part"])] += a["qe"] - a["qs"]
        best_I = max(((f, n) for (f, p), n in cov.items() if p == "I"),
                     key=lambda x: x[1], default=(None, 0))
        best_L = max(((f, n) for (f, p), n in cov.items() if p == "LTR"),
                     key=lambda x: x[1], default=(None, 0))
        # internal consensus coverage (union of ref intervals on the main internal ref)
        icov = 0.0
        if best_I[0]:
            main = Counter(a["ref"] for a in ch if a["fam"] == best_I[0] and a["part"] == "I")
            mref = main.most_common(1)[0][0]
            iv = sorted((a["rs"], a["re"]) for a in ch if a["ref"] == mref)
            tot, cs, ce = 0, None, None
            for s, e in iv:
                if cs is None or s > ce:
                    if cs is not None:
                        tot += ce - cs
                    cs, ce = s, e
                else:
                    ce = max(ce, e)
            tot += ce - cs
            icov = tot / next(a["rlen"] for a in ch if a["ref"] == mref)
        n_j = Counter()
        for s1, s2 in zip(ch, ch[1:]):
            cls, hom, size = classify(s1, s2, args.min_gap)
            if cls is None:
                continue
            n_j[cls] += 1
            mh = ""
            if cls == "deletion_clean":
                seq = refseq.get(s1["ref"], "")
                a, b = ((s1["re"], s2["rs"]) if s1["strand"] == "+" else (s2["re"], s1["rs"]))
                if seq and 0 < a < b < len(seq):
                    mh = microhomology(seq, a, b)
                    for _ in range(20):                      # null: same distance, random place
                        a0 = rng.randrange(1, len(seq) - (b - a) - 1)
                        mh_null.append(microhomology(seq, a0, a0 + (b - a)))
            rel = None
            if s1["part"] == "I":
                rel = (s1["re"] if s1["strand"] == "+" else s1["rs"]) / s1["rlen"]
            junctions.append({"query": q, "class": cls, "homology": hom,
                              "size": size if size is not None else "",
                              "fam1": s1["fam"], "part1": s1["part"],
                              "fam2": s2["fam"], "part2": s2["part"],
                              "qpos": s1["qe"], "rel_pos_in_consensus1":
                              f"{rel:.3f}" if rel is not None else "",
                              "ident1": f"{s1['ident']:.3f}",
                              "ident2": f"{s2['ident']:.3f}",
                              "piece1_bp": s1["qe"] - s1["qs"], "piece2_bp": s2["qe"] - s2["qs"],
                              "microhomology_exact": mh})
        elements.append({"query": q, "qlen": qlen, "n_pieces": len(ch),
                         "internal_family": best_I[0] or "", "ltr_family": best_L[0] or "",
                         "internal_bp": best_I[1], "ltr_bp": best_L[1],
                         "internal_consensus_cov": f"{icov:.3f}",
                         **{f"n_{k}": n_j.get(k, 0) for k in
                            ("deletion_clean", "deletion_with_insert",
                             "internal_family_switch", "internal_family_switch_gapped",
                             "ltr_internal_mismatch", "ltr_family_switch", "inversion",
                             "duplication")}})

    for name, rows in (("junctions.tsv", junctions), ("elements.tsv", elements)):
        if not rows:
            continue
        cols = list(rows[0])
        with open(os.path.join(out, name), "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")

    jc = Counter(j["class"] for j in junctions)
    print(f"[junction:{args.tag}] queries aligned: {len(elements)}; junctions: {dict(jc)}")
    if mh_null:
        with open(os.path.join(out, "microhomology_null.txt"), "w") as f:
            f.write("\n".join(map(str, mh_null)) + "\n")
        obs = np.array([j["microhomology_exact"] for j in junctions
                        if j["class"] == "deletion_clean" and j["microhomology_exact"] != ""])
        nul = np.array(mh_null)
        print(f"  exact microhomology at clean deletions: median {np.median(obs):.0f} bp "
              f"(null {np.median(nul):.0f}); >=2bp {np.mean(obs>=2)*100:.0f}% vs null "
              f"{np.mean(nul>=2)*100:.0f}%; >=4bp {np.mean(obs>=4)*100:.0f}% vs "
              f"{np.mean(nul>=4)*100:.0f}%")
    for cls in ("deletion_clean", "deletion_with_insert", "internal_family_switch",
                "internal_family_switch_gapped", "ltr_internal_mismatch", "inversion",
                "duplication"):
        hs = np.array([j["homology"] for j in junctions if j["class"] == cls])
        if len(hs):
            print(f"  {cls:>18}: n={len(hs)}  junction homology median {np.median(hs):.0f} bp; "
                  f"<=10 bp (micro/none) {np.mean(hs <= 10)*100:.0f}%, "
                  f">=50 bp (long) {np.mean(hs >= 50)*100:.0f}%")


if __name__ == "__main__":
    main()
