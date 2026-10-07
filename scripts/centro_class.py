#!/usr/bin/env python3
"""Centrophilic / centrophobic classification of ATHILA elements and families.

Element level (per accession, per chromosome):
  core        = the CEN178 array cluster with most repeats (arrays merged if <= --merge-gap)
  centromeric = element overlaps the core span
  pericentromeric = outside the core, within --peri-bp of it
  arm         = further away

Family level: in each accession, a family's centromeric fraction is compared with
that accession's overall ATHILA centromeric fraction (log2 ratio). Families whose
ratio is consistently > 0 across accessions (sign test, BH-corrected) are
CENTROPHILIC; consistently < 0 are CENTROPHOBIC; otherwise NEUTRAL.

Also tests whether recombinant mosaics are enriched in centromeres.

Usage:
    centro_class.py --arrays cen178_arrays.tsv --elements elements.tsv \
        --primary primary_family.tsv --mosaics kmer_recombinants.tsv --out-dir OUT
"""

from __future__ import annotations

import argparse
import math
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import binomtest, fisher_exact


def short(f):
    return f[6:] if f.startswith("ATHILA") else f


def cores(path, merge_gap):
    arr = defaultdict(list)
    for l in open(path).readlines()[1:]:
        a, c, s, e, n, _ = l.rstrip("\n").split("\t")
        arr[(a, c)].append((int(s), int(e), int(n)))
    out = {}
    for k, xs in arr.items():
        xs.sort()
        cl, cur = [], [xs[0][0], xs[0][1], xs[0][2]]
        for s, e, n in xs[1:]:
            if s - cur[1] <= merge_gap:
                cur[1] = max(cur[1], e); cur[2] += n
            else:
                cl.append(cur); cur = [s, e, n]
        cl.append(cur)
        out[k] = max(cl, key=lambda x: x[2])
    return out


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for r, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--mosaics", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--merge-gap", type=int, default=500_000)
    ap.add_argument("--peri-bp", type=int, default=2_000_000)
    ap.add_argument("--min-per-acc", type=int, default=3)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    core = cores(args.arrays, args.merge_gap)
    fam = {l.split("\t")[0]: short(l.rstrip("\n").split("\t")[1])
           for l in open(args.primary).readlines()[1:]}
    mosaic = {l.split("\t")[0] for l in open(args.mosaics).readlines()[1:]}

    elems = {}
    for l in open(args.elements).readlines()[1:]:
        f = l.rstrip("\n").split("\t")
        ns, acc, chrom, s, e = f[0], f[1], f[3], f[4], f[5]
        if (acc, chrom) not in core or not s or ns not in fam:
            continue
        s, e = int(s), int(e)
        cs, ce, _ = core[(acc, chrom)]
        dist = 0 if (e >= cs and s <= ce) else min(abs(s - ce), abs(cs - e))
        cls = ("centromeric" if dist == 0 else
               "pericentromeric" if dist <= args.peri_bp else "arm")
        elems[ns] = {"acc": acc, "chrom": chrom, "fam": fam[ns], "dist": dist,
                     "class": cls, "mosaic": ns in mosaic}

    with open(os.path.join(args.out_dir, "element_centro_class.tsv"), "w") as f:
        f.write("ns_id\taccession\tchrom\tfamily\tdist_to_core_bp\tclass\tmosaic\n")
        for ns, r in elems.items():
            f.write(f"{ns}\t{r['acc']}\t{r['chrom']}\t{r['fam']}\t{r['dist']}\t"
                    f"{r['class']}\t{int(r['mosaic'])}\n")

    # family classification across accessions
    by_acc = defaultdict(list)
    for r in elems.values():
        by_acc[r["acc"]].append(r)
    fams = sorted({r["fam"] for r in elems.values()})
    rows = []
    for fm in fams:
        lr, pos, neg = [], 0, 0
        n_all = sum(1 for r in elems.values() if r["fam"] == fm)
        c_all = sum(1 for r in elems.values() if r["fam"] == fm and r["class"] == "centromeric")
        a_all = sum(1 for r in elems.values() if r["fam"] == fm and r["class"] == "arm")
        for acc, rs in by_acc.items():
            mine = [r for r in rs if r["fam"] == fm]
            if len(mine) < args.min_per_acc:
                continue
            fc = (sum(r["class"] == "centromeric" for r in mine) + 0.5) / (len(mine) + 1)
            bg = (sum(r["class"] == "centromeric" for r in rs) + 0.5) / (len(rs) + 1)
            x = math.log2(fc / bg)
            lr.append(x); pos += x > 0; neg += x < 0
        n_acc = pos + neg
        p_phil = binomtest(pos, n_acc, 0.5, alternative="greater").pvalue if n_acc else 1
        p_phob = binomtest(neg, n_acc, 0.5, alternative="greater").pvalue if n_acc else 1
        rows.append({"family": fm, "n_elements": n_all,
                     "frac_centromeric": c_all / n_all, "frac_arm": a_all / n_all,
                     "n_acc_tested": len(lr),
                     "median_log2_enrich": float(np.median(lr)) if lr else float("nan"),
                     "acc_enriched": pos, "acc_depleted": neg,
                     "p_phil": p_phil, "p_phob": p_phob})
    qp = bh([r["p_phil"] for r in rows]); qn = bh([r["p_phob"] for r in rows])
    for r, a, b in zip(rows, qp, qn):
        r["q_phil"], r["q_phob"] = a, b
        r["class"] = ("centrophilic" if a < 0.05 else
                      "centrophobic" if b < 0.05 else
                      "neutral" if r["n_acc_tested"] else "untested")

    cols = ["family", "class", "n_elements", "frac_centromeric", "frac_arm",
            "median_log2_enrich", "n_acc_tested", "acc_enriched", "acc_depleted",
            "q_phil", "q_phob"]
    with open(os.path.join(args.out_dir, "family_centro_class.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: -r["frac_centromeric"]):
            f.write("\t".join(f"{r[c]:.3g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")

    # are recombinants enriched in centromeres?
    tot = Counter(r["class"] for r in elems.values())
    mos = Counter(r["class"] for r in elems.values() if r["mosaic"])
    lines = ["# Centrophilic / centrophobic classification\n",
             f"- elements placed: {len(elems)} in {len(by_acc)} accessions",
             f"- classes: " + ", ".join(f"{k} {v} ({100*v/len(elems):.1f}%)" for k, v in tot.most_common()),
             "\n## Families\n",
             "| family | class | n | % centromeric | % arm | median log2 enrich | accessions +/- | q phil | q phob |",
             "|---|---|---:|---:|---:|---:|---|---:|---:|"]
    for r in sorted(rows, key=lambda r: -r["frac_centromeric"]):
        lines.append(f"| {r['family']} | **{r['class']}** | {r['n_elements']} | "
                     f"{100*r['frac_centromeric']:.0f} | {100*r['frac_arm']:.0f} | "
                     f"{r['median_log2_enrich']:+.2f} | {r['acc_enriched']}/{r['acc_depleted']} | "
                     f"{r['q_phil']:.2g} | {r['q_phob']:.2g} |")
    lines.append("\n## Are recombinant mosaics enriched in centromeres?\n")
    lines.append("| class | all elements | mosaics | % mosaic |\n|---|---:|---:|---:|")
    for k in ("centromeric", "pericentromeric", "arm"):
        lines.append(f"| {k} | {tot[k]} | {mos[k]} | {100*mos[k]/max(1,tot[k]):.2f} |")
    a, b = mos["centromeric"], tot["centromeric"] - mos["centromeric"]
    c = sum(mos.values()) - a
    d = (len(elems) - tot["centromeric"]) - c
    orr, p = fisher_exact([[a, b], [c, d]])
    lines.append(f"\nCentromeric vs not: odds ratio {orr:.2f}, Fisher p = {p:.2g}")
    open(os.path.join(args.out_dir, "summary.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
