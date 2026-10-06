#!/usr/bin/env python3
"""k-mer marker mosaic scan for Athila inter-family recombination.

Idea (per J.G.): build family-specific k-mer marker sets (a k-mer is a marker for
family F if it occurs in F's consensus/exemplars and in NO other family -- the
set-difference / KMC-intersection logic), then slide along each element and record
which family each marker k-mer belongs to. A spatially segregated switch from one
family's markers to another's = a recombination mosaic. This resolves even >90%-
similar families (every ~31-mer spans a diagnostic difference) and needs no
alignment and no bitscore margin, so it has no high-homology blind spot.

Usage:
    kmer_recomb.py --exemplars EX.fasta --fulllength DIR --out-dir OUT \
        [--k 31] [--nbins 24] [--min-seg 8] [--min-total 20] [--dom-frac 0.6]
"""

from __future__ import annotations

import argparse
import glob
import os
from collections import Counter, defaultdict

import lib_ids as L
from athila_recomb_nt import read_fasta, family_of, revcomp


def build_markers(exemplars, k):
    """family -> set(kmers); then markers unique to one family (both strands)."""
    fam_k = defaultdict(set)
    for h, s in read_fasta(exemplars).items():
        if "_I" in h and "_LTR" not in h:
            fam = family_of(h)
            s = s.upper()
            for i in range(len(s) - k + 1):
                km = s[i:i + k]
                if set(km) <= set("ACGT"):
                    fam_k[fam].add(km)
    fams = sorted(fam_k)
    marker = {}            # kmer -> family (specific markers only), both strands
    for f in fams:
        others = set().union(*[fam_k[g] for g in fams if g != f])
        for km in fam_k[f] - others:
            marker[km] = f
            marker[revcomp(km)] = f
    return marker, fams


def scan_element(seq, marker, k, nbins, min_seg, min_total, dom_frac):
    """Return (primary_family, mosaic_segments) for one element sequence."""
    seq = seq.upper()
    hits = []                      # (pos, family)
    for i in range(len(seq) - k + 1):
        f = marker.get(seq[i:i + k])
        if f:
            hits.append((i, f))
    if len(hits) < min_total:
        prim = Counter(f for _, f in hits).most_common(1)
        return (prim[0][0] if prim else None), None

    counts = Counter(f for _, f in hits)
    primary = counts.most_common(1)[0][0]

    # bin along the element; per-bin dominant family
    L_ = len(seq)
    binsz = max(1, L_ // nbins)
    binfam = {}
    bins = defaultdict(Counter)
    for pos, f in hits:
        bins[pos // binsz][f] += 1
    for b, c in bins.items():
        tot = sum(c.values())
        fam, n = c.most_common(1)[0]
        if n >= 0.5 * tot:         # a clear local majority
            binfam[b] = fam

    # collapse consecutive bins into segments
    segs = []
    for b in sorted(binfam):
        f = binfam[b]
        if segs and segs[-1][0] == f and b == segs[-1][2] + 1:
            segs[-1][1] += bins[b][f]; segs[-1][2] = b
        else:
            segs.append([f, bins[b][f], b, b])   # [family, markers, first_bin, last_bin]

    strong = [s for s in segs if s[1] >= min_seg]
    distinct = []
    for s in strong:
        if not distinct or distinct[-1][0] != s[0]:
            distinct.append(s)
    if len({s[0] for s in distinct}) >= 2:
        return primary, distinct
    return primary, None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--markers", default=None,
                    help="robust marker TSV (kmer<TAB>family) from kmc_markers.py")
    ap.add_argument("--tag", default="kmer_recomb", help="output subdirectory name")
    ap.add_argument("--clip-dir", default=None,
                    help="TEsorter dir with *.dom.gff3: scan only the internal domain span")
    ap.add_argument("--fulllength", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--min-seg", type=int, default=8)
    ap.add_argument("--min-total", type=int, default=20)
    ap.add_argument("--dom-frac", type=float, default=0.6)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, args.tag)
    os.makedirs(out, exist_ok=True)

    if args.markers:
        marker = {}
        with open(args.markers) as fh:
            next(fh)
            for line in fh:
                km, fam = line.rstrip("\n").split("\t")
                km = km.upper()
                marker[km] = fam
                marker[revcomp(km)] = fam
        fams = sorted(set(marker.values()))
    else:
        marker, fams = build_markers(args.exemplars, args.k)
    print(f"[kmer] families: {len(fams)}; marker k-mers (both strands): {len(marker)}")

    span = {}
    if args.clip_dir:
        for gff in glob.glob(os.path.join(args.clip_dir, "*.dom.gff3")):
            acc = L.accession_from_path(gff)
            for line in open(gff):
                f = line.split("\t")
                if len(f) < 9:
                    continue
                ns = L.nsid(acc, f[0])
                a, b = int(f[3]), int(f[4])
                lo, hi = span.get(ns, (a, b))
                span[ns] = (min(lo, a), max(hi, b))
        print(f"[kmer] clipping to internal domain span for {len(span)} elements")

    primary = {}            # ns_id -> family (abundance)
    recs = []               # (ns_id, segments)
    n_scanned = 0
    for flf in sorted(glob.glob(os.path.join(args.fulllength,
                                             "*DP_FULLLENGTH_RENAMED.fasta"))):
        acc = L.accession_from_path(flf)
        for el, seq in read_fasta(flf).items():
            ns = L.nsid(acc, el)
            if args.clip_dir:
                if ns not in span:
                    continue
                lo, hi = span[ns]
                seq = seq[lo - 1:hi]
            prim, segs = scan_element(seq, marker, args.k, args.nbins,
                                      args.min_seg, args.min_total, args.dom_frac)
            n_scanned += 1
            if prim:
                primary[ns] = prim
            if segs:
                recs.append((ns, segs))

    # write primary family (abundance) + recombinants
    with open(os.path.join(out, "primary_family.tsv"), "w") as f:
        f.write("ns_id\tprimary_family\n")
        for ns, fam in primary.items():
            f.write(f"{ns}\t{fam}\n")

    pair_cnt = Counter()
    odd = Counter()
    with open(os.path.join(out, "kmer_recombinants.tsv"), "w") as f:
        f.write("ns_id\taccession\tchrom\tstart\tend\tn_segments\t"
                "mosaic\tfamilies\tbreakpoint_bins\n")
        for ns, segs in sorted(recs):
            acc, el = L.split_nsid(ns)
            c = L.element_coords(el) or {}
            mosaic = ">".join(f"{s[0][6:]}({s[1]})" for s in segs)
            famset = sorted({s[0] for s in segs})
            bps = ";".join(str(segs[i][3]) for i in range(len(segs) - 1))
            for a in range(len(famset)):
                for b in range(a + 1, len(famset)):
                    pair_cnt[(famset[a], famset[b])] += 1
            f.write(f"{ns}\t{acc}\t{c.get('chrom','')}\t{c.get('start','')}\t"
                    f"{c.get('end','')}\t{len(segs)}\t{mosaic}\t"
                    f"{','.join(x[6:] for x in famset)}\t{bps}\n")

    ab = Counter(primary.values())
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("# k-mer marker mosaic scan — Athila inter-family recombination\n\n")
        f.write(f"- elements scanned: {n_scanned}\n")
        f.write(f"- elements with a primary family (>= {args.min_total} markers): "
                f"{len(primary)}\n")
        f.write(f"- **mosaic (recombinant) elements: {len(recs)}**\n\n")
        f.write("## Primary-family abundance\n\n| family | n |\n|---|---:|\n")
        for fam, n in ab.most_common():
            f.write(f"| {fam[6:]} | {n} |\n")
        f.write("\n## Top recombining family pairs (mosaic)\n\n| A | B | n |\n|---|---|---:|\n")
        for (a, b), n in pair_cnt.most_common(20):
            f.write(f"| {a[6:]} | {b[6:]} | {n} |\n")

    print(f"[kmer] scanned {n_scanned}; primary-typed {len(primary)}; "
          f"mosaics {len(recs)}")
    print(f"[kmer] top pairs: " +
          ", ".join(f"{a[6:]}-{b[6:]}:{n}" for (a, b), n in pair_cnt.most_common(6)))
    print(f"[kmer] wrote {out}/")


if __name__ == "__main__":
    main()
