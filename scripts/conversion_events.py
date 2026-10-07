#!/usr/bin/env python3
"""Per-event look at conversion-like (A>B>A) calls.

For each candidate:
  * marker islands along the element (independent sites, family-labelled);
  * per-event shuffle p-value: shuffle the element's own island labels N times and
    count how often a donor tract with >= the observed number of sites, flanked by
    host runs, appears (same classifier, same min-sites);
  * an orthogonal check: % identity of the element to the host and donor TAIR12
    internal exemplars in sliding windows (minimap2 --cs); a real tract should match
    the donor better than the host inside the tract and not outside;
  * context: centromeric class, LTR identity (age), copies of the same event in
    other accessions.

Usage: conversion_events.py --calls calls.tsv --min-sites 2 --markers M --fulllength DIR
           --clip-dir DIR --exemplars EX --results R --out-dir OUT
"""
import argparse
import glob
import os
import random
import re
import subprocess
import tempfile
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import lib_ids as L
from athila_recomb_nt import read_fasta
from island_scan_lib import islands, path_sites
from tract_scan import load_markers, load_spans, classify

ap = argparse.ArgumentParser()
for a in ("--calls", "--markers", "--fulllength", "--clip-dir", "--exemplars", "--results", "--out-dir"):
    ap.add_argument(a, required=True)
ap.add_argument("--min-sites", type=int, default=2)
ap.add_argument("--perms", type=int, default=10000)
ap.add_argument("--k", type=int, default=31)
ap.add_argument("--window", type=int, default=150)
ap.add_argument("--minimap2", default=os.path.expanduser("~/minimap2/minimap2"))
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
rng = random.Random(0)

calls = []
for l in open(args.calls).readlines()[1:]:
    f = l.rstrip("\n").split("\t")
    if int(f[0]) == args.min_sites and f[2] == "conversion":
        calls.append(f[1])
marker = load_markers(args.markers)
span = load_spans(args.clip_dir)
ex = read_fasta(args.exemplars)
ex_I = defaultdict(list)
for h, sq in ex.items():
    if "_I" in h and "_LTR" not in h:
        fam = re.match(r"ATHILA(\d+[abc]?)", h).group(1)
        ex_I[fam].append((h, sq.upper()))

age = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{args.results}/age/element_age.tsv").readlines()[1:]}
ctx = {l.split("\t")[0]: l.split("\t")[5] for l in open(f"{args.results}/centro/element_centro_class.tsv").readlines()[1:]}

seqs = {}
for flf in glob.glob(os.path.join(args.fulllength, "*DP_FULLLENGTH_RENAMED.fasta")):
    acc = L.accession_from_path(flf)
    want = [c for c in calls if c.startswith(acc + "::")]
    if want:
        fa = read_fasta(flf)
        for c in want:
            lo, hi = span[c]
            seqs[c] = fa[L.split_nsid(c)[1]][lo - 1:hi].upper()


def identity_track(seq, refs, win):
    """Per-window identity of seq to the best of several refs (minimap2 --cs; unaligned = NaN)."""
    with tempfile.TemporaryDirectory() as t:
        open(f"{t}/q.fa", "w").write(f">q\n{seq}\n")
        with open(f"{t}/r.fa", "w") as o:
            for i, (h, r) in enumerate(refs):
                o.write(f">r{i}\n{r}\n")
        out = subprocess.run([args.minimap2, "-x", "asm20", "-k", "11", "-w", "5", "-s", "40",
                              "-c", "--cs", "-N", "50", "-p", "0.1", f"{t}/r.fa", f"{t}/q.fa"],
                             capture_output=True, text=True).stdout
    best = None
    tracks = []
    for line in out.splitlines():
        state = np.full(len(seq), np.nan)
        f = line.split("\t")
        qs, qe, strand = int(f[2]), int(f[3]), f[4]
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        ops = re.findall(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+)", cs)
        q = qs if strand == "+" else qe - 1
        step = 1 if strand == "+" else -1
        for op in ops:
            if op[0] == ":":
                for _ in range(int(op[1:])):
                    if np.isnan(state[q]): state[q] = 1
                    q += step
            elif op[0] == "*":
                if np.isnan(state[q]): state[q] = 0
                q += step
            elif op[0] == "+":
                for _ in range(len(op) - 1):
                    if np.isnan(state[q]): state[q] = 0
                    q += step
        tracks.append(state)
    x = np.arange(0, len(seq) - win, win // 3)
    ys = []
    for state in tracks:
        ys.append([np.nanmean(state[i:i + win]) if np.isfinite(state[i:i + win]).sum() > win * 0.5 else np.nan
                   for i in x])
    y = np.nanmax(np.array(ys), axis=0) if ys else np.full(len(x), np.nan)
    return x + win / 2, y * 100


rows = []
for c in calls:
    s = seqs[c]
    hits = [(i, marker[s[i:i + args.k]]) for i in range(len(s) - args.k + 1) if s[i:i + args.k] in marker]
    isl = islands(hits)
    path = path_sites(isl, args.min_sites)
    host, donor = path[0][0], path[1][0]
    tract_sites = path[1][1]
    t_start, t_end = path[1][2], path[1][3] + args.k
    # per-event shuffle p-value
    labs = [x[0] for x in isl]; posn = [(a, b) for _, a, b in isl]
    hit = 0
    for _ in range(args.perms):
        sh = labs[:]; rng.shuffle(sh)
        p = path_sites([(f, a, b) for f, (a, b) in zip(sh, posn)], args.min_sites)
        if classify(p) == "conversion" and p[1][1] >= tract_sites:
            hit += 1
    pval = (hit + 1) / (args.perms + 1)
    n_host = sum(1 for f in labs if f == host); n_donor = sum(1 for f in labs if f == donor)
    # identity tracks
    xh, yh = identity_track(s, ex_I[host], args.window)
    xd, yd = identity_track(s, ex_I[donor], args.window)
    ins = (xh >= t_start) & (xh <= t_end)
    d_in = np.nanmean(yd[ins] - yh[ins]) if ins.any() else np.nan
    d_out = np.nanmean(yd[~ins] - yh[~ins])
    acc, el = L.split_nsid(c)
    co = L.element_coords(el) or {}
    rows.append({"element": c, "accession": acc, "chrom": co.get("chrom"), "start": co.get("start"),
                 "host": host, "donor": donor, "sites_host": n_host, "sites_donor": n_donor,
                 "tract_sites": tract_sites, "tract_bp": t_end - t_start,
                 "tract_start_in_span": t_start, "span_bp": len(s), "shuffle_p": pval,
                 "donor_minus_host_identity_in_tract": d_in,
                 "donor_minus_host_identity_outside": d_out,
                 "ltr_identity": age.get(c, "NA"), "context": ctx.get(c, "NA")})
    # figure
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 4.2), sharex=True, gridspec_kw={"height_ratios": [1, 2.3]},
                                 facecolor="#fcfcfb")
    col = {host: "#9a9890", donor: "#e34948"}
    for f, a, b in isl:
        a1.plot([a, b + args.k], [1 if f == donor else 0] * 2, color=col.get(f, "#2a78d6"), lw=6, solid_capstyle="butt")
    a1.set_yticks([0, 1], [f"host ATHILA{host}", f"donor ATHILA{donor}"]); a1.set_ylim(-0.6, 1.6)
    a1.axvspan(t_start, t_end, color="#e34948", alpha=0.12)
    a1.set_title(f"{c}   tract {tract_sites} sites / {t_end - t_start} bp   shuffle p = {pval:.2g}", loc="left", fontsize=9.5)
    a2.plot(xh, yh, color="#52514e", lw=1.6, label=f"identity to ATHILA{host} (host)")
    a2.plot(xd, yd, color="#e34948", lw=1.6, label=f"identity to ATHILA{donor} (donor)")
    a2.axvspan(t_start, t_end, color="#e34948", alpha=0.12)
    a2.set_ylabel("% identity (150 bp windows)"); a2.set_xlabel("position in element internal span (bp)")
    a2.legend(frameon=False, fontsize=8, loc="lower left")
    for ax in (a1, a2):
        ax.set_facecolor("#fcfcfb")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    a2.grid(axis="y", color="#e6e5e0", lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, re.sub(r"[^A-Za-z0-9_.-]", "_", c) + ".png"), dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

cols = list(rows[0])
with open(os.path.join(args.out_dir, "conversion_events.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(f"{r[k]:.3g}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
for r in rows:
    print(f"{r['element'][:62]:62s} {r['host']:>3}>{r['donor']:<3} tract {r['tract_sites']:>2} sites {r['tract_bp']:>5} bp | "
          f"p={r['shuffle_p']:.3g} | donor-host id in tract {r['donor_minus_host_identity_in_tract']:+.1f} "
          f"outside {r['donor_minus_host_identity_outside']:+.1f} | {r['context']}, LTR id {r['ltr_identity']}")
