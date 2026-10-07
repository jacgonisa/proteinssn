#!/usr/bin/env python3
"""Internal region (between the LTRs) of every Athilafinder intact element.

Full-length sequences are in element orientation (minus-strand elements reverse-
complemented), so internal = [LTR5_length, length - LTR3_length).
Writes one FASTA per genome with the same names as the full-length files.

Usage: internal_regions.py SUMMARY_DIR FULLLENGTH_DIR OUT_DIR
"""
import glob, os, sys
import lib_ids as L
from athila_recomb_nt import read_fasta
summ, fl_dir, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
n = skipped = 0
for p in sorted(glob.glob(os.path.join(summ, "*SUMMARY_TABLE.txt"))):
    acc = L.accession_from_path(p)
    flf = glob.glob(os.path.join(fl_dir, f"{acc}.fasta.DP_FULLLENGTH_RENAMED.fasta"))
    if not flf:
        continue
    fl = read_fasta(flf[0])
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    h = rows[0]
    with open(os.path.join(out, os.path.basename(flf[0])), "w") as o:
        for r in rows[1:]:
            d = dict(zip(h, r))
            if d.get("quality") != "intact" or d["TE_ID"] not in fl:
                continue
            try:
                l5, l3 = int(d["LTR5_length"]), int(d["LTR3_length"])
            except ValueError:
                skipped += 1
                continue
            s = fl[d["TE_ID"]]
            seg = s[l5:len(s) - l3]
            if len(seg) >= 300:
                o.write(f">{d['TE_ID']}\n{seg}\n"); n += 1
            else:
                skipped += 1
print(f"internal regions: {n} (skipped {skipped})")
