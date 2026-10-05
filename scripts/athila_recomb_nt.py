#!/usr/bin/env python3
"""Within-Athila recombination at NUCLEOTIDE resolution.

Protein domains are too conserved to tell ATHILA families apart (best/2nd family
bitscore ~1.05x). Nucleotide domains are not (INT margin median ~1.7x). So we:

  1. extract each element's core-domain NUCLEOTIDE sequence (dom.gff3 coords on the
     full-length element, reverse-complemented on the - strand);
  2. assign each domain to its best ATHILA family by nucleotide search against the
     TAIR12 ATHILA internal exemplars;
  3. flag elements whose LONG domains (GAG/RT/RH/INT) confidently map to DIFFERENT
     families = within-Athila inter-family recombinants.

Usage:
    athila_recomb_nt.py --data-dir DIR --fulllength DIR --exemplars EX.fasta \
        --out-dir OUT [--min-bits 60] [--margin 1.3] [--min-len 150] [--threads 8]
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import tempfile
from collections import Counter, defaultdict

import lib_ids as L

LONG = ["GAG", "RT", "RH", "INT"]            # PROT too short/conserved even in nt
_TAB = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    d, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                d[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        d[name] = "".join(buf)
    return d


def revcomp(s):
    return s.translate(_TAB)[::-1]


def family_of(header):
    """ATHILA6a_I_nonauto#.. -> ATHILA6a ; ATHILA4c.1_I -> ATHILA4c."""
    h = header.split("#")[0]
    m = re.match(r"(ATHILA\d+[abc]?)", h)
    return m.group(1) if m else h


def build_refs(exemplars, workdir):
    ref = os.path.join(workdir, "refs.fa")
    n = 0
    with open(ref, "w") as o:
        for h, s in read_fasta(exemplars).items():
            if "_I" in h and "_LTR" not in h:      # internal coding sequences only
                o.write(f">{h}\n{s}\n")
                n += 1
    return ref, n


def extract_domains(data_dir, fl_dir, workdir, min_len):
    """Write query FASTA of every Athila element's core-domain NT."""
    q = os.path.join(workdir, "query.fa")
    n_el = n_dom = 0
    with open(q, "w") as out:
        for gff in sorted(glob.glob(os.path.join(data_dir, "*.dom.gff3"))):
            acc = L.accession_from_path(gff)
            flf = glob.glob(os.path.join(fl_dir, f"{acc}*DP_FULLLENGTH_RENAMED.fasta"))
            if not flf:
                continue
            seqs = read_fasta(flf[0])
            seen = set()
            for line in open(gff):
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                dom = attrs.get("gene")
                if dom not in L.CORE_DOMAINS:
                    continue
                el = f[0]
                if el not in seqs:
                    continue
                s, e, strand = int(f[3]), int(f[4]), f[6]
                sub = seqs[el][s - 1:e]
                if strand == "-":
                    sub = revcomp(sub)
                if len(sub) < min_len:
                    continue
                ns = L.nsid(acc, el)
                key = (ns, dom)
                if key in seen:          # keep the first (longest-listed) per domain
                    continue
                seen.add(key)
                out.write(f">{ns}@@{dom}\n{sub}\n")
                n_dom += 1
            n_el += len(set(k[0] for k in seen))
    print(f"[nt] domain NT queries: {n_dom}")
    return q


def search(query, ref, workdir, threads):
    res = os.path.join(workdir, "res.m8")
    subprocess.run(["mmseqs", "easy-search", query, ref, res,
                    os.path.join(workdir, "mm"), "--search-type", "3",
                    "-e", "1e-3", "--max-seqs", "50", "--threads", str(threads),
                    "--format-output", "query,target,bits", "-v", "1"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = defaultdict(dict)   # query -> {family: best_bits}
    with open(res) as fh:
        for line in fh:
            q, t, b = line.rstrip("\n").split("\t")
            fam = family_of(t)
            b = float(b)
            if b > best[q].get(fam, 0):
                best[q][fam] = b
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--fulllength", required=True)
    ap.add_argument("--exemplars", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-bits", type=float, default=60.0)
    ap.add_argument("--margin", type=float, default=1.3)
    ap.add_argument("--min-len", type=int, default=150)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    out = os.path.join(args.out_dir, "athila_recomb_nt")
    os.makedirs(out, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        ref, nref = build_refs(args.exemplars, tmp)
        print(f"[nt] family internal refs: {nref}")
        query = extract_domains(args.data_dir, args.fulllength, tmp, args.min_len)
        best = search(query, ref, tmp, args.threads)

    # assign family per (element, domain) with confidence
    assign = {}
    for q, fams in best.items():
        ns, dom = q.split("@@")
        ranked = sorted(fams.items(), key=lambda kv: -kv[1])
        fam, b = ranked[0]
        b2 = ranked[1][1] if len(ranked) > 1 else 0.0
        conf = b >= args.min_bits and (b2 == 0 or b / b2 >= args.margin)
        assign[(ns, dom)] = (fam, b, b2, conf)

    with open(os.path.join(out, "domain_family_nt.tsv"), "w") as f:
        f.write("ns_id\tdomain\tfamily\tbits\tsecond_bits\tmargin\tconfident\n")
        for (ns, dom), (fam, b, b2, conf) in sorted(assign.items()):
            mg = b / b2 if b2 else 999
            f.write(f"{ns}\t{dom}\t{fam}\t{b:.0f}\t{b2:.0f}\t{mg:.2f}\t{int(conf)}\n")

    by_elem = defaultdict(dict)
    for (ns, dom), (fam, b, b2, conf) in assign.items():
        if conf and dom in LONG:
            by_elem[ns][dom] = fam

    recs = []
    for ns, dm in by_elem.items():
        if len(dm) >= 2 and len(set(dm.values())) >= 2:
            recs.append((ns, dm))
    tested = sum(1 for d in by_elem.values() if len(d) >= 2)
    print(f"[nt] elements with >=2 confident LONG domains: {tested}")
    print(f"[nt] within-Athila inter-family recombinants: {len(recs)}")

    odd = Counter()
    pair = Counter()
    with open(os.path.join(out, "recombinants_nt.tsv"), "w") as f:
        f.write("ns_id\taccession\tchrom\tstart\tend\tn_long_domains\t"
                "architecture\tfamilies\todd_domain\n")
        for ns, dm in sorted(recs, key=lambda x: -len(x[1])):
            acc, el = L.split_nsid(ns)
            c = L.element_coords(el) or {}
            famc = Counter(dm.values())
            maj = famc.most_common(1)[0][0]
            odds = [d for d in LONG if d in dm and dm[d] != maj]
            for d in odds:
                odd[d] += 1
            fs = sorted(set(dm.values()))
            for i in range(len(fs)):
                for j in range(i + 1, len(fs)):
                    pair[(fs[i], fs[j])] += 1
            arch = ";".join(f"{d}:{dm[d]}" for d in LONG if d in dm)
            f.write(f"{ns}\t{acc}\t{c.get('chrom','')}\t{c.get('start','')}\t"
                    f"{c.get('end','')}\t{len(dm)}\t{arch}\t{','.join(fs)}\t"
                    f"{','.join(odds) or 'NA'}\n")

    with open(os.path.join(out, "summary_nt.md"), "w") as f:
        f.write("# Within-Athila recombination (nucleotide per-domain)\n\n")
        f.write(f"- elements with >=2 confident long domains: {tested}\n")
        f.write(f"- **inter-family recombinants: {len(recs)}**\n\n")
        f.write("## Odd-domain tally\n\n")
        for d in LONG:
            f.write(f"- {d}: {odd.get(d,0)}\n")
        f.write("\n## Top recombining family pairs\n\n| A | B | n |\n|---|---|---:|\n")
        for (a, b), n in pair.most_common(15):
            f.write(f"| {a} | {b} | {n} |\n")

    print(f"[nt] odd-domain tally: {dict(odd)}")
    print(f"[nt] top pairs: {pair.most_common(5)}")
    print(f"[nt] wrote {out}/recombinants_nt.tsv, summary_nt.md")


if __name__ == "__main__":
    main()
