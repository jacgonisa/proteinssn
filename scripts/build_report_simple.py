#!/usr/bin/env python3
"""Short HTML report focused on crossovers between ATHILA families (k-mer method).

Usage: build_report_simple.py RESULTS_DIR OUT.html
"""
import base64
import io
import sys
from collections import Counter

from PIL import Image

R, OUT = sys.argv[1:3]
C = f"{R}/crossovers"


def img(path, alt, maxw=1500):
    im = Image.open(path)
    if im.width > maxw:
        im = im.resize((maxw, int(im.height * maxw / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.convert("RGB").save(buf, "PNG", optimize=True)
    return f'<img src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}" alt="{alt}">'


def fig(path, alt, cap):
    return f'<figure><div class="plate">{img(path, alt)}</div><figcaption>{cap}</figcaption></figure>'


def table(head, rows, num=()):
    th = "".join(f'<th class="{"n" if i in num else ""}">{x}</th>' for i, x in enumerate(head))
    tb = "".join("<tr>" + "".join(f'<td class="{"n" if i in num else ""}">{x}</td>' for i, x in enumerate(r)) + "</tr>"
                 for r in rows)
    return f'<div class="tw"><table><thead><tr>{th}</tr></thead><tbody>{tb}</tbody></table></div>'


def tsv(p):
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    return rows[0], rows[1:]


h, calls = tsv(f"{C}/crossover_calls.tsv")
H = {c: i for i, c in enumerate(h)}
n_calls = len(calls)
n_cross = sum(r[H["class"]] == "crossover" for r in calls)
n_nest = n_calls - n_cross
flags = Counter(x for r in calls for x in r[H["nesting_flags"]].split(";") if x)
h, loci = tsv(f"{C}/crossover_loci.tsv")
rep = [r for r in loci if int(r[6]) >= 2]
n_rep_el = sum(int(r[5]) for r in rep)
h, pairs = tsv(f"{C}/crossover_pairs.tsv")

loci_rows = [[f"ATHILA{r[0].replace('-', ' × ')}", f"{r[7].replace('>', ' → ')}", r[1],
              f"{r[2]}–{r[3]}" if r[2] != r[3] else r[2], r[4], r[6]] for r in rep]
pair_rows = [[f"ATHILA{r[0].replace('-', ' × ')}", r[1], r[2], r[3], f"{r[4]}%",
              f"{r[5]} / {r[6]}"] for r in pairs]

CSS = """
:root{
  /* Layout: one reading column, tables and figures at full width. */
  --bg:#f5f6f4; --surface:#ffffff; --ink:#161a1a; --muted:#596061; --line:#dde1de;
  --accent:#0f6b61; --plate:#fcfcfb; --chip:#eef1ef;
  --display:"Newsreader", Georgia, serif; --body:"Public Sans", system-ui, sans-serif;
  --mono:"JetBrains Mono", ui-monospace, Menlo, monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --bg:#111413; --surface:#181c1b; --ink:#e6ebe9; --muted:#a3acaa; --line:#2b3230; --accent:#5cc2b4; --chip:#222826; color-scheme:dark;}}
:root[data-theme="dark"]{--bg:#111413; --surface:#181c1b; --ink:#e6ebe9; --muted:#a3acaa; --line:#2b3230; --accent:#5cc2b4; --chip:#222826; color-scheme:dark;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15.5px/1.62 var(--body);padding-inline:18px;padding-block:0 60px}
main{max-width:980px;margin:0 auto}
.col{max-width:70ch}
header{padding-block:48px 22px;border-bottom:1px solid var(--line)}
.eyebrow{font:600 11.5px/1 var(--mono);letter-spacing:.09em;text-transform:uppercase;color:var(--accent)}
h1{font:600 clamp(30px,4.6vw,44px)/1.1 var(--display);margin:12px 0;text-wrap:balance}
h2{font:600 25px/1.2 var(--display);margin:0 0 6px;text-wrap:balance}
.lede{font-size:17.5px;color:var(--muted);margin:0;max-width:66ch}
section{padding-block:34px 4px;border-top:1px solid var(--line);margin-top:24px}
section:first-of-type{border-top:0}
p{margin:10px 0} ul{padding-left:20px} li{margin:4px 0}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin:18px 0}
.kv div{background:var(--surface);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.kv .v{font:600 26px/1.1 var(--display);font-variant-numeric:tabular-nums}
.kv .k{font-size:12.5px;color:var(--muted);margin-top:4px}
figure{margin:18px 0}
.plate{background:var(--plate);border:1px solid var(--line);border-radius:10px;padding:10px;overflow-x:auto}
.plate img{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:13.5px;color:var(--muted);margin-top:8px}
.tw{overflow-x:auto;margin:14px 0;border:1px solid var(--line);border-radius:10px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{padding:7px 12px;border-bottom:1px solid var(--line);text-align:left}
th{font:600 11.5px/1.3 var(--mono);text-transform:uppercase;letter-spacing:.05em;color:var(--muted);background:var(--chip)}
tbody tr:last-child td{border-bottom:0}
td.n,th.n{text-align:right;font-family:var(--mono);font-size:12.5px;font-variant-numeric:tabular-nums}
.note{border-left:3px solid var(--accent);padding:8px 14px;background:var(--surface);border-radius:0 8px 8px 0;font-size:14.5px;margin:14px 0}
code{font:13px var(--mono);background:var(--chip);padding:1px 5px;border-radius:4px}
footer{margin-top:36px;padding-top:16px;border-top:1px solid var(--line);font-size:13px;color:var(--muted)}
"""

html = f"""<title>ATHILA Crossovers</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Newsreader:opsz,wght@6..72,600&family=Public+Sans:wght@400;600;700&display=swap">
<style>{CSS}</style>
<main>
<header>
<div class="eyebrow">Athilafinder intact elements · 153 Arabidopsis thaliana genomes · all ATHILA families</div>
<h1>Crossovers between ATHILA families</h1>
<p class="lede">Which ATHILA elements carry the 5′ part of one family and the 3′ part of another, after excluding nested insertions? Every intact element called by Athilafinder was scanned with family-specific k-mers, whatever domains it carries.</p>
</header>

<section><div class="col">
<h2>Result</h2></div>
<div class="kv">
<div><div class="v">20,248</div><div class="k">intact elements scanned (all families)</div></div>
<div><div class="v">{n_calls}</div><div class="k">elements with a family switch</div></div>
<div><div class="v">{n_nest}</div><div class="k">of these flagged as nested or merged</div></div>
<div><div class="v">{n_cross}</div><div class="k">crossover elements</div></div>
<div><div class="v">{len(loci)}</div><div class="k">crossover loci</div></div>
<div><div class="v">{len(rep)}</div><div class="k">loci seen in ≥3 accessions ({n_rep_el} elements)</div></div>
</div>
<div class="col"><p>Crossovers are rare (about 1% of intact elements) and mostly inherited: a handful of events, each present at the same place in many genomes. The largest is ATHILA2 × 4c on Chr3 (76 accessions); an independent ATHILA2 × 4c crossover with a different breakpoint sits on Chr1 (21 accessions).</p></div>
</section>

<section><div class="col"><h2>Crossover loci seen in several genomes</h2>
<p>Calls with the same family pair, chromosome and breakpoint interval, within 1 Mb of each other across accessions, are treated as one locus. Each of these carries one element per accession.</p></div>
{table(["family pair", "5′ → 3′", "chromosome", "position (Mb)", "switch", "accessions"], loci_rows, num=(5,))}
</section>

<section><div class="col"><h2>By family pair</h2></div>
{fig(f"{C}/fig_pairs.png", "Crossover loci per family pair", "Grey: all loci. Blue: loci seen in two or more accessions.")}
{table(["family pair", "elements", "loci", "replicated loci", "family homology", "centromere class"], pair_rows, num=(1, 2, 3, 4))}
<div class="col"><p>Family homology is the identity of the two TAIR12 consensus sequences. Centromere classes come from the CEN178 analysis (each family's share of copies inside the centromeric satellite core, tested across accessions).</p></div>
</section>

<section><div class="col"><h2>Where the switch falls</h2></div>
{fig(f"{C}/fig_breakpoints.png", "Breakpoint position relative to domains", "Breakpoints placed on each element's own TEsorter domains. Outside annotated domains usually means the 3′ region after INT, or elements whose domains are not all annotated.")}
</section>

<section><div class="col"><h2>How it was done</h2>
<ul>
<li><b>Region.</b> The internal region of each intact element, between the 5′ and 3′ LTRs given by Athilafinder.</li>
<li><b>Family markers.</b> 31-mers present in at least 10% of a family's clean elements and in less than 5% of every other family's (KMC). ATHILA7, 8a and 8b have no clean elements, so their markers come from the TAIR12 consensus.</li>
<li><b>Switch test.</b> Marker hits sharing one variant are merged into independent sites. A changepoint likelihood ratio compares "one family throughout" with "family A, then family B". The threshold was set on held-out clean elements.</li>
<li><b>Calibration.</b> At a 0.5% false-positive rate on held-out clean elements, the test recovers 63% of synthetic crossovers (1% → 75%, 0.1% → 13%). A simple run-length rule reached only 37% at the same false-positive rate.</li>
<li><b>Nesting.</b> An element is set aside as nested or merged if its internal region contains a block of LTR k-mers (from Athilafinder solo LTRs, absent from internal regions; {flags.get('ltr_inside', 0)} elements) or if its TEsorter domains are duplicated, on both strands or out of order ({flags.get('duplicated', 0) + flags.get('mixed_strand', 0) + flags.get('out_of_order', 0)} elements).</li>
</ul></div></section>

<section><div class="col"><h2>Caveats</h2>
<ul>
<li>Replication across accessions rules out random noise. It cannot rule out an inherited element whose two halves were misassigned the same way in every genome.</li>
<li>Within internal regions ATHILA4c is only weakly distinct from ATHILA4 and 4a, so crossovers involving 4c rest on few markers.</li>
<li>Recovery is about 60%, so some crossovers are missed, especially between near-identical families such as ATHILA6a and 6b.</li>
<li>None of these calls has yet been confirmed by phylogenetic incongruence or nucleotide breakpoint analysis.</li>
</ul></div></section>
<footer>Code: <code>crossover_scan.py</code>, <code>calibrate_markers.py</code>, <code>crossover_summary.py</code> on the <code>athila-analysis</code> branch of github.com/jacgonisa/proteinssn. Analysis by Jacob Gonzalez with Claude.</footer>
</main>"""
open(OUT, "w").write(html)
print("wrote", OUT)
