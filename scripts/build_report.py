#!/usr/bin/env python3
"""Assemble the self-contained HTML report from result files and figures.

Usage: build_report.py RESULTS_DIR EXAMPLES_DIR OUT.html
"""
import base64
import io
import sys
from collections import Counter

from PIL import Image

R, EX, OUT = sys.argv[1:4]


def img(path, alt, maxw=1700):
    im = Image.open(path)
    if im.width > maxw:
        im = im.resize((maxw, int(im.height * maxw / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.convert("RGB").save(buf, "PNG", optimize=True)
    b = base64.b64encode(buf.getvalue()).decode()
    return f'<img src="data:image/png;base64,{b}" alt="{alt}" loading="lazy">'


def fig(path, alt, caption, wide=True):
    cls = "fig wide" if wide else "fig"
    return f'<figure class="{cls}"><div class="plate">{img(path, alt)}</div><figcaption>{caption}</figcaption></figure>'


def tsv(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return rows[0], rows[1:]


def table(head, rows, num=(), caption=None):
    h = "".join(f'<th class="{"n" if i in num else ""}">{c}</th>' for i, c in enumerate(head))
    b = "".join("<tr>" + "".join(f'<td class="{"n" if i in num else ""}">{c}</td>'
                                 for i, c in enumerate(r)) + "</tr>" for r in rows)
    cap = f"<caption>{caption}</caption>" if caption else ""
    return f'<div class="tw"><table>{cap}<thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def pill(cls):
    return f'<span class="pill {cls}">{cls}</span>'


# ---------------------------------------------------------------- tables ---- #
# family centro classes
h, rows = tsv(f"{R}/centro/family_centro_class.tsv")
cen_rows = [[f"ATHILA{r[0]}", pill(r[1]), r[2], f"{float(r[3])*100:.0f}%", f"{float(r[4])*100:.0f}%",
             f"{float(r[5]):+.2f}", f"{r[7]} / {r[8]}"] for r in rows]
cen_table = table(["family", "class", "copies", "centromeric", "arm", "median log₂ enrichment",
                   "accessions enriched / depleted"], cen_rows, num=(2, 3, 4, 5))

# permutation (loci level)
h, rows = tsv(f"{R}/perm_pairs/pair_permutation.tsv")
H = {c: i for i, c in enumerate(h)}
perm = []
for r in rows:
    if r[0] != "loci":
        continue
    tags = []
    if float(r[H["q_enrich_abundance"]]) < 0.05: tags.append("enriched vs abundance")
    if float(r[H["q_enrich_degree"]]) < 0.05: tags.append("enriched vs propensity")
    if float(r[H["q_deplete_abundance"]]) < 0.05: tags.append("depleted vs abundance")
    if float(r[H["q_deplete_degree"]]) < 0.05: tags.append("depleted vs propensity")
    if tags and (int(r[H["observed"]]) > 0 or "depleted vs propensity" in tags):
        perm.append([r[H["pair"]], r[H["observed"]], f'{float(r[H["oe_abundance"]]):.2g}',
                     f'{float(r[H["oe_degree"]]):.2g}', "; ".join(tags)])
perm.sort(key=lambda x: -int(x[1]))
perm_table = table(["family pair", "distinct loci", "O/E abundance", "O/E propensity", "significant (q < 0.05)"],
                   perm[:14], num=(1, 2, 3))

# robust k-mer mosaics per pair (internal span)
h, rows = tsv(f"{R}/bias_kmer/bias_kmer_pairs.tsv")
H = {c: i for i, c in enumerate(h)}
pairs = [[r[H["pair"]], f'{float(r[H["homology"]]):.0f}%', r[H["raw"]], r[H["loci"]]]
         for r in rows if int(r[H["raw"]]) > 0]
pairs_table = table(["family pair", "family homology", "mosaic elements", "distinct loci"], pairs, num=(1, 2, 3))

# inter-clade chimeras
h, rows = tsv(f"{EX}/recombination_events.tsv")
inter = Counter((r[9], r[8]) for r in rows if r[0] == "inter_clade")
inter_table = table(["partners", "foreign domain", "elements"],
                    [[p.replace("Athila+", "ATHILA + "), s.replace(";", " + "), str(n)]
                     for (p, s), n in inter.most_common()], num=(2,))

# census
h, rows = tsv(f"{R}/fragments/copies_with_context.tsv")
lab = Counter(r[6] for r in rows)
accs = len({r[1] for r in rows})

# ------------------------------------------------------------------- page --- #
CSS = """
:root{
  /* Layout: one 70ch reading column; figures and tables break out to 1100px. */
  --bg:#f5f6f4; --surface:#ffffff; --ink:#161a1a; --muted:#596061; --line:#dde1de;
  --accent:#0f6b61; --cen:#c43d36; --away:#2a6cc4; --plate:#fcfcfb; --chip:#eef1ef;
  --display:"Newsreader", Georgia, "Times New Roman", serif;
  --body:"Public Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono:"JetBrains Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#111413; --surface:#181c1b; --ink:#e6ebe9; --muted:#a3acaa; --line:#2b3230;
    --accent:#5cc2b4; --cen:#ef6b62; --away:#6ea7ef; --chip:#222826; color-scheme:dark;
  }
}
:root[data-theme="dark"]{
  --bg:#111413; --surface:#181c1b; --ink:#e6ebe9; --muted:#a3acaa; --line:#2b3230;
  --accent:#5cc2b4; --cen:#ef6b62; --away:#6ea7ef; --chip:#222826; color-scheme:dark;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15.5px/1.62 var(--body);
  padding-inline:18px;padding-block:0 64px}
main{max-width:1100px;margin:0 auto}
.col{max-width:70ch}
header{padding-block:56px 28px;border-bottom:1px solid var(--line);margin-bottom:8px}
.eyebrow{font:600 11.5px/1 var(--mono);letter-spacing:.09em;text-transform:uppercase;color:var(--accent)}
h1{font:600 clamp(32px,5vw,50px)/1.08 var(--display);margin:14px 0 14px;text-wrap:balance;letter-spacing:-.01em}
h2{font:600 28px/1.2 var(--display);margin:0 0 6px;text-wrap:balance}
h3{font:600 16px/1.35 var(--body);margin:26px 0 6px}
.lede{font-size:18px;color:var(--muted);max-width:68ch;margin:0}
.meta{display:flex;flex-wrap:wrap;gap:8px 22px;margin-top:20px;font:12.5px/1.4 var(--mono);color:var(--muted)}
.meta b{color:var(--ink);font-weight:600}
section{padding-block:40px 8px;border-top:1px solid var(--line);margin-top:28px}
section:first-of-type{border-top:0}
.sec-num{font:500 12px/1 var(--mono);color:var(--muted);letter-spacing:.06em;display:block;margin-bottom:10px}
p{margin:10px 0}
a{color:var(--accent);text-underline-offset:2px}
a:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:2px}
code,.mono{font:13px/1.4 var(--mono);background:var(--chip);padding:1px 5px;border-radius:4px}
nav.toc{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:4px 26px;
  margin:26px 0 6px;padding:16px 18px;background:var(--surface);border:1px solid var(--line);border-radius:10px}
nav.toc a{text-decoration:none;color:var(--ink);font-size:14px;display:flex;gap:10px;padding:3px 0}
nav.toc a span{font:12px/1.6 var(--mono);color:var(--muted);min-width:22px}
nav.toc a:hover{color:var(--accent)}
.findings{list-style:none;padding:0;margin:22px 0 0;display:grid;gap:12px}
.findings li{padding:14px 16px;background:var(--surface);border:1px solid var(--line);border-radius:10px}
.findings b{display:block;margin-bottom:2px}
.fig{margin:22px 0 10px}
.plate{background:var(--plate);border:1px solid var(--line);border-radius:10px;padding:10px;overflow-x:auto}
.plate img{display:block;max-width:100%;height:auto;margin:0 auto}
figcaption{font-size:13.5px;color:var(--muted);margin-top:8px;max-width:90ch}
.tw{overflow-x:auto;margin:16px 0;border:1px solid var(--line);border-radius:10px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:13.5px}
caption{text-align:left;padding:10px 12px 0;color:var(--muted);font-size:13px}
th,td{padding:7px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{font:600 11.5px/1.3 var(--mono);text-transform:uppercase;letter-spacing:.05em;color:var(--muted);background:var(--chip)}
tbody tr:last-child td{border-bottom:0}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono);font-size:12.5px}
.pill{display:inline-block;font:600 11px/1 var(--mono);padding:4px 8px;border-radius:99px;white-space:nowrap;
  color:var(--ink);background:var(--chip);border:1px solid var(--line)}
.pill.centrophilic{color:#fff;background:var(--cen);border-color:var(--cen)}
.pill.centrophobic{color:#fff;background:var(--away);border-color:var(--away)}
.note{border-left:3px solid var(--accent);padding:8px 14px;margin:16px 0;background:var(--surface);border-radius:0 8px 8px 0;font-size:14.5px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:18px}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:18px 0}
.kv div{background:var(--surface);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.kv .v{font:600 24px/1.1 var(--display);font-variant-numeric:tabular-nums}
.kv .k{font-size:12.5px;color:var(--muted);margin-top:4px}
ul.plain{padding-left:20px}
ul.plain li{margin:4px 0}
footer{margin-top:40px;padding-top:18px;border-top:1px solid var(--line);font-size:13px;color:var(--muted)}
@media (max-width:520px){body{font-size:15px} header{padding-block:36px 20px}}
"""

S = []
S.append(f"""<title>ATHILA Recombination Atlas</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Public+Sans:wght@400;600;700&display=swap">
<style>{CSS}</style>
<main>
<header>
  <div class="eyebrow">Athila LTR retrotransposons · 153 Arabidopsis thaliana genomes</div>
  <h1>How ATHILA families recombine, and where</h1>
  <p class="lede">Do ATHILA domains stay with their family along the whole element, or do elements carry pieces from different families? This report collects every analysis run so far: detection, validation, the families that recombine, how that relates to homology and to the centromere, and what happens to the copies that are chopped up.</p>
  <div class="meta"><span><b>153</b> genomes</span><span><b>19,592</b> TEsorter-classified elements</span><span><b>{len(rows):,}</b> ATHILA-derived copies genome-wide</span><span><b>148</b> genomes with CEN178 annotation</span><span>Built 7 Oct 2026</span></div>
</header>

<nav class="toc" aria-label="Sections">
  <a href="#findings"><span>—</span>Key findings</a>
  <a href="#data"><span>1</span>Data and definitions</a>
  <a href="#detect"><span>2</span>Detecting recombinants</a>
  <a href="#inter"><span>3</span>Chimeras with other LTR clades</a>
  <a href="#within"><span>4</span>Recombination between ATHILA families</a>
  <a href="#type"><span>5</span>Crossover-like or conversion-like</a>
  <a href="#homology"><span>6</span>Recombination and homology</a>
  <a href="#partners"><span>7</span>Who recombines with whom</a>
  <a href="#centro"><span>8</span>Centrophilic and centrophobic families</a>
  <a href="#age"><span>9</span>Age: two regimes</a>
  <a href="#all"><span>10</span>All copies: chopping and junctions</a>
  <a href="#caveats"><span>11</span>Caveats and next steps</a>
  <a href="#files"><span>12</span>Files and code</a>
</nav>

<section id="findings"><div class="col">
<h2>Key findings</h2>
<ul class="findings">
<li><b>Recombination between ATHILA families is real but needs nucleotide-level markers to see.</b>Protein domains are too similar across families to assign them (best vs second-best family differ by ~5%). Family-specific k-mer markers resolve them and recover every planted and known case.</li>
<li><b>619 recombinant mosaics (3.3% of typed elements), almost all single-switch.</b>Crossover-like switches dominate; only two convincing conversion-like tracts were found, though short tracts are hard to detect.</li>
<li><b>More similar families recombine more.</b>Across all 91 family pairs, recombination rises with family homology after normalising for abundance and detectability (Spearman ρ = +0.35, permutation p ≈ 0.0015; odds ratio 4.5 per +10% homology).</li>
<li><b>Partner choice is not random.</b>ATHILA6a–6b, 1–6, 2–4c, 7–7a and 0–3 recombine more than expected; ATHILA1 and ATHILA2, the two most common families, recombine with each other four times less than abundance predicts.</li>
<li><b>Centrophilic families (ATHILA5, 6b, 1) never recombine with each other</b> (0 events vs ~10.6 expected); their recombination partners are neutral or centrophobic families.</li>
<li><b>Two regimes by age.</b>The youngest recombinants are ATHILA6a×6b chimeras with near-identical LTRs, half of them centromeric, consistent with recombination during retrotransposition. The oldest are ATHILA2×4c and 2×6b in pericentromeres.</li>
<li><b>Most ATHILA sequence is chopped.</b>A typical genome has 128 intact elements, 87 solo LTRs and about 400 internal fragments. Centrophobic families are 4–25 times more fragmented than centrophilic ones; copies inside centromeres are the least fragmented and form the fewest solo LTRs.</li>
<li><b>Deletion junctions show little microhomology.</b>Clean internal deletions are no more microhomology-rich than random breakpoints (a weak excess at ≥4 bp in fragments), unlike family switches, which sit in long shared sequence.</li>
</ul></div></section>
""")

S.append(f"""<section id="data"><span class="sec-num">1</span><div class="col">
<h2>Data and definitions</h2>
<p>Input is Athilafinder full-length ATHILA calls and TEsorter domain annotation for 153 accessions on the CSD3 RDS, plus Athilafinder's solo-LTR calls, the genome assemblies themselves, and TRASH CEN178 satellite annotation for 148 of them. ATHILA families follow the pipeline's own definitions: TAIR12 exemplar internal and LTR sequences for ATHILA0–9.</p>
<ul class="plain">
<li><b>Element IDs</b> are reused across genomes, so every ID is prefixed with its accession.</li>
<li><b>Full-5 elements</b> carry all of GAG, PROT, RT, RH and INT (3,860 of 19,592). After collapsing near-identical copies (95% identity, 80% coverage) this is about 510 distinct elements, the effective sample size for the domain-network analyses.</li>
<li><b>Centromeric context</b>: an element is <i>centromeric</i> if it lies inside the main CEN178 cluster of its own chromosome in its own genome (median core 3.0 Mb), <i>pericentromeric</i> within 2 Mb of it, otherwise <i>arm</i>. Red and blue mean the same thing in every figure: red is centromere-associated, blue is away from it.</li>
<li><b>Distinct loci</b>: the same insertion appears in many accessions. Wherever counts could be inflated by inheritance, they are also reported per distinct locus (same chromosome, start within 500 kb).</li>
</ul></div></section>
""")

S.append(f"""<section id="detect"><span class="sec-num">2</span><div class="col">
<h2>Detecting recombinants</h2>
<p>The first approach compared per-domain protein similarity networks. Every full-5 element was alignable to every other in every domain, so the networks were complete graphs, and a test with planted chimeras showed it could detect a domain from an unrelated sequence (3 of 5) but not one swapped from another ATHILA family (0 of 5). TEsorter calls all of them simply "Athila", so family resolution has to come from elsewhere.</p>
<p>Three refinements, each validated before use:</p>
<ul class="plain">
<li><b>Nucleotide family typing.</b> Each domain's DNA was matched to the TAIR12 family exemplars. For INT the best vs second-best family margin went from 1.05× (protein) to a median 1.68× (DNA).</li>
<li><b>Family-specific k-mer markers.</b> 31-mers present in at least 25% of a family's clean members and in no other family (built with KMC, subtracting all other families). Even sister families have thousands of distinguishing k-mers (ATHILA6a vs 6b: about 4,900 each way).</li>
<li><b>Scanning only the internal coding span</b> (GAG start to INT end), so LTR and env-region sequence cannot create false end switches. With these markers 3.3% of elements are mosaics, against 31% with exemplar-only markers.</li>
</ul>
<p>An earlier whole-domain method suggested recombination <i>decreases</i> with homology (ρ = −1). That was a detection artefact: it could not resolve swaps between similar families. Section 6 has the corrected result.</p>
</div></section>
""")

S.append(f"""<section id="inter"><span class="sec-num">3</span><div class="col">
<h2>Chimeras with other LTR clades</h2>
<p>33 elements classified as Athila carry a domain that TEsorter assigns to another Ty3/gypsy clade. The clearest is an Athila element whose GAG and PROT are Athila but whose integrase and RNase H-like domain are Retand, found at the same Chr5 region in 13 accessions. The Retand module sits on the opposite strand, so it may be a nested antisense insertion rather than a fused element.</p>
</div>
{inter_table}
{fig(f"{EX}/chimera_9994_Chr5.png", "Domain architecture and per-domain networks of an Athila element with a Retand integrase", "Accession 9994, Chr5:16,533,472. A: domain architecture coloured by clade; B: the element's integrase sits in the Retand cluster of the INT network; C: its GAG sits in the Athila cluster.")}
</section>
""")

S.append(f"""<section id="within"><span class="sec-num">4</span><div class="col">
<h2>Recombination between ATHILA families</h2>
<p>With robust k-mer markers, 619 of 18,610 typed elements are mosaics of two ATHILA families. They resolve into the pairs below. Recombination sits at the ends of the coding region (GAG or INT); the RT–RH core is never the swapped part.</p>
</div>
<div class="grid2">{pairs_table}<div class="col">
<h3>Two well-supported examples</h3>
<p><b>GAG from ATHILA2, pol from ATHILA4c.</b> 80 copies, all on Chr3 between 16.1 and 19.8 Mb, in different accessions. They collapse to 9 loci, so this is one or a few old events inherited many times.</p>
<p><b>Toufl-1 and Elh-2</b> are the only complete (full-5) recombinants: GAG from ATHILA7a, with PROT, RT, RH and INT all from ATHILA7.</p>
</div></div>
{fig(f"{EX}/recomb_Toufl1_full5_ALL.png", "Per-domain nucleotide networks of the Toufl-1 recombinant", "Toufl-1, Chr5:18,239,639. Nucleotide networks of all confidently typed copies (GAG 11,155, RT 4,118, INT 4,291 nodes). The starred element sits with ATHILA7a for GAG and with ATHILA7 for RT and INT.")}
{fig(f"{EX}/recomb_9994_allfamilies.png", "Per-domain nucleotide networks of a Chr3 ATHILA2/ATHILA4c recombinant", "9994, Chr3:17,724,292, one of the 80 Chr3 copies: GAG clusters with ATHILA2, RT and RH with ATHILA4c. Its INT is degraded, which is why it is not full-5.")}
</section>
""")

S.append(f"""<section id="type"><span class="sec-num">5</span><div class="col">
<h2>Crossover-like or conversion-like</h2>
<p>A mosaic can switch family once (A→B, crossover-like) or carry a tract of one family inside another (A→B→A, conversion-like). Counting independent variant sites rather than overlapping k-mers, and comparing with a shuffle of those sites:</p>
</div>
{table(["min. independent sites per run", "crossover-like (expected by chance)", "conversion-like (expected by chance)"],
       [["2", "297 (40)", "7 (61)"], ["3", "180 (13)", "2 (15)"], ["5", "37 (1.5)", "0 (0.6)"]], num=(0,))}
<div class="col"><p>Two conversion-like tracts are well supported, both in accession Had-6b on Chr4: a 1.3 kb ATHILA9 tract inside ATHILA4c (13 sites) and an 0.8 kb ATHILA1 tract inside ATHILA4c (4 sites). Planting synthetic tracts into clean elements shows the limit: 250–500 bp tracts between divergent families are recovered 20–60% of the time, 100 bp tracts almost never, and nothing between ATHILA6a and 6b. Within what can be seen, crossover-like recombination dominates.</p></div>
{table(["host > donor", "crossover recovered", "conversion 100 bp", "250 bp", "500 bp", "1 kb"],
       [["2 > 6b", "97%", "3%", "32%", "63%", "17%"], ["1 > 2", "85%", "12%", "30%", "33%", "8%"],
        ["2 > 4c", "88%", "0%", "3%", "10%", "17%"], ["6 > 1", "40%", "0%", "20%", "37%", "17%"],
        ["6a > 6b", "23%", "0%", "2%", "0%", "0%"]], num=(1, 2, 3, 4, 5),
       caption="Detection power from planted events (0% false positives on unmodified elements)")}
</section>
""")

S.append(f"""<section id="homology"><span class="sec-num">6</span><div class="col">
<h2>Recombination and homology</h2>
<p>All 91 family pairs were tested, including the 76 with no recombinants. Every normalisation gives the same positive association between family homology and recombination: raw counts, distinct loci, per opportunity (abundance of A × abundance of B), per element of the rarer family, and adjusted for how many distinguishing markers a pair has.</p>
</div>
<div class="kv"><div><div class="v">+0.35</div><div class="k">Spearman ρ, homology vs recombination (all pairs)</div></div>
<div><div class="v">0.0015</div><div class="k">permutation p (5,000 shuffles)</div></div>
<div><div class="v">4.5×</div><div class="k">odds of any recombination per +10% homology (95% CI 2.0–10)</div></div>
<div><div class="v">+0.27 to +0.43</div><div class="k">ρ when any one family is dropped</div></div></div>
{fig(f"{EX}/bias_kmer.png", "Recombination against family homology under four normalisations", "Each point is a family pair. Open grey circles: no recombinant found. The exception is ATHILA2–4c (53% homology), which reduces to 9 loci.")}
</section>
""")

S.append(f"""<section id="partners"><span class="sec-num">7</span><div class="col">
<h2>Who recombines with whom</h2>
<p>Family labels were permuted in two ways. The abundance null redraws partners in proportion to how common each family is. The propensity null keeps how often each family recombines but shuffles its partners. Counts are distinct loci.</p>
</div>{perm_table}</section>
""")

S.append(f"""<section id="centro"><span class="sec-num">8</span><div class="col">
<h2>Centrophilic and centrophobic families</h2>
<p>A family is centrophilic when its copies are more often inside the CEN178 core than ATHILA as a whole, consistently across accessions (sign test, each genome counted once). Overall 26% of ATHILA copies are centromeric, 65% pericentromeric and 8% on arms.</p>
</div>{cen_table}
{fig(f"{EX}/centro_network.png", "ATHILA family recombination network coloured by centromere preference", "Nodes are families (area ∝ copies; percentage = share inside the CEN178 core). Edges are distinct recombinant loci: black = more than expected, grey = as expected, dashed = fewer than expected. No edges join the three centrophilic families.")}
<div class="col"><p>Recombinant mosaics are slightly more common inside centromeres (3.6%) than in pericentromeres (3.2%) and on arms (1.8%); centromeric versus the rest gives p = 0.053.</p></div>
</section>
""")

S.append(f"""<section id="age"><span class="sec-num">9</span><div class="col">
<h2>Age: two regimes</h2>
<p>Athilafinder reports 5′–3′ LTR identity for intact elements. LTRs are identical when an element inserts, so identity measures time since insertion. Centromeric copies are much younger (median 98.5%) than pericentromeric ones (95.1%).</p>
</div>
{fig(f"{R}/report_figs/age.png", "Element age by context and recombinant frequency by age", "Left: LTR identity by context. Right: share of mosaics per age quintile.")}
<div class="col"><p>Mosaics are most common among the oldest and the youngest elements. The youngest group (LTRs at least 98.9% identical; 172 mosaics, 57 loci) is almost entirely ATHILA6a × 6b, half of it centromeric. Elements that young have had no time to recombine after inserting, which points to template switching during reverse transcription: they were born chimeric. The oldest group is ATHILA2 × 4c and 2 × 6b, nearly all pericentromeric.</p></div>
</section>
""")

S.append(f"""<section id="all"><span class="sec-num">10</span><div class="col">
<h2>All copies: chopping and junctions</h2>
<p>To include truncated copies, every chromosome of every genome was aligned to the ATHILA exemplars and the hits merged into copies. Copies matching an Athilafinder intact element or solo LTR keep that label; everything else is a fragment.</p>
</div>
<div class="kv"><div><div class="v">128</div><div class="k">intact elements per genome (median)</div></div>
<div><div class="v">87</div><div class="k">solo LTRs per genome</div></div>
<div><div class="v">397</div><div class="k">internal fragments per genome</div></div>
<div><div class="v">188</div><div class="k">LTR fragments per genome</div></div></div>
{fig(f"{R}/report_figs/fragments.png", "Internal fragments per intact element by family and by context", "Centrophobic families are the most fragmented, centrophilic the least. Copies inside centromeres are the least fragmented, partly because they are younger.")}
<div class="col"><h3>Solo LTRs: LTR–LTR recombination</h3>
<p>A solo LTR is left when the two LTRs of an element recombine and delete the internal region. Solo LTRs are rarer inside centromeres (32% of LTR-bearing copies) than in pericentromeres (45%) or on arms (49%); centromeric versus the rest gives an odds ratio of 0.55, p = 6 × 10⁻¹⁰⁸. The same holds within most families, for example ATHILA5 at 11% inside the centromere against 39% in the pericentromere. Centromeric copies are younger, so part of this is time; solo LTRs cannot be dated, so the two cannot be fully separated.</p></div>
{fig(f"{R}/report_figs/solo_ltr.png", "Solo-LTR fraction by LTR family and context", "LTR families are groups where TAIR12 exemplars share an LTR: ATHILA6 = 6/6a/6b, ATHILA3 = 0/3, ATHILA7 = 7/7a.")}
<div class="col"><h3>Junctions inside copies</h3>
<p>Each copy was split-aligned to the family consensus and every junction between consecutive pieces classified. Fragments carry about 18 times more inversions and 4 times more duplications and LTR family switches than intact elements.</p></div>
{fig(f"{R}/report_figs/junction_rates.png", "Rearrangement junctions per 100 copies in intact elements and fragments", "LTR–internal family mismatches are not shown: several TAIR12 LTR exemplars are shared between families, so they are not evidence of recombination.")}
<div class="col"><p>Microhomology was measured exactly from the consensus at each clean deletion and compared with random breakpoint pairs at the same distance. After collapsing deletions inherited across genomes (599 fragment deletions are 86 distinct junctions), there is no enrichment in intact elements and only a weak excess at ≥4 bp in fragments (5.8% vs 2.1%). By contrast, clean family switches have a median of 60 bp of shared sequence at the junction, as expected for homologous recombination.</p></div>
{fig(f"{R}/report_figs/microhomology.png", "Microhomology at clean deletion junctions against random breakpoints", "Unique deletion junctions only.")}
{fig(f"{R}/report_figs/truncation.png", "Positions of fragment truncation ends along the family consensus", "Fragment ends are spread almost uniformly along the consensus: a mild excess 3′ of INT and between GAG and PROT, fewer between PROT and RT. Domain positions are projected from TEsorter annotation of full-length elements.")}
</section>
""")

S.append("""<section id="caveats"><span class="sec-num">11</span><div class="col">
<h2>Caveats and next steps</h2>
<ul class="plain">
<li>All recombinants are candidates. None has yet been confirmed with per-domain phylogenies (AU test) or nucleotide breakpoint methods such as GARD or RDP.</li>
<li>Family homology is one consensus-alignment value per pair. Identity at each breakpoint would be a more direct predictor.</li>
<li>Sister families can share pericentromeric niches, which raises local opportunity independently of homology.</li>
<li>ATHILA7, 8a and 8b have no clean elements, so their markers come from exemplars only.</li>
<li>Young ATHILA6a×6b mosaics could also be an intermediate ATHILA6 lineage absent from the exemplars. Building a consensus for that group would separate the two.</li>
<li>Centromeric copies are younger, which partly explains their lower solo-LTR and fragment rates.</li>
<li>Copy discovery uses minimap2 at ≥70% identity, so very old, highly diverged remnants are missed.</li>
</ul>
<p>Useful next steps: phylogenetic confirmation of the 6a×6b and 2×4c recombinants; breakpoint-level homology instead of family-level; an ATHILA6 subfamily consensus; and a nested-insertion filter for the "deletion with insert" junctions.</p>
</div></section>

<section id="files"><span class="sec-num">12</span><div class="col">
<h2>Files and code</h2>
<p>Everything is on the <code>athila-analysis</code> branch of <code>github.com/jacgonisa/proteinssn</code>. Main scripts:</p>
<ul class="plain">
<li><code>kmc_markers.py</code>, <code>kmer_recomb.py</code>: family markers and mosaic scan</li>
<li><code>island_scan.py</code>, <code>conv_power.py</code>: crossover vs conversion and detection power</li>
<li><code>recomb_bias_kmer.py</code>, <code>perm_pairs.py</code>: homology association and partner permutation</li>
<li><code>centro_class.py</code>, <code>centro_network.py</code>: centromere classification and network</li>
<li><code>genome_copies.py</code>, <code>junction_scan.py</code>, <code>fragment_analysis.py</code>, <code>solo_ltr.py</code>, <code>age_effects.py</code>: all copies, junctions, solo LTRs and age</li>
</ul>
</div></section>
<footer>Analysis by Jacob Gonzalez with Claude. Data: Athilafinder and TEsorter outputs on the pangenome RDS; TRASH CEN178 annotation; TAIR12 ATHILA exemplars.</footer>
</main>""")

open(OUT, "w").write("\n".join(S))
print("wrote", OUT)
