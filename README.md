# proteinssn

Build and visualise a **protein Sequence Similarity Network (SSN)** from a FASTA
file with a single command. It runs an all-vs-all [DIAMOND](https://github.com/bbuchfink/diamond)
search, keeps the alignments that pass your identity / coverage / e-value
thresholds, groups proteins into clusters (connected components), and draws the
network for you.

```
proteins.faa  ──►  DIAMOND all-vs-all  ──►  filter edges  ──►  clusters  ──►  plots
```

![Example SSN](docs/example_network.png)

*Example: 60 laccase proteins at 30% identity / 50% coverage. Nodes are proteins,
coloured by cluster; two main families separate out, with a few singletons.*

## Install

DIAMOND is the only external dependency. The easiest way to get everything:

```bash
git clone https://github.com/jacgonisa/proteinssn
cd proteinssn
conda env create -f environment.yml     # installs diamond + python deps
conda activate proteinssn
pip install -e .
```

Already have DIAMOND on your `PATH`? Then just:

```bash
pip install git+https://github.com/jacgonisa/proteinssn
```

## Quick start

```bash
ssn build examples/laccase_demo.faa -i 30 -c 50
```

That's it. You get a folder `laccase_demo_ssn/` with:

| file | what it is |
|------|-----------|
| `*.edges.tsv`     | the network: `source  target  %identity  e-value` |
| `*.clusters.tsv`  | which cluster each protein belongs to |
| `*.network.png`   | static figure, nodes coloured by cluster |
| `*.network.html`  | **interactive** network — open in a browser, pan/zoom/hover |
| `*.diamond.tsv`   | raw DIAMOND hits (reuse it, see below) |

## The four commands

```bash
ssn build   proteins.faa            # full pipeline: search -> network -> plots
ssn plot    proteins.edges.tsv      # re-draw an existing edge list (no re-search)
ssn stats   proteins.clusters.tsv   # print the cluster size distribution
ssn sweep   proteins.faa            # try a range of identity thresholds
```

Run `ssn <command> -h` for all options. The ones you'll actually use:

- `-i, --identity` — minimum % identity for two proteins to be connected (default 30)
- `-c, --coverage` — minimum % of *both* sequences covered by the alignment (off by default)
- `-e, --evalue` — maximum e-value for an edge (off by default)
- `-t, --threads` — DIAMOND threads
- `-o, --outdir` — where to write results
- `--viz {both,static,interactive}` — which plots to make

### Choosing a threshold

Not sure what identity cutoff separates your protein families? Sweep it:

```bash
ssn sweep examples/laccase_demo.faa --start 20 --stop 90 --step 10
```

You get a table + plot of how the number of clusters and the largest cluster
change with the threshold — the "elbow" is usually a sensible cutoff. If you
already ran `ssn build`, reuse its DIAMOND table so nothing is recomputed:

```bash
ssn sweep --table laccase_demo_ssn/laccase_demo.diamond.tsv
```

### Re-plotting without re-searching

The all-vs-all search is the slow part. Once you have an `.edges.tsv` you can
recolour / re-render instantly:

```bash
ssn plot laccase_demo_ssn/laccase_demo.edges.tsv --viz interactive
```

## Using it in Cytoscape

Prefer to lay the network out yourself? Import `*.edges.tsv` as a network
(`source`/`target`, with `pident` as an edge weight) and `*.clusters.tsv` as a
node table to colour nodes by cluster.

## How clustering works

Proteins are nodes; a surviving DIAMOND hit is an edge. Clusters are the
**connected components** of that graph, labelled `cluster01`, `cluster02`, …
from largest to smallest. Proteins with no surviving edges become singleton
clusters. Raise `--identity` (and/or add `--coverage`) to split loosely related
families apart; lower it to merge them.

## Requirements

- Python ≥ 3.8
- DIAMOND ≥ 2.1
- networkx, pandas, matplotlib, pyvis (installed automatically)

## Acknowledgements & disclaimer

The core SSN logic here — the all-vs-all comparison, reciprocal-hit cleaning,
threshold-based edge filtering and connected-component clustering — is adapted
from **[SSNetworks](https://github.com/MiguelMSandin/SSNetworks)** by
**Miguel M. Sandin**. Some of the code was taken from that repository and then
adapted by me (Jacobo González) to:

- use **DIAMOND** instead of BLAST for the all-vs-all search (much faster on
  large protein sets),
- fit the **[MetEOr](https://github.com/jacgonisa/MetEOr)** workflow it grew out
  of, and
- be a **friendly, self-contained CLI** with built-in static + interactive
  visualisation and a threshold sweep.

Please cite the original work if you use this tool:

> Miguel M. Sandin. *SSNetworks: A brief, quick and dirty introduction to
> Sequence Similarity Networks.* GitHub repository.
> https://github.com/MiguelMSandin/SSNetworks

Note that the upstream SSNetworks repository is released under **GPL-3.0**.

## References

- Sandin, M. M. *SSNetworks.* https://github.com/MiguelMSandin/SSNetworks
- Buchfink, B., Reuter, K. & Drost, H-G. (2021) Sensitive protein alignments at
  tree-of-life scale using DIAMOND. *Nature Methods* 18, 366–368.
- Hagberg, A., Schult, D. & Swart, P. (2008) Exploring network structure,
  dynamics, and function using NetworkX. *Proc. 7th Python in Science Conf.*
