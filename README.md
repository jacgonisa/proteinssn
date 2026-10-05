# Athila domain-wise SSN analysis

Do Athila LTR-retrotransposon families hold together along the whole element, or do
different protein domains of the same element place it in different families?
Cross-layer discordance between per-domain sequence-similarity networks is the
candidate signal for recombination / template switching / chimeric assembly. This
is treated as a **multiplex network** problem (same nodes = elements, one layer per
domain), not network alignment.

Input: a TEsorter run per genome for 153 *A. thaliana* genomes (Athilafinder
full-length Athila elements), on CSD3.

## What the data looks like (verified)

- 19,592 Athila elements; 98.6% `Athila` clade.
- Domain occupancy: GAG 18,966 · PROT 17,780 · **RT 4,598 · RH 4,475 · INT 4,400**.
  Only **3,860** elements are "full-5" (all of GAG+PROT+RT+RH+INT).
- No `env`/ORF2 domain is emitted → layers = GAG, PROT, RT, RH, INT.
- Element IDs are reused across genomes → every id is namespaced with its accession
  (parsed from the filename) before pooling.

## Primary vs secondary analysis

- **Primary (5-layer multiplex):** the ~3,860 full-5 elements, deduplicated to
  representatives (mmseqs 95/80 on the concatenated peptide). This is where the
  chimera test lives.
- **Secondary (GAG–PROT 2-layer):** the full fragment set, a weaker screen.

## Pipeline (Phases 0–5 runnable; 6–8 scaffolded)

| phase | script | output |
|---|---|---|
| 0 inventory | `p0_inventory.py` | `inventory/elements.tsv`, `summary.md` |
| 1 clade screen | `p1_clade_disagree.py` | `inventory/clade_disagreement.md` |
| 2 sequence sets | `p2_build_sets.py` | `phase2/node_table.tsv`, `analysis_nodes.tsv`, per-layer FASTAs |
| 3 all-vs-all | `p3_allvsall.py` | `phase3/<LAYER>.edges.tsv` (bitscore-ratio normalised) |
| 4 concordance | `p4_concordance.py` | `phase4/candidates.tsv`, `concordance_scores.tsv` |
| 5 partitions | `p5_partition.py` | `phase5/` NMI/ARI, Sankey flows, multiplex switch-list |
| 6–8 | `p6_filter.py` / `p7_confirm.py` / `p8_pangenome.py` | scaffolds |

`scripts/lib_ids.py` is the single parser for all TEsorter formats.

## Run

```bash
mamba env create -f env/athila_ssn.yml
conda activate athila_ssn
# edit config.yaml paths if needed, then:
bash run_all.sh 16          # or: snakemake --configfile config.yaml -c16
```

Inputs are read-only (owned by another user); all outputs go to `out_dir`
(`/rds/user/jg2070/hpc-work/athila_ssn/results`), never in place.

## Ground rules

- Discordance is **not** called "recombination" anywhere — it is discordance until
  Phase 7 (domain-specific rate shifts / selection produce the same network pattern).
- Report Phase 0 + Phase 1 numbers before Phase 2; checkpoint again after Phase 4
  before the expensive Phase 7.
- Every file-format assumption is verified against a real file.

Reuses the [`proteinssn`](https://github.com/jacgonisa/proteinssn) DIAMOND wrapper.
