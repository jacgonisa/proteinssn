"""Athila domain-wise SSN pipeline — Phases 0-5 (6-8 are scaffolds).

Run on CSD3 inside the `athila_ssn` conda env:
    snakemake -s Snakefile --configfile config.yaml -c16
Dry run:
    snakemake -s Snakefile --configfile config.yaml -n
"""

import os

configfile: "config.yaml"

OUT      = config["out_dir"]
IN       = config["input_dir"]
LAYERS   = config["layers"]
SCRIPTS  = os.path.join(workflow.basedir, "scripts")
PY       = "python"

INV   = os.path.join(OUT, "inventory")
P2    = os.path.join(OUT, "phase2")
P3    = os.path.join(OUT, "phase3")
P4    = os.path.join(OUT, "phase4")
P5    = os.path.join(OUT, "phase5")


rule all:
    input:
        os.path.join(INV, "summary.md"),
        os.path.join(INV, "clade_disagreement.md"),
        os.path.join(P2,  "analysis_nodes.tsv"),
        expand(os.path.join(P3, "{layer}.edges.tsv"), layer=LAYERS),
        os.path.join(P4, "candidates.tsv"),
        os.path.join(P5, "phase5_summary.md"),


# ---- Phase 0 --------------------------------------------------------------- #
rule p0_inventory:
    output:
        summary = os.path.join(INV, "summary.md"),
        table   = os.path.join(INV, "elements.tsv"),
    params:
        counts = os.path.join(config["parent_dir"], "Athila_counts.tsv"),
    shell:
        "{PY} {SCRIPTS}/p0_inventory.py --input-dir {IN} --out-dir {OUT} "
        "--counts {params.counts}"


# ---- Phase 1 --------------------------------------------------------------- #
rule p1_clade:
    output:
        os.path.join(INV, "clade_disagreement.md"),
    shell:
        "{PY} {SCRIPTS}/p1_clade_disagree.py --input-dir {IN} --out-dir {OUT}"


# ---- Phase 2 --------------------------------------------------------------- #
rule p2_build:
    output:
        nodes   = os.path.join(P2, "analysis_nodes.tsv"),
        table   = os.path.join(P2, "node_table.tsv"),
        primary = expand(os.path.join(P2, "layers", "{layer}.primary.faa"),
                         layer=LAYERS),
    params:
        lf   = config["length_frac_min"],
        mid  = config["mmseqs_min_id"],
        mcov = config["mmseqs_min_cov"],
        cmode = config["mmseqs_cov_mode"],
    threads: config["threads"]
    shell:
        "{PY} {SCRIPTS}/p2_build_sets.py --input-dir {IN} --out-dir {OUT} "
        "--length-frac-min {params.lf} --min-id {params.mid} "
        "--min-cov {params.mcov} --cov-mode {params.cmode} --threads {threads}"


# ---- Phase 3 (primary: one network per layer on the full-5 node set) -------- #
rule p3_layer:
    input:
        faa = os.path.join(P2, "layers", "{layer}.primary.faa"),
    output:
        edges = os.path.join(P3, "{layer}.edges.tsv"),
    params:
        prefix = os.path.join(P3, "{layer}"),
        sens   = config["diamond_sensitivity"],
        mts    = config["diamond_max_target_seqs"],
        ev     = config["diamond_evalue"],
    threads: config["threads"]
    shell:
        "{PY} {SCRIPTS}/p3_allvsall.py --faa {input.faa} "
        "--out-prefix {params.prefix} --sensitivity {params.sens} "
        "--max-target-seqs {params.mts} --evalue {params.ev} --threads {threads}"


# ---- Phase 4 --------------------------------------------------------------- #
rule p4_concordance:
    input:
        edges = expand(os.path.join(P3, "{layer}.edges.tsv"), layer=LAYERS),
        nodes = os.path.join(P2, "analysis_nodes.tsv"),
    output:
        os.path.join(P4, "candidates.tsv"),
    params:
        k     = config["knn_k"],
        perms = config["null_permutations"],
    shell:
        "{PY} {SCRIPTS}/p4_concordance.py --edges-dir {P3} "
        "--layers {LAYERS} --nodes {input.nodes} --out-dir {OUT} "
        "--k {params.k} --perms {params.perms}"


# ---- Phase 5 --------------------------------------------------------------- #
rule p5_partition:
    input:
        edges = expand(os.path.join(P3, "{layer}.edges.tsv"), layer=LAYERS),
        nodes = os.path.join(P2, "analysis_nodes.tsv"),
    output:
        os.path.join(P5, "phase5_summary.md"),
    params:
        res = config["leiden_resolution"],
    shell:
        "{PY} {SCRIPTS}/p5_partition.py --edges-dir {P3} --layers {LAYERS} "
        "--nodes {input.nodes} --out-dir {OUT} --resolution {params.res}"
