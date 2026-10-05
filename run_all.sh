#!/usr/bin/env bash
# Run the Athila SSN pipeline (Phases 0-5) on CSD3.
set -euo pipefail

cd "$(dirname "$0")"

# Activate the env (create it first: mamba env create -f env/athila_ssn.yml)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate athila_ssn

THREADS="${1:-16}"

echo "== dry run =="
snakemake -s Snakefile --configfile config.yaml -n

echo "== run =="
snakemake -s Snakefile --configfile config.yaml -c"${THREADS}" -p --rerun-incomplete
