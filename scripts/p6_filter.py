#!/usr/bin/env python3
"""Phase 6 — artefact filtering of Phase-4 candidates.  [SCAFFOLD — not run yet]

For every candidate element, check and record (pass/fail + reason):
  1. Same element?  All domains within one annotated element's coordinates, correct
     strand, Ty3/Gypsy order (GAG-PROT-RT-RH-INT). Domains kb apart with an LTR
     between them = two nested elements, not a chimera.
  2. Nested insertion?  Intersect with the full TE annotation. Athila sits in
     pericentromeric heterochromatin where nesting is the norm — the largest
     false-positive source.
  3. Assembly quality.  Flag candidates in collapsed/low-coverage regions, near
     contig ends, or in poorly assembled accessions.
  4. Element integrity.  ORF intact (no frameshift/stop breaking the polyprotein);
     LTR-LTR identity + TSD presence as independent evidence of one insertion.
  5. Pipeline artefacts.  Did Athilafinder merge adjacent elements?

Inputs available on CSD3 (verified in Phase 0):
  * `*.DP_FULLLENGTH_RENAMED.fasta` (parent dir) — full-length element nucleotides.
  * `*.cls.lib` — per-element classified nucleotide library.
  * `*.dom.gff3` — domain coordinates/strand for the within-element order check.
  * (locate) dedicated LTR coordinates in the Athilafinder output for check 4.

Output: candidates_filtered.tsv with a column per check and a short reason string.
Expect most candidates to die here — that is the correct outcome.
"""

import sys

if __name__ == "__main__":
    sys.exit("p6_filter.py is a scaffold; implement after the Phase-4 checkpoint.")
