#!/usr/bin/env python3
"""Phase 7 — phylogenetic confirmation of survivors.  [SCAFFOLD — not run yet]

Take the top ~20-50 Phase-6 survivors plus representatives of both implicated
lineages and confirm properly:
  * Per-domain alignments: `mafft --maxiterate 1000 --localpair` -> `IQ-TREE -m MFP -B 1000`.
  * Topological incongruence: AU test (`iqtree -z` constraint trees, `-au`).
  * Breakpoint positions on nucleotide alignments of full elements where parent
    lineages are alignable: `GARD` (HyPhy) and `RDP5`.

A convincing chimera needs all three: one domain confidently placed in a different
lineage, a breakpoint between domains, and an intact element structure.

Tools (install when needed): mafft, iqtree, hyphy (bioconda).
"""

import sys

if __name__ == "__main__":
    sys.exit("p7_confirm.py is a scaffold; implement after the Phase-6 checkpoint.")
