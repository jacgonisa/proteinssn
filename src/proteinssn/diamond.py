"""Thin wrapper around the DIAMOND aligner (makedb + all-vs-all blastp)."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from .network import DIAMOND_COLUMNS

SENSITIVITY_FLAGS = {
    "fast": [],
    "sensitive": ["--sensitive"],
    "more-sensitive": ["--more-sensitive"],
    "very-sensitive": ["--very-sensitive"],
    "ultra-sensitive": ["--ultra-sensitive"],
}


def find_diamond() -> str:
    """Return the DIAMOND executable path or exit with a helpful message."""
    exe = shutil.which("diamond")
    if exe:
        return exe
    sys.exit(
        "error: 'diamond' was not found on your PATH.\n"
        "  Install it with either of:\n"
        "    conda install -c bioconda diamond\n"
        "    mamba install -c bioconda diamond\n"
        "  or activate an environment that already has it, then re-run."
    )


def count_sequences(fasta: str) -> int:
    n = 0
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
    return n


def fasta_ids(fasta: str) -> list[str]:
    """First whitespace-delimited token of every FASTA header."""
    ids = []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def _run(cmd: list[str], quiet: bool) -> None:
    if not quiet:
        print("  $ " + " ".join(cmd), file=sys.stderr)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        sys.exit(f"error: DIAMOND failed (exit {exc.returncode}).")


def all_vs_all(
    fasta: str,
    workdir: Path,
    threads: int = 4,
    sensitivity: str = "sensitive",
    max_target_seqs: int = 3000,
    evalue: float = 1e-5,
    iterate: bool = True,
    quiet: bool = False,
) -> Path:
    """Build a DIAMOND db and run an all-vs-all blastp. Returns the table path."""
    diamond = find_diamond()
    workdir.mkdir(parents=True, exist_ok=True)
    stem = Path(fasta).stem
    db = workdir / f"{stem}.dmnd"
    table = workdir / f"{stem}.diamond.tsv"
    tmpdir = workdir / "tmp"
    tmpdir.mkdir(exist_ok=True)

    if not quiet:
        print("Step 1/2  building DIAMOND database", file=sys.stderr)
    _run([diamond, "makedb", "--in", fasta, "--db", str(db), "--quiet"], quiet)

    if not quiet:
        print("Step 2/2  all-vs-all blastp", file=sys.stderr)
    cmd = [
        diamond, "blastp",
        "-q", fasta,
        "--db", str(db),
        "--out", str(table),
        "--outfmt", "6", *DIAMOND_COLUMNS,
        "--max-target-seqs", str(max_target_seqs),
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--tmpdir", str(tmpdir),
    ]
    cmd += SENSITIVITY_FLAGS.get(sensitivity, ["--sensitive"])
    if iterate:
        cmd.append("--iterate")
    if quiet:
        cmd.append("--quiet")
    _run(cmd, quiet)
    return table
