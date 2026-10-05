"""Single source of truth for IDs and TEsorter file parsing.

Formats were verified against the real 153-genome TEsorter run; do not hard-code a
separator that has not been seen there. Key facts:

* One TEsorter run per genome. The genome ACCESSION is only in the filename
  (e.g. ``100620.Chr_scaffolds.fasta.cls.tsv`` -> ``100620.Chr_scaffolds``).
* Element IDs (``Atha_Chr1.13127937-13136618_P_Athila``) are reused verbatim across
  genomes, so they MUST be namespaced with the accession before pooling.
* ``cls.tsv``  cols: ``#TE Order Superfamily Clade Complete Strand Domains``
  where ``Domains`` is space-separated ``DOMAIN|CLADE`` tokens.
* ``dom.faa`` header:
  ``>ELEMENT|...:Ty3-GAG ID=...;gene=GAG;clade=Athila;coverage=;evalue=;probability=``
* ``dom.tsv``  cols: ``#id length evalue coverge probability score`` where
  ``id`` is ``ELEMENT|...:Ty3-DOMAIN``.
"""

from __future__ import annotations

import re
from pathlib import Path

# The five REXdb domains that make up the Athila (Ty3/Gypsy) polyprotein, in order.
CORE_DOMAINS = ["GAG", "PROT", "RT", "RH", "INT"]
POL_DOMAINS = ["RT", "RH", "INT"]

NS_SEP = "::"  # accession <-> element separator; absent from both sides

_FASTA_SUFFIX = re.compile(r"\.fasta.*$")
_DOMAIN_FROM_ID = re.compile(r":[^:|]*?-([A-Za-z0-9]+)$")  # ...:Ty3-GAG -> GAG
_ATTR = re.compile(r"(\w+)=([^;]+)")
# Atha_Chr1.13127937-13136618_P_Athila -> (Chr1, 13127937, 13136618, P, Athila)
_ELEM_COORDS = re.compile(r"\.?([A-Za-z0-9]+)\.(\d+)-(\d+)_([A-Za-z]+)_(.+)$")


def element_coords(element: str):
    """Parse chrom/start/end/tag/clade from an element id. Returns a dict or None."""
    m = _ELEM_COORDS.search(element)
    if not m:
        return None
    return {"chrom": m.group(1), "start": int(m.group(2)), "end": int(m.group(3)),
            "tag": m.group(4), "clade_suffix": m.group(5)}


def accession_from_path(path: str | Path) -> str:
    """``.../100620.Chr_scaffolds.fasta.cls.tsv`` -> ``100620.Chr_scaffolds``."""
    return _FASTA_SUFFIX.sub("", Path(path).name)


def nsid(accession: str, element: str) -> str:
    return f"{accession}{NS_SEP}{element}"


def split_nsid(namespaced: str) -> tuple[str, str]:
    acc, _, element = namespaced.partition(NS_SEP)
    return acc, element


def element_of(seqid_or_header: str) -> str:
    """Element id = everything before the first ``|`` (strip a leading ``>``)."""
    s = seqid_or_header[1:] if seqid_or_header.startswith(">") else seqid_or_header
    return s.split()[0].split("|", 1)[0]


def domain_of(seqid: str) -> str | None:
    """Pull the domain name out of a ``...:Ty3-GAG`` style id."""
    m = _DOMAIN_FROM_ID.search(seqid.split()[0])
    return m.group(1) if m else None


# --------------------------------------------------------------------------- #
# cls.tsv
# --------------------------------------------------------------------------- #
def parse_cls_tsv(path: str | Path):
    """Yield one dict per element with its per-domain clade calls.

    Keys: element, order, superfamily, clade, complete, strand, domains (dict
    domain->clade in the order TEsorter listed them).
    """
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            domains: dict[str, str] = {}
            for tok in f[6].split():
                dom, _, clade = tok.partition("|")
                if dom:
                    domains[dom] = clade
            yield {
                "element": f[0], "order": f[1], "superfamily": f[2],
                "clade": f[3], "complete": f[4], "strand": f[5],
                "domains": domains,
            }


# --------------------------------------------------------------------------- #
# dom.faa
# --------------------------------------------------------------------------- #
def parse_dom_faa(path: str | Path):
    """Yield (element, domain, attrs: dict, seq) for every domain peptide.

    ``attrs`` includes gene, clade, coverage, evalue, probability when present.
    Domain is taken from ``gene=`` (falling back to the id suffix).
    """
    element = domain = None
    attrs: dict[str, str] = {}
    seq: list[str] = []

    def _emit():
        if element is not None and seq:
            return (element, domain, attrs, "".join(seq))
        return None

    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                rec = _emit()
                if rec:
                    yield rec
                head = line[1:].rstrip("\n")
                first = head.split(None, 1)
                seqid = first[0]
                rest = first[1] if len(first) > 1 else ""
                element = element_of(seqid)
                attrs = {k: v for k, v in _ATTR.findall(rest)}
                domain = attrs.get("gene") or domain_of(seqid)
                seq = []
            else:
                seq.append(line.strip())
        rec = _emit()
        if rec:
            yield rec


# --------------------------------------------------------------------------- #
# dom.tsv
# --------------------------------------------------------------------------- #
def parse_dom_tsv(path: str | Path):
    """Yield dict per domain hit: element, domain, length, evalue, coverage,
    probability, score."""
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            yield {
                "element": element_of(f[0]),
                "domain": domain_of(f[0]),
                "length": int(f[1]) if f[1].isdigit() else None,
                "evalue": _to_float(f[2]),
                "coverage": _to_float(f[3]),
                "probability": _to_float(f[4]),
                "score": _to_float(f[5]),
            }


def _to_float(x: str):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def is_athila(rec: dict) -> bool:
    """Keep LTR / Gypsy / Athila-clade elements only."""
    return (rec["order"] == "LTR" and rec["superfamily"] == "Gypsy"
            and rec["clade"] == "Athila")
