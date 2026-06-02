from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from ga4gh.core.identifiers import ga4gh_identify
from ga4gh.vrs.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
    SequenceReference,
)
from ga4gh.vrs.models import sha512t24u


def _read_fasta_contig_sequence(reference_fasta: str | Path, contig: str) -> str:
    """
    Minimal FASTA reader sufficient for the committed micro-reference fixture.
    """
    p = Path(reference_fasta)
    want = f">{contig}"
    seq_lines: list[str] = []
    in_contig = False
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if in_contig:
                break
            in_contig = line.strip() == want
            continue
        if in_contig:
            seq_lines.append(line.strip())
    if not seq_lines:
        raise ValueError(f"Contig {contig!r} not found in reference FASTA: {p}")
    return "".join(seq_lines).upper()


@lru_cache(maxsize=256)
def contig_refget_accession(reference_fasta: str, contig: str) -> str:
    """
    Return the refget accession (sha512t24u) for the contig sequence.
    For MVP tests we compute from the micro-reference FASTA.
    """
    seq = _read_fasta_contig_sequence(reference_fasta, contig)
    return "SQ." + sha512t24u(seq.encode("ascii"))


def vrs_allele_id(*, chrom: str, pos: int, ref: str, alt: str, reference_fasta: str) -> str:
    """
    Compute GA4GH VRS Allele identifier for a normalized VCF biallelic allele.

    VRS uses 0-based inter-residue coordinates:
      start = POS-1
      end   = POS-1 + len(REF)
    """
    if pos < 1:
        raise ValueError("VCF POS must be 1-based and >= 1")
    if not ref or not alt or alt == "-":
        raise ValueError("Invalid VCF allele (empty REF/ALT or ALT='-')")

    refget = contig_refget_accession(str(reference_fasta), chrom)
    loc = SequenceLocation(
        sequenceReference=SequenceReference(refgetAccession=refget),
        start=pos - 1,
        end=pos - 1 + len(ref),
    )
    allele = Allele(location=loc, state=LiteralSequenceExpression(sequence=alt))
    curie = ga4gh_identify(allele, in_place="never")
    if curie is None:
        raise ValueError("Unable to compute GA4GH identifier for allele")
    return curie

