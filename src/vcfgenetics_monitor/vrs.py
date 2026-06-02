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


_GRCH38_CONTIG_TO_SQ: dict[str, str] = {
    # Source: widely-used GRCh38 contig -> refget (SQ.*) mapping (e.g. gnomAD tooling).
    # This avoids hashing entire chromosome sequences at runtime.
    "chr1": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
    "chr2": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
    "chr3": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
    "chr4": "SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc",
    "chr5": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
    "chr6": "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
    "chr7": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
    "chr8": "SQ.209Z7zJ-mFypBEWLk4rNC6S_OxY5p7bs",
    "chr9": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI",
    "chr10": "SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB",
    "chr11": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
    "chr12": "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
    "chr13": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
    "chr14": "SQ.eK4D2MosgK_ivBkgi6FVPg5UXs1bYESm",
    "chr15": "SQ.AsXvWL1-2i5U_buw6_niVIxD6zTbAuS6",
    "chr16": "SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
    "chr17": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
    "chr18": "SQ.vWwFhJ5lQDMhh-czg06YtlWqu0lvFAZV",
    "chr19": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
    "chr20": "SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
    "chr21": "SQ.5ZUqxCmDDgN4xTRbaSjN8LwgZironmB8",
    "chr22": "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
    "chrX": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
    "chrY": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
}


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
    p = Path(reference_fasta)
    # For tiny committed micro-references, hashing is fast and ensures tests are self-contained.
    if p.exists() and p.stat().st_size < 5_000_000:
        seq = _read_fasta_contig_sequence(reference_fasta, contig)
        return "SQ." + sha512t24u(seq.encode("ascii"))

    key = contig if contig.startswith("chr") else f"chr{contig}"
    if key in _GRCH38_CONTIG_TO_SQ:
        return _GRCH38_CONTIG_TO_SQ[key]

    raise ValueError(
        f"Unknown GRCh38 contig refget accession for contig={contig!r}. "
        "Provide a small reference for tests, or use a primary GRCh38 contig."
    )


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

