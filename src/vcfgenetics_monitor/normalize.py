from __future__ import annotations

from pathlib import Path

from .bcftools import BcftoolsRequiredError, bcftools_norm_split_leftalign_vcf_text
from .ingest import read_vcf_records, to_vcf_text
from .models import NormalizedAllele
from .vrs import vrs_allele_id


def _parse_norm_vcf(vcf_text: str) -> list[tuple[str, int, str, str, str]]:
    """
    Return list of (chrom, pos, id, ref, alt) from bcftools-normalized VCF.
    """
    out: list[tuple[str, int, str, str, str]] = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        chrom, pos, vid, ref, alt = line.split("\t", 5)[:5]
        out.append((chrom, int(pos), vid, ref, alt))
    return out


def normalize_watchlist_vcf(
    watchlist_vcf: str | Path, *, reference_fasta: str | Path
) -> list[NormalizedAllele]:
    """
    bcftools-based split + left-align. Any normalization failure becomes a technical error entry.
    """
    records = read_vcf_records(watchlist_vcf)
    alleles: list[NormalizedAllele] = []
    for r in records:
        try:
            raw_vcf = to_vcf_text([r])
            norm_vcf = bcftools_norm_split_leftalign_vcf_text(
                raw_vcf, reference_fasta=str(reference_fasta)
            )
            for chrom, pos, vid, ref, alt in _parse_norm_vcf(norm_vcf):
                alleles.append(
                    NormalizedAllele(
                        source_id=vid,
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        vrs_id=vrs_allele_id(
                            chrom=chrom,
                            pos=pos,
                            ref=ref,
                            alt=alt,
                            reference_fasta=str(reference_fasta),
                        ),
                        technical_error=None,
                    )
                )
        except BcftoolsRequiredError:
            raise
        except Exception as e:
            alleles.append(
                NormalizedAllele(
                    source_id=r.id,
                    chrom=r.chrom,
                    pos=r.pos,
                    ref=r.ref,
                    alt=r.alt,
                    vrs_id=None,
                    technical_error=str(e),
                )
            )
    return alleles

