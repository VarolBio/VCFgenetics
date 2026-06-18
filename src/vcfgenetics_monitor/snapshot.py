from __future__ import annotations

from pathlib import Path

from .bcftools import bcftools_norm_split_leftalign_vcf_text
from .clinvar import clnsig_to_tier, parse_clnsig_raw, parse_geneinfo
from .ingest import read_vcf_records, to_vcf_text
from .models import SnapshotHit
from .normalize import _parse_norm_vcf
from .store import load_snapshot_hits, open_store, query_snapshot
from .vrs import vrs_allele_id


def _snapshot_date_from_path(path: str | Path) -> str:
    p = Path(path)
    stem = p.name
    if stem.endswith(".vcf.gz"):
        stem = stem[: -len(".vcf.gz")]
    elif stem.endswith(".vcf"):
        stem = stem[: -len(".vcf")]
    if stem.startswith("clinvar_"):
        return stem[len("clinvar_") :]
    return stem


def build_snapshot_hits(snapshot_path: str | Path, *, reference_fasta: str) -> dict[str, SnapshotHit]:
    rows = read_vcf_records(snapshot_path)
    hits: dict[str, SnapshotHit] = {}
    snap_label = str(snapshot_path)

    for r in rows:
        raw_vcf = to_vcf_text([r])
        norm_text = bcftools_norm_split_leftalign_vcf_text(
            raw_vcf, reference_fasta=reference_fasta
        )
        for chrom, pos, _vid, ref, alt in _parse_norm_vcf(norm_text):
            vrs_id = vrs_allele_id(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                reference_fasta=reference_fasta,
            )
            clnsig_raw = parse_clnsig_raw(r.info)
            tier = clnsig_to_tier(clnsig_raw)
            geneinfo = parse_geneinfo(r.info)
            hits.setdefault(
                vrs_id,
                SnapshotHit(
                    snapshot=snap_label,
                    vrs_id=vrs_id,
                    clnsig_raw=clnsig_raw,
                    tier=tier,
                    geneinfo=geneinfo,
                ),
            )
    return hits


def load_snapshot(
    snapshot_path: str | Path,
    *,
    reference_fasta: str,
    db_path: str | Path | None = None,
) -> dict[str, SnapshotHit]:
    hits = build_snapshot_hits(snapshot_path, reference_fasta=reference_fasta)
    if db_path is None:
        return hits

    snapshot_date = _snapshot_date_from_path(snapshot_path)
    conn = open_store(db_path)
    try:
        load_snapshot_hits(conn, snapshot_date=snapshot_date, hits=hits)
        return query_snapshot(conn, snapshot_date=snapshot_date)
    finally:
        conn.close()
