from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .models import BiallelicVariant


def _parse_info(info_str: str) -> dict[str, str]:
    if info_str == ".":
        return {}
    out: dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "true"
    return out


def read_vcf_records(path: str | Path) -> list[BiallelicVariant]:
    """
    Minimal VCF parser for the MVP fixtures: emits one BiallelicVariant per VCF record line.
    Multi-allelic ALT values remain comma-separated here; splitting is delegated to bcftools.
    """
    p = Path(path)
    rows: list[BiallelicVariant] = []
    for line in p.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            raise ValueError(f"Malformed VCF line (expected >=8 fields): {line}")
        chrom, pos, vid, ref, alt, qual, flt, info = fields[:8]
        rows.append(
            BiallelicVariant(
                chrom=chrom,
                pos=int(pos),
                ref=ref,
                alt=alt,
                id=vid,
                info=_parse_info(info),
            )
        )
    return rows


def to_vcf_text(records: Iterable[BiallelicVariant]) -> str:
    recs = list(records)
    contigs = sorted({r.chrom for r in recs})
    info_keys = sorted({k for r in recs for k in r.info.keys()})
    header_lines = ["##fileformat=VCFv4.2"]
    for c in contigs:
        header_lines.append(f"##contig=<ID={c}>")
    for k in info_keys:
        header_lines.append(f"##INFO=<ID={k},Number=1,Type=String,Description=\"{k}\">")
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    lines = ["\n".join(header_lines)]
    for r in recs:
        info = "." if not r.info else ";".join(f"{k}={v}" for k, v in r.info.items())
        lines.append(f"{r.chrom}\t{r.pos}\t{r.id}\t{r.ref}\t{r.alt}\t.\t.\t{info}")
    return "\n".join(lines) + "\n"

