from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class Verdict(str, Enum):
    CHANGED = "CHANGED"
    UNCHANGED = "UNCHANGED"
    UNRESOLVABLE_TECHNICAL = "UNRESOLVABLE_TECHNICAL"
    UNRESOLVABLE_SCIENTIFIC = "UNRESOLVABLE_SCIENTIFIC"


class ClinSigTier(str, Enum):
    BENIGN = "BENIGN"
    VUS_OR_UNKNOWN = "VUS_OR_UNKNOWN"
    PATHOGENIC = "PATHOGENIC"
    OTHER = "OTHER"


@dataclass(frozen=True)
class BiallelicVariant:
    chrom: str
    pos: int  # 1-based VCF POS
    ref: str
    alt: str
    id: str
    info: dict[str, str]


@dataclass(frozen=True)
class NormalizedAllele:
    source_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    vrs_id: Optional[str]
    technical_error: Optional[str] = None


@dataclass(frozen=True)
class SnapshotHit:
    snapshot: str
    vrs_id: str
    clnsig_raw: str
    tier: ClinSigTier
    geneinfo: str | None = None


@dataclass(frozen=True)
class DiffRow:
    allele: NormalizedAllele
    verdict: Verdict
    old: SnapshotHit | None
    new: SnapshotHit | None

