from __future__ import annotations

from .models import ClinSigTier


def parse_clnsig_raw(info: dict[str, str]) -> str:
    return info.get("CLNSIG", "")


def parse_geneinfo(info: dict[str, str]) -> str | None:
    return info.get("GENEINFO")


def clnsig_to_tier(clnsig_raw: str) -> ClinSigTier:
    # ClinVar can be multi-valued; for MVP we treat presence of a high-confidence
    # term as enough to assign the tier.
    s = clnsig_raw.replace(" ", "_")
    parts = [p for p in s.split(",") if p]
    if any(p in {"Pathogenic", "Likely_pathogenic"} for p in parts):
        return ClinSigTier.PATHOGENIC
    if any(p in {"Benign", "Likely_benign"} for p in parts):
        return ClinSigTier.BENIGN
    if any(p in {"Uncertain_significance"} for p in parts):
        return ClinSigTier.VUS_OR_UNKNOWN
    if not clnsig_raw:
        return ClinSigTier.VUS_OR_UNKNOWN
    return ClinSigTier.OTHER

