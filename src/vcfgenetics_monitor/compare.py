from __future__ import annotations

from dataclasses import dataclass

from .models import ClinSigTier, Verdict


def is_confident_tier(t: ClinSigTier) -> bool:
    return t in {ClinSigTier.BENIGN, ClinSigTier.PATHOGENIC}


def verdict_for(*, old_tier: ClinSigTier | None, new_tier: ClinSigTier | None) -> Verdict:
    """
    Four-verdict model.

    - Scientific-unresolvable: absent in new snapshot, or present only as VUS/unknown/other.
    - Changed: meaningful tier movement across benign/VUS/pathogenic boundaries.
    - Unchanged: confident tier stable.
    """
    if new_tier is None:
        return Verdict.UNRESOLVABLE_SCIENTIFIC
    if not is_confident_tier(new_tier):
        return Verdict.UNRESOLVABLE_SCIENTIFIC

    if old_tier is None:
        # Newly classified into a confident tier; for MVP treat as unchanged baseline.
        return Verdict.UNCHANGED

    if old_tier != new_tier:
        return Verdict.CHANGED
    return Verdict.UNCHANGED

