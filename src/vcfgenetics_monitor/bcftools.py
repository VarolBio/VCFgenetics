from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Iterable, Sequence


class BcftoolsRequiredError(RuntimeError):
    pass


@dataclass(frozen=True)
class BcftoolsResult:
    returncode: int
    stdout: str
    stderr: str


def require_bcftools() -> None:
    if shutil.which("bcftools") is None:
        raise BcftoolsRequiredError(
            "bcftools is required for normalization (no Python fallback). "
            "Install bcftools and ensure it is on PATH."
        )


def run_bcftools(args: Sequence[str], *, input_text: str | None = None) -> BcftoolsResult:
    require_bcftools()
    proc = subprocess.run(
        ["bcftools", *args],
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return BcftoolsResult(proc.returncode, proc.stdout, proc.stderr)


def bcftools_norm_split_leftalign_vcf_text(vcf_text: str, *, reference_fasta: str) -> str:
    """
    Run: bcftools norm -f <ref> -m-  (split multiallelic + left-align).
    Input/Output is VCF text.
    """
    res = run_bcftools(
        ["norm", "-f", reference_fasta, "-m-", "-"],
        input_text=vcf_text,
    )
    if res.returncode != 0:
        raise RuntimeError(f"bcftools norm failed: {res.stderr.strip()}")
    return res.stdout

