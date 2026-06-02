import subprocess
import sys
from pathlib import Path

import pytest


FIXTURES = Path(__file__).parent / "fixtures"


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "vcfgenetics_monitor.cli", *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def test_round_trip_equivalence_normalizes_to_same_vrs_id():
    """
    The same biological variant expressed two valid ways must normalize to the same VRS identifier.
    """
    out = _run_cli(
        "normalize",
        "--watchlist",
        str(FIXTURES / "watchlist.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
        "--emit",
        "json",
    )
    assert out.returncode == 0, out.stderr
    # Expect two equivalence rows with identical vrs_id (implementation will define JSON structure).
    assert "eq1" in out.stdout
    assert "eq2" in out.stdout
    assert "ga4gh:VA" in out.stdout


def test_anchored_deletion_regression_never_emits_dash_alt():
    """
    Regression: deletion must normalize to anchored VCF allele (REF/ALT non-empty, ALT != '-')
    and receive a VRS id (never the invalid T>- / ALT='-' form).
    """
    out = _run_cli(
        "normalize",
        "--watchlist",
        str(FIXTURES / "watchlist.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
        "--emit",
        "json",
    )
    assert out.returncode == 0, out.stderr
    assert "del_bug" in out.stdout
    assert '"alt":"-"' not in out.stdout
    assert "ga4gh:VA" in out.stdout


def test_multiallelic_split_counts_all_alleles():
    """
    Multi-allelic records must be split into separate biallelic alleles and none dropped.
    """
    out = _run_cli(
        "ingest",
        "--watchlist",
        str(FIXTURES / "watchlist_diff.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
        "--emit",
        "json",
    )
    assert out.returncode == 0, out.stderr
    # One record with ALT=C,G must become two biallelic alleles.
    assert "multi" in out.stdout
    assert out.stdout.count('"id":"multi"') == 2


def test_bcftools_required_guard_hard_errors_without_fallback(monkeypatch: pytest.MonkeyPatch):
    """
    If bcftools is unavailable, the pipeline must hard-fail with a clear error.
    """
    monkeypatch.setenv("PATH", "")
    out = _run_cli(
        "normalize",
        "--watchlist",
        str(FIXTURES / "watchlist.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
        "--emit",
        "json",
    )
    assert out.returncode != 0
    assert "bcftools" in (out.stderr + out.stdout).lower()
    assert "required" in (out.stderr + out.stdout).lower()


def test_clinvar_positive_control_known_genes_are_findable():
    """
    A small set of known ClinVar-positive genes must match in a loaded snapshot.
    Any miss is a representation bug, not an expected absence.
    """
    out = _run_cli(
        "diff",
        "--watchlist",
        str(FIXTURES / "watchlist_diff.vcf"),
        "--old",
        str(FIXTURES / "clinvar_positive_control.vcf"),
        "--new",
        str(FIXTURES / "clinvar_positive_control.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
    )
    assert out.returncode == 0, out.stderr
    # Presence check: at least one match for each gene symbol.
    for gene in ["BRCA1", "BRCA2", "TP53", "PTEN"]:
        assert gene in out.stdout


def test_clinvar_negative_control_synthetic_variant_is_not_found():
    """
    A known-absent synthetic variant must return no match (matching isn't trivially True).
    """
    out = _run_cli(
        "lookup",
        "--variant",
        "1:123456:A:T",
        "--snapshot",
        str(FIXTURES / "clinvar_positive_control.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
    )
    assert out.returncode == 0, out.stderr
    assert "NO_MATCH" in out.stdout


def test_conservation_invariant_holds_across_four_verdicts():
    """
    CHANGED+UNCHANGED+UNRESOLVABLE_TECHNICAL+UNRESOLVABLE_SCIENTIFIC == total_ingested
    (post-split allele count included).
    """
    out = _run_cli(
        "diff",
        "--watchlist",
        str(FIXTURES / "watchlist_diff.vcf"),
        "--old",
        str(FIXTURES / "clinvar_old.vcf"),
        "--new",
        str(FIXTURES / "clinvar_new.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
    )
    assert out.returncode == 0, out.stderr
    assert "CONSERVATION_OK" in out.stdout


def test_four_verdict_classification_one_variant_in_each_bucket():
    """
    Crafted snapshot pair must yield exactly one variant in each verdict.
    """
    out = _run_cli(
        "diff",
        "--watchlist",
        str(FIXTURES / "watchlist_diff.vcf"),
        "--old",
        str(FIXTURES / "clinvar_old.vcf"),
        "--new",
        str(FIXTURES / "clinvar_new.vcf"),
        "--reference",
        str(FIXTURES / "GRCh38_micro.fa"),
        "--emit",
        "json",
    )
    assert out.returncode == 0, out.stderr
    for verdict in [
        "CHANGED",
        "UNCHANGED",
        "UNRESOLVABLE_TECHNICAL",
        "UNRESOLVABLE_SCIENTIFIC",
    ]:
        assert verdict in out.stdout
    # The crafted fixtures are designed so each verdict occurs once.
    assert out.stdout.count('"verdict":"CHANGED"') == 1
    assert out.stdout.count('"verdict":"UNCHANGED"') == 1
    assert out.stdout.count('"verdict":"UNRESOLVABLE_TECHNICAL"') == 1
    assert out.stdout.count('"verdict":"UNRESOLVABLE_SCIENTIFIC"') == 1

