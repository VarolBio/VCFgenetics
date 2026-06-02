# ClinVar snapshot diff MVP (CLI-only)

## The product (one sentence)
**A variant goes in → it is canonically normalized (real `bcftools norm`, GA4GH VRS identifier) → baseline classification is stored → two dated ClinVar weekly snapshots are diffed → each watched variant emits exactly one of four verdicts, and CHANGED variants are printed as alert candidates.**

## Architecture (snapshot diff, not live polling)
- **ClinVar input**: bulk **weekly snapshot VCF** (`vcf_GRCh38/clinvar_YYYYMMDD.vcf.gz`). No per-variant API calls.
- **Detection mechanism**: diff **older vs newer** snapshot using the same code path used for MVP tests.

## Four-verdict model (never drop variants)
Every ingested (post-split) allele must land in exactly one verdict:
- **CHANGED**: meaningful tier change between snapshots (e.g. VUS → Pathogenic).
- **UNCHANGED**: stable, confidently-known tier.
- **UNRESOLVABLE_TECHNICAL**: our failure (normalization / parse / match bug). Must be loud.
- **UNRESOLVABLE_SCIENTIFIC**: pipeline worked but ClinVar has no confident answer yet (absent or only VUS/unknown).

Conservation invariant must always hold:
\[
CHANGED + UNCHANGED + UNRESOLVABLE\_TECHNICAL + UNRESOLVABLE\_SCIENTIFIC = total\_ingested
\]

## Hard requirements
- **`bcftools` is required**. If absent, the pipeline aborts (no Python fallback).
- **GRCh38 consistency**: reference FASTA contig naming must match the VCFs (ClinVar GRCh38 VCFs use numeric contigs like `13`, `17`, etc.).

## Quickstart (dev)
Create a venv and install:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e ".[dev]"
pytest
```

## Prove on real ClinVar snapshots (dated weekly releases)
1) Download two dated ClinVar GRCh38 snapshot VCFs:

```bash
source .venv/bin/activate
python scripts/fetch_clinvar_snapshots.py --old 20260517 --new 202606?? --out data/clinvar
```

2) Run a diff against your watchlist VCF (VCF input is your registry/baseline side):

```bash
source .venv/bin/activate
python -m vcfgenetics_monitor.cli diff \
  --watchlist path/to/watchlist.vcf \
  --old data/clinvar/clinvar_20260517.vcf.gz \
  --new data/clinvar/clinvar_YYYYMMDD.vcf.gz \
  --reference /path/to/GRCh38.fa.gz
```

## What this proves / what it does not prove
- **Proves**: bcftools-based normalization robustness (including anchored deletions + multi-allelic splitting), deterministic VRS canonical keys, and snapshot-diff reclassification detection with a conservation invariant.
- **Does not prove**: clinical validity, completeness of ClinVar capture, ACMG/AMP interpretation, or that this is clinical-grade software. This is not a clinical system.

