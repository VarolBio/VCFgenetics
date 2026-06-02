## Status: proven vs assumed (MVP)

### Proven (by `pytest`)
- **bcftools requirement is enforced**: absence hard-fails; no silent fallback.
- **Anchored deletion regression is blocked**: no invalid dash ALT emitted.
- **Multi-allelic records are split and counted**: no silent drops.
- **VRS canonicalization works**: equivalent representations normalize to the same `ga4gh:VA.*` id.
- **Four-verdict output is exhaustive**: conservation invariant holds; each allele lands in exactly one verdict.
- **Positive/negative controls**: known-present fixture variants are findable; synthetic absent is not.

### Assumed (not yet proven end-to-end on full ClinVar)
- **Reference FASTA compatibility**: user-supplied GRCh38 must match ClinVar contig naming and sequence content.
- **ClinVar semantics**: MVP tiering uses `CLNSIG` only; review-status nuance is not modeled.
- **Scaling**: full weekly ClinVar loads and diff performance characteristics are not optimized yet.

