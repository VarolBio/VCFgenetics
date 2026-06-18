from __future__ import annotations

import argparse
import json
import sys

from .bcftools import BcftoolsRequiredError
from .compare import verdict_for
from .normalize import normalize_watchlist_vcf
from .models import DiffRow, Verdict
from .snapshot import load_snapshot
from .store import lookup_vrs, open_store


def _cmd_normalize(args: argparse.Namespace) -> int:
    alleles = normalize_watchlist_vcf(args.watchlist, reference_fasta=args.reference)
    if args.emit == "json":
        rows = [
            {
                "id": a.source_id,
                "chrom": a.chrom,
                "pos": a.pos,
                "ref": a.ref,
                "alt": a.alt,
                "vrs_id": a.vrs_id,
                "technical_error": a.technical_error,
            }
            for a in alleles
        ]
        print(json.dumps(rows, sort_keys=True, separators=(",", ":")))
    else:
        for a in alleles:
            print(f"{a.source_id}\t{a.chrom}:{a.pos}\t{a.ref}>{a.alt}\t{a.vrs_id or 'NO_VRS'}")
    return 0


def _cmd_ingest(args: argparse.Namespace) -> int:
    alleles = normalize_watchlist_vcf(args.watchlist, reference_fasta=args.reference)
    if args.emit == "json":
        rows = [
            {
                "id": a.source_id,
                "chrom": a.chrom,
                "pos": a.pos,
                "ref": a.ref,
                "alt": a.alt,
                "vrs_id": a.vrs_id,
                "technical_error": a.technical_error,
            }
            for a in alleles
        ]
        print(json.dumps(rows, sort_keys=True, separators=(",", ":")))
    else:
        for a in alleles:
            print(f"{a.source_id}\t{a.chrom}:{a.pos}\t{a.ref}>{a.alt}\t{a.vrs_id or 'NO_VRS'}")
    return 0


def _cmd_lookup(args: argparse.Namespace) -> int:
    # --variant CHR:POS:REF:ALT
    try:
        chrom, pos, ref, alt = args.variant.split(":", 3)
        pos_i = int(pos)
    except Exception:
        print("ERROR: --variant must be CHR:POS:REF:ALT", file=sys.stderr)
        return 2

    from .vrs import vrs_allele_id

    vrs_id = vrs_allele_id(
        chrom=chrom,
        pos=pos_i,
        ref=ref,
        alt=alt,
        reference_fasta=args.reference,
    )
    if args.db:
        from .snapshot import _snapshot_date_from_path

        conn = open_store(args.db)
        try:
            hit = lookup_vrs(
                conn,
                snapshot_date=_snapshot_date_from_path(args.snapshot),
                vrs_id=vrs_id,
            )
        finally:
            conn.close()
        if hit is None:
            snap = load_snapshot(
                args.snapshot, reference_fasta=args.reference, db_path=args.db
            )
            hit = snap.get(vrs_id)
    else:
        snap = load_snapshot(args.snapshot, reference_fasta=args.reference)
        hit = snap.get(vrs_id)
    if hit is None:
        print("NO_MATCH")
        return 0
    print(f"MATCH\t{hit.vrs_id}\t{hit.tier}\t{hit.clnsig_raw}\t{hit.geneinfo or ''}")
    return 0


def _cmd_diff(args: argparse.Namespace) -> int:
    alleles = normalize_watchlist_vcf(args.watchlist, reference_fasta=args.reference)
    old = load_snapshot(args.old, reference_fasta=args.reference, db_path=args.db)
    new = load_snapshot(args.new, reference_fasta=args.reference, db_path=args.db)

    rows: list[DiffRow] = []
    counts = {v.value: 0 for v in Verdict}
    changed_rows = []
    tech_rows = []
    sci_rows = []

    for a in alleles:
        if a.technical_error or not a.vrs_id:
            verdict = Verdict.UNRESOLVABLE_TECHNICAL
            row = DiffRow(allele=a, verdict=verdict, old=None, new=None)
            rows.append(row)
            counts[verdict.value] += 1
            tech_rows.append(row)
            continue

        old_hit = old.get(a.vrs_id)
        new_hit = new.get(a.vrs_id)
        verdict = verdict_for(
            old_tier=(old_hit.tier if old_hit else None),
            new_tier=(new_hit.tier if new_hit else None),
        )
        row = DiffRow(allele=a, verdict=verdict, old=old_hit, new=new_hit)
        rows.append(row)
        counts[verdict.value] += 1
        if verdict is Verdict.CHANGED:
            changed_rows.append(row)
        elif verdict is Verdict.UNRESOLVABLE_SCIENTIFIC:
            sci_rows.append(row)

    total = len(alleles)
    ok = sum(counts.values()) == total

    if args.emit == "json":
        out_rows = []
        for r in rows:
            out_rows.append(
                {
                    "id": r.allele.source_id,
                    "vrs_id": r.allele.vrs_id,
                    "verdict": r.verdict.value,
                    "old_tier": (r.old.tier.value if r.old else None),
                    "new_tier": (r.new.tier.value if r.new else None),
                }
            )
        print(json.dumps(out_rows, sort_keys=True, separators=(",", ":")))
    else:
        print("COUNTS", json.dumps(counts, sort_keys=True))
        print("TOTAL_INGESTED", total)
        print("CONSERVATION_OK" if ok else "CONSERVATION_BROKEN")
        if changed_rows:
            print("CHANGED")
            for r in changed_rows:
                print(
                    f"- {r.allele.source_id}\t{r.allele.vrs_id}\t{(r.old.tier.value if r.old else 'NA')}->{(r.new.tier.value if r.new else 'NA')}"
                )
        if tech_rows:
            print("UNRESOLVABLE_TECHNICAL")
            for r in tech_rows:
                print(f"- {r.allele.source_id}\t{r.allele.technical_error}")
        if sci_rows:
            print("UNRESOLVABLE_SCIENTIFIC")
            for r in sci_rows:
                gi = r.new.geneinfo if r.new else (r.old.geneinfo if r.old else None)
                print(f"- {r.allele.source_id}\t{r.allele.vrs_id}\t{gi or ''}")

        # Positive-control friendliness: echo any gene symbols seen.
        genes = set()
        for hit in list(old.values()) + list(new.values()):
            if hit.geneinfo:
                genes.add(hit.geneinfo.split(":")[0])
        if genes:
            print("GENES", ",".join(sorted(genes)))

    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="vcfgenetics-monitor")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_norm = sub.add_parser("normalize")
    p_norm.add_argument("--watchlist", required=True)
    p_norm.add_argument("--reference", required=True)
    p_norm.add_argument("--emit", choices=["json", "text"], default="text")
    p_norm.set_defaults(func=_cmd_normalize)

    # Placeholders for other commands required by the spec suite.
    p_ing = sub.add_parser("ingest")
    p_ing.add_argument("--watchlist", required=True)
    p_ing.add_argument("--reference", required=True)
    p_ing.add_argument("--emit", choices=["json", "text"], default="text")
    p_ing.set_defaults(func=_cmd_ingest)

    p_diff = sub.add_parser("diff")
    p_diff.add_argument("--watchlist", required=True)
    p_diff.add_argument("--old", required=True)
    p_diff.add_argument("--new", required=True)
    p_diff.add_argument("--reference", required=True)
    p_diff.add_argument("--db", help="Optional SQLite path to persist snapshot rows")
    p_diff.add_argument("--emit", choices=["json", "text"], default="text")
    p_diff.set_defaults(func=_cmd_diff)

    p_lookup = sub.add_parser("lookup")
    p_lookup.add_argument("--variant", required=True)
    p_lookup.add_argument("--snapshot", required=True)
    p_lookup.add_argument("--reference", required=True)
    p_lookup.add_argument("--db", help="Optional SQLite path to persist snapshot rows")
    p_lookup.set_defaults(func=_cmd_lookup)

    ns = p.parse_args(argv)
    try:
        return int(ns.func(ns))
    except BcftoolsRequiredError as e:
        print(str(e), file=sys.stderr)
        return 2
    except NotImplementedError as e:
        print(f"NOT_IMPLEMENTED: {e}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())

