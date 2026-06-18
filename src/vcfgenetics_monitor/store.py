from __future__ import annotations

import sqlite3
from pathlib import Path

from .models import ClinSigTier, SnapshotHit


_SCHEMA = """
CREATE TABLE IF NOT EXISTS clinvar (
    snapshot_date TEXT NOT NULL,
    vrs_id        TEXT NOT NULL,
    clnsig_raw    TEXT NOT NULL,
    clnsig_tier   TEXT NOT NULL,
    geneinfo      TEXT,
    PRIMARY KEY (snapshot_date, vrs_id)
);
CREATE INDEX IF NOT EXISTS idx_clinvar_vrs
    ON clinvar (snapshot_date, vrs_id);
"""


def open_store(path: str | Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(path))
    conn.executescript(_SCHEMA)
    return conn


def load_snapshot_hits(
    conn: sqlite3.Connection,
    *,
    snapshot_date: str,
    hits: dict[str, SnapshotHit],
) -> None:
    conn.execute("DELETE FROM clinvar WHERE snapshot_date = ?", (snapshot_date,))
    conn.executemany(
        """
        INSERT OR REPLACE INTO clinvar
            (snapshot_date, vrs_id, clnsig_raw, clnsig_tier, geneinfo)
        VALUES (?, ?, ?, ?, ?)
        """,
        [
            (
                snapshot_date,
                h.vrs_id,
                h.clnsig_raw,
                h.tier.value,
                h.geneinfo,
            )
            for h in hits.values()
        ],
    )
    conn.commit()


def query_snapshot(conn: sqlite3.Connection, *, snapshot_date: str) -> dict[str, SnapshotHit]:
    rows = conn.execute(
        """
        SELECT vrs_id, clnsig_raw, clnsig_tier, geneinfo
        FROM clinvar
        WHERE snapshot_date = ?
        """,
        (snapshot_date,),
    ).fetchall()
    out: dict[str, SnapshotHit] = {}
    for vrs_id, clnsig_raw, tier_raw, geneinfo in rows:
        out[vrs_id] = SnapshotHit(
            snapshot=snapshot_date,
            vrs_id=vrs_id,
            clnsig_raw=clnsig_raw,
            tier=ClinSigTier(tier_raw),
            geneinfo=geneinfo,
        )
    return out


def lookup_vrs(
    conn: sqlite3.Connection, *, snapshot_date: str, vrs_id: str
) -> SnapshotHit | None:
    row = conn.execute(
        """
        SELECT vrs_id, clnsig_raw, clnsig_tier, geneinfo
        FROM clinvar
        WHERE snapshot_date = ? AND vrs_id = ?
        """,
        (snapshot_date, vrs_id),
    ).fetchone()
    if row is None:
        return None
    vrs_id, clnsig_raw, tier_raw, geneinfo = row
    return SnapshotHit(
        snapshot=snapshot_date,
        vrs_id=vrs_id,
        clnsig_raw=clnsig_raw,
        tier=ClinSigTier(tier_raw),
        geneinfo=geneinfo,
    )
