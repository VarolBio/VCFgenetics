from __future__ import annotations

import argparse
from pathlib import Path

import requests


BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--old", required=True, help="ClinVar snapshot date YYYYMMDD")
    p.add_argument("--new", required=True, help="ClinVar snapshot date YYYYMMDD")
    p.add_argument("--out", default="data/clinvar", help="Output directory")
    ns = p.parse_args(argv)

    out = Path(ns.out)
    for d in [ns.old, ns.new]:
        for suffix in [".vcf.gz", ".vcf.gz.tbi"]:
            fn = f"clinvar_{d}{suffix}"
            url = f"{BASE}/{fn}"
            dest = out / fn
            if dest.exists() and dest.stat().st_size > 0:
                continue
            print(f"Downloading {url} -> {dest}")
            _download(url, dest)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

