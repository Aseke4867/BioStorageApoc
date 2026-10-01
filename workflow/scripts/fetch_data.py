"""Download one input (or a set of S3 objects) and write a provenance record next to it.

  fetch_data.py url  --url URL --out FILE [--md5 HEX] --key NAME --accession ACC --note TXT
  fetch_data.py s3cat --bucket URL --prefix P --out FILE --key NAME   (concatenate every object
                                                                    under P, e.g. 96 FASTA files)
  fetch_data.py merge --prov A.prov.tsv B.prov.tsv ... --out downloads.tsv

A URL that starts with local: is copied from disk instead (used by the test dataset).
Each download streams to FILE.part, is checked against the archive's MD5 when one is known,
and only then renamed, so an interrupted or corrupt download never looks finished.
The .prov.tsv sidecar holds: key, accession, url, bytes, md5, sha256, retrieved_utc, note.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

UA = {"User-Agent": "aitu-bioinfo-course-project/1.0"}
PROV_COLS = ["key", "accession", "url", "bytes", "md5", "sha256", "retrieved_utc", "note"]


def hashes(path: Path) -> tuple[str, str]:
    m, s = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            m.update(chunk)
            s.update(chunk)
    return m.hexdigest(), s.hexdigest()


def stream(url: str, dest: Path, retries: int = 5) -> None:
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120, headers=UA) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            return
        except requests.RequestException as exc:
            print(f"attempt {attempt}/{retries} failed: {exc}", file=sys.stderr)
            if attempt == retries:
                raise
            time.sleep(2 ** attempt)


def s3_keys(bucket: str, prefix: str) -> list[str]:
    keys, token = [], None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        xml = requests.get(bucket + "/", params=params, timeout=60, headers=UA).text
        keys += re.findall(r"<Key>([^<]+)</Key>", xml)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            return keys
        token = m.group(1)


def write_prov(out: Path, rec: dict) -> None:
    with open(str(out) + ".prov.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PROV_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(rec)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["url", "s3cat", "merge"])
    ap.add_argument("--url")
    ap.add_argument("--bucket")
    ap.add_argument("--prefix")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--md5", default="")
    ap.add_argument("--key", default="")
    ap.add_argument("--accession", default="")
    ap.add_argument("--note", default="")
    ap.add_argument("--prov", type=Path, nargs="*")
    ap.add_argument("--reuse", help="file already downloaded from --url; verified, not re-downloaded")
    a = ap.parse_args(argv)

    if a.mode == "merge":
        rows = []
        for p in a.prov:
            with open(p, newline="") as fh:
                rows += list(csv.DictReader(fh, delimiter="\t"))
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=PROV_COLS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        return 0

    a.out.parent.mkdir(parents=True, exist_ok=True)
    part = a.out.with_name(a.out.name + ".part")
    if a.mode == "url" and a.reuse and Path(a.reuse).exists():
        # an earlier manual download of the same URL: verify it instead of downloading again
        shutil.move(a.reuse, part)
        src = a.url
    elif a.mode == "url":
        if a.url.startswith("local:"):
            src_path = a.url[len("local:"):]
            if src_path.endswith(".gz") and not str(a.out).endswith(".gz"):
                with gzip.open(src_path, "rb") as fi, open(part, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copyfile(src_path, part)
        else:
            stream(a.url, part)
        src = a.url
    else:
        keys = s3_keys(a.bucket, a.prefix)
        if not keys:
            raise SystemExit(f"no objects under {a.bucket}/{a.prefix}")
        with open(part, "wb") as fo:
            for k in keys:
                r = requests.get(f"{a.bucket}/{k}", timeout=60, headers=UA)
                r.raise_for_status()
                fo.write(r.content if r.content.endswith(b"\n") else r.content + b"\n")
        src = f"{a.bucket}/{a.prefix} ({len(keys)} objects)"
    md5, sha = hashes(part)
    if a.md5 and md5 != a.md5:
        part.unlink()
        raise SystemExit(f"MD5 mismatch for {a.out}: got {md5}, archive says {a.md5}")
    part.replace(a.out)
    write_prov(a.out, {"key": a.key, "accession": a.accession, "url": src, "bytes": a.out.stat().st_size,
                       "md5": md5, "sha256": sha,
                       "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                       "note": a.note + ("" if not a.md5 else "; md5 verified against archive")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
