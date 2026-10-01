"""Quality-score binning for FASTQ (task 4), plus the two measurements that justify it.

Binning replaces the per-base Phred score with a few representative values. Qualities are the
least compressible part of a FASTQ, so fewer distinct symbols means a much smaller archive; the
cost is that variant callers lose information they use to weigh evidence. This script only
rewrites qualities; downstream cost is measured by aligning + calling on each binned file.

Schemes (Phred value ranges -> representative value):
  orig      unchanged
  illumina8 Illumina's 8-level scheme (white paper "Reducing whole-genome data storage
            footprint", 2014): 0-1 unchanged, 2-9->6, 10-19->15, 20-24->22, 25-29->27,
            30-34->33, 35-39->37, >=40->40
  bin4      NovaSeq-style 4 levels: 0-2->2, 3-14->12, 15-30->23, >=31->37
  bin2      two levels: <20->12, >=20->37 (a 'good / bad' flag)
  noqual    every base gets Q30 (equivalent to throwing qualities away; what CRAM
            lossy modes and many archives effectively do)

Modes
  --scheme S --out F       write a binned FASTQ (input may be .gz; output .gz if it ends with .gz)
  --profile TSV            write the quality-value histogram and the per-stream compressed size
                           (read names / bases / qualities, each zstd -19) for the input
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

PHRED_OFFSET = 33
SCHEMES = {
    "orig": None,
    "illumina8": [(0, 1, None), (2, 9, 6), (10, 19, 15), (20, 24, 22), (25, 29, 27), (30, 34, 33),
                  (35, 39, 37), (40, 93, 40)],
    "bin4": [(0, 2, 2), (3, 14, 12), (15, 30, 23), (31, 93, 37)],
    "bin2": [(0, 19, 12), (20, 93, 37)],
    "noqual": [(0, 93, 30)],
}


def translation_table(scheme: str) -> bytes | None:
    rules = SCHEMES[scheme]
    if rules is None:
        return None
    table = bytearray(range(256))
    for lo, hi, rep in rules:
        for q in range(lo, hi + 1):
            if rep is not None:
                table[q + PHRED_OFFSET] = rep + PHRED_OFFSET
    return bytes(table)


def open_any(path: Path, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, compresslevel=1) if "w" in mode else gzip.open(path, mode)
    return open(path, mode)


def rebin(src: Path, dst: Path, scheme: str) -> int:
    table = translation_table(scheme)
    n = 0
    with open_any(src, "rb") as fi, open_any(dst, "wb") as fo:
        for i, line in enumerate(fi):
            if i % 4 == 3 and table is not None:
                line = line.translate(table)
                n += 1
            fo.write(line)
    return n


def profile(src: Path, out: Path, threads: int) -> None:
    """Histogram of quality values + compressed size of each FASTQ stream on its own."""
    hist: Counter = Counter()
    with tempfile.TemporaryDirectory() as td:
        streams = {k: open(Path(td) / k, "wb") for k in ("names", "bases", "qualities")}
        with open_any(src, "rb") as fi:
            for i, line in enumerate(fi):
                k = i % 4
                if k == 0:
                    streams["names"].write(line)
                elif k == 1:
                    streams["bases"].write(line)
                elif k == 3:
                    streams["qualities"].write(line)
                    hist.update(line[:-1])
        for fh in streams.values():
            fh.close()
        rows = []
        for k in streams:
            p = Path(td) / k
            raw = p.stat().st_size
            comp = len(subprocess.run(["zstd", "-19", f"-T{threads}", "-q", "-c", str(p)],
                                      capture_output=True, check=True).stdout)
            rows.append(("stream", k, raw, comp))
    total_q = sum(hist.values())
    with open(out, "w") as fo:
        fo.write("kind\tkey\tvalue_a\tvalue_b\n")
        for kind, k, raw, comp in rows:
            fo.write(f"{kind}\t{k}\t{raw}\t{comp}\n")       # value_a raw bytes, value_b zstd-19 bytes
        for sym, c in sorted(hist.items()):
            fo.write(f"qual\t{sym - PHRED_OFFSET}\t{c}\t{c / total_q:.6f}\n")  # count, share


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", type=Path, required=True)
    ap.add_argument("--scheme", choices=SCHEMES)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--profile", type=Path)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    if a.profile:
        profile(a.src, a.profile, a.threads)
    if a.scheme:
        if not a.out:
            ap.error("--scheme needs --out")
        rebin(a.src, a.out, a.scheme)
    if not (a.profile or a.scheme):
        ap.error("nothing to do: give --scheme/--out and/or --profile")
    return 0


if __name__ == "__main__":
    sys.exit(main())
