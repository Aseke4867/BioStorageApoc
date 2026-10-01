"""Benchmark compression codecs on one sequencing dataset (tasks 4 and 5).

A dataset is one or two uncompressed FASTQ files (paired-end = two). Every codec is run as an
external command under /usr/bin/time, so wall time and peak memory are measured the same way
for all of them. Every result is verified: the decompressed output must have the same MD5 as
the input (for Spring: after dropping the redundant read-name copy on the '+' line), otherwise
the row is marked lossless=False and the script exits non-zero.

General-purpose codecs compress each file separately. Spring is FASTQ-aware and gets both mates
together (that is where its reordering gain comes from); -l switches it to long-read mode.

Output TSV columns
  dataset codec level threads files in_bytes out_bytes ratio bits_per_base
  comp_s decomp_s comp_MBps decomp_MBps comp_rss_MB decomp_rss_MB lossless
"""
from __future__ import annotations

import argparse
import hashlib
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

log = logging.getLogger("bench_compress")

# name -> (level label, compress template, decompress template, extension)
# {i} input file, {o} output file, {t} threads
CODECS = {
    "gzip-6":   ("6",  "gzip -6 -c {i} > {o}",             "gzip -dc {i} > {o}",            ".gz"),
    "pigz-9":   ("9",  "pigz -9 -p {t} -c {i} > {o}",      "pigz -dc -p {t} {i} > {o}",     ".gz"),
    "bzip2-9":  ("9",  "bzip2 -9 -c {i} > {o}",            "bzip2 -dc {i} > {o}",           ".bz2"),
    "xz-6":     ("6",  "xz -6 -T {t} -c {i} > {o}",        "xz -dc -T {t} {i} > {o}",       ".xz"),
    "zstd-3":   ("3",  "zstd -3 -T{t} -q -c {i} > {o}",    "zstd -dc -q {i} > {o}",         ".zst"),
    "zstd-19":  ("19", "zstd -19 -T{t} -q -c {i} > {o}",   "zstd -dc -q {i} > {o}",         ".zst"),
    "zstd-19-long": ("19 --long=27", "zstd -19 --long=27 -T{t} -q -c {i} > {o}",
                     "zstd -dc -q --long=27 {i} > {o}", ".zst"),
}
TIME_RE = {"wall": re.compile(r"Elapsed \(wall clock\) time.*: (.+)$"),
           "rss": re.compile(r"Maximum resident set size \(kbytes\): (\d+)")}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_plusline_normalised(path: Path) -> str:
    """MD5 with each record's optional '+' line text dropped. fasterq-dump repeats the read name
    there; Spring keeps names once and writes a bare '+'. That repeat carries no information,
    so for Spring 'lossless' means: identical names, bases and qualities."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for i, line in enumerate(fh):
            h.update(b"+\n" if i % 4 == 2 else line)
    return h.hexdigest()


def parse_wall(s: str) -> float:
    parts = [float(x) for x in s.strip().split(":")]
    return sum(p * 60 ** k for k, p in enumerate(reversed(parts)))


def timed(cmd: str) -> tuple[float, float]:
    """Run a shell command under GNU time; return (wall seconds, peak RSS MB)."""
    p = subprocess.run(["/usr/bin/time", "-v", "bash", "-o", "pipefail", "-c", cmd],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd}\n{p.stderr[-2000:]}")
    wall = rss = float("nan")
    for line in p.stderr.splitlines():
        if m := TIME_RE["wall"].search(line):
            wall = parse_wall(m.group(1))
        elif m := TIME_RE["rss"].search(line):
            rss = int(m.group(1)) / 1024
    return wall, rss


def count_bases(fastqs: list[Path]) -> int:
    """Total bases = sum of sequence-line lengths (every 4th line, offset 1). Streamed."""
    n = 0
    for f in fastqs:
        with open(f, "rb") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    n += len(line) - 1
    return n


def bench_general(name: str, files: list[Path], work: Path, threads: int) -> dict:
    level, ctpl, dtpl, ext = CODECS[name]
    comp_s = decomp_s = 0.0
    crss = drss = 0.0
    out_bytes, ok = 0, True
    for f in files:
        c = work / (f.name + ext)
        d = work / (f.name + ".roundtrip")
        w, r = timed(ctpl.format(i=f, o=c, t=threads))
        comp_s, crss = comp_s + w, max(crss, r)
        w, r = timed(dtpl.format(i=c, o=d, t=threads))
        decomp_s, drss = decomp_s + w, max(drss, r)
        out_bytes += c.stat().st_size
        ok &= md5(d) == md5(f)
        c.unlink()
        d.unlink()
    return {"codec": name, "level": level, "comp_s": comp_s, "decomp_s": decomp_s,
            "comp_rss_MB": crss, "decomp_rss_MB": drss, "out_bytes": out_bytes, "lossless": ok}


def bench_spring(files: list[Path], work: Path, threads: int, long_reads: bool) -> dict:
    arc = work / "data.spring"
    flag = " -l" if long_reads else ""
    inputs = " ".join(str(f) for f in files)
    w1, r1 = timed(f"cd {work} && spring -c -i {inputs} -o {arc} -t {threads} -w {work}{flag}")
    outs = [work / f"rt_{k}.fastq" for k in range(len(files))]
    w2, r2 = timed(f"cd {work} && spring -d -i {arc} -o {' '.join(map(str, outs))} -t {threads} -w {work}")
    ok = all(md5_plusline_normalised(o) == md5_plusline_normalised(f) for o, f in zip(outs, files))
    size = arc.stat().st_size
    for p in outs + [arc]:
        p.unlink()
    return {"codec": "spring" + ("-long" if long_reads else ""), "level": "lossless",
            "comp_s": w1, "decomp_s": w2, "comp_rss_MB": r1, "decomp_rss_MB": r2,
            "out_bytes": size, "lossless": ok}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", type=Path, nargs="+", required=True, help="uncompressed FASTQ file(s)")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--codecs", default=",".join(CODECS) + ",spring")
    ap.add_argument("--long-reads", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workdir", type=Path, default=None, help="scratch dir on a fast local disk")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    in_bytes = sum(f.stat().st_size for f in a.fastq)
    bases = count_bases(a.fastq)
    work = Path(tempfile.mkdtemp(prefix="bench_", dir=a.workdir))
    rows = []
    try:
        for name in a.codecs.split(","):
            log.info("%s: %s", a.dataset, name)
            if name == "spring":
                r = bench_spring(a.fastq, work, a.threads, a.long_reads)
            else:
                r = bench_general(name, a.fastq, work, a.threads)
            rows.append(r)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    df = pd.DataFrame(rows)
    df.insert(0, "dataset", a.dataset)
    df["threads"] = a.threads
    df["files"] = len(a.fastq)
    df["in_bytes"] = in_bytes
    df["bases"] = bases
    df["ratio"] = in_bytes / df["out_bytes"]
    df["bits_per_base"] = 8 * df["out_bytes"] / bases
    df["comp_MBps"] = in_bytes / 1e6 / df["comp_s"]
    df["decomp_MBps"] = in_bytes / 1e6 / df["decomp_s"]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.round(4).to_csv(a.out, sep="\t", index=False)
    log.info("\n%s", df[["codec", "ratio", "bits_per_base", "comp_MBps", "decomp_MBps", "lossless"]]
             .round(2).to_string(index=False))
    return 0 if df["lossless"].all() else 4


if __name__ == "__main__":
    sys.exit(main())
