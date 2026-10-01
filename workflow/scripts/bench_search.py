"""Is a compressed archive still searchable, and what does a query cost? (task 6)

The same Illumina run is stored six ways, and the same questions are asked of each:

  format          how a query is answered
  fastq.gz        stream-decompress everything, scan for the gene's k-mers (seqkit grep)
  fastq.zst       same, zstd stream
  spring          decompress the whole archive to disk first (no random access), then scan
  BAM / CRAM      the alignment IS the index: ask for the gene's coordinates (samtools view)
  BLAST db        2-bit packed sequence + word index: blastn the gene against the reads

Queries
  gene    reads that come from one gene (lacZ). Recall is measured against the BAM region
          answer, so speed is never reported without accuracy.
  random  N random 1 kb windows, per-query latency (BAM vs CRAM only; FASTQ formats have no
          coordinates at all, which is itself the result)
  scan    count every read (full decompression throughput)

Also records bytes on disk and bits per base for every format.
"""
from __future__ import annotations

import argparse
import random
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd


def sh(cmd: str) -> tuple[float, str]:
    t0 = time.perf_counter()
    p = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    dt = time.perf_counter() - t0
    if p.returncode != 0:
        raise RuntimeError(f"failed: {cmd}\n{p.stderr[-1500:]}")
    return dt, p.stdout


def names_from_hits(text: str) -> set[str]:
    """Lines are 'header<TAB>sequence' (FASTQ records flattened with paste)."""
    return {l[1:].split()[0] for l in text.splitlines() if l.startswith("@")}


def kmers(seq: str, k: int, step: int) -> list[str]:
    """k-mers of the gene on both strands (reads come from either strand)."""
    fwd = [seq[i:i + k] for i in range(0, len(seq) - k + 1, step)]
    comp = str.maketrans("ACGT", "TGCA")
    return fwd + [x.translate(comp)[::-1] for x in fwd]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-gz", nargs=2, required=True)
    ap.add_argument("--fastq-zst", nargs=2, required=True)
    ap.add_argument("--spring", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--cram", required=True)
    ap.add_argument("--cram-archive", required=True)
    ap.add_argument("--blastdb", required=True, help="BLAST db prefix")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--bases", type=int, required=True, help="sequenced bases, for bits/base")
    ap.add_argument("--random-regions", type=int, default=200)
    ap.add_argument("--region-len", type=int, default=1000)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--out", type=Path, required=True, help="query results TSV")
    ap.add_argument("--sizes", type=Path, required=True, help="bytes per format TSV")
    a = ap.parse_args(argv)
    t = a.threads
    rows, sizes = [], []

    def size_of(fmt, paths):
        b = sum(Path(p).stat().st_size for p in paths)
        sizes.append({"format": fmt, "bytes": b, "bits_per_base": 8 * b / a.bases})

    blast_files = list(Path(a.blastdb).parent.glob(Path(a.blastdb).name + ".*"))
    size_of("fastq.gz (as archived)", a.fastq_gz)
    size_of("fastq.zst (-19)", a.fastq_zst)
    size_of("spring", [a.spring])
    size_of("BAM", [a.bam])
    size_of("CRAM 3.1", [a.cram])
    size_of("CRAM 3.1 archive", [a.cram_archive])
    size_of("BLAST db", blast_files)

    # ---------------- gene query
    _, gene_fa = sh(f"samtools faidx {a.ref} {a.region}")
    gene = "".join(gene_fa.splitlines()[1:]).upper()
    work = Path(tempfile.mkdtemp(dir=a.workdir))
    (work / "gene.fa").write_text(f">gene\n{gene}\n")
    (work / "kmers.txt").write_text("\n".join(kmers(gene, 31, 20)) + "\n")

    def rec(fmt, query, times, found, truth=None):
        r = {"format": fmt, "query": query, "median_s": statistics.median(times), "min_s": min(times),
             "n_runs": len(times), "reads_found": len(found) if found is not None else None}
        if truth is not None and found is not None:
            r["recall_vs_bam"] = len(found & truth) / len(truth) if truth else float("nan")
            r["extra_vs_bam"] = len(found - truth)
        rows.append(r)

    def repeat(cmd):
        out, ts = None, []
        for _ in range(a.repeats):
            dt, out = sh(cmd)
            ts.append(dt)
        return ts, out

    ts, out = repeat(f"samtools view -F 0x904 {a.bam} {a.region} | cut -f1")
    truth = set(out.split())
    rec("BAM", "gene", ts, truth, truth)
    ts, out = repeat(f"samtools view -F 0x904 -T {a.ref} {a.cram} {a.region} | cut -f1")
    rec("CRAM 3.1", "gene", ts, set(out.split()), truth)
    ts, out = repeat(f"samtools view -F 0x904 -T {a.ref} {a.cram_archive} {a.region} | cut -f1")
    rec("CRAM 3.1 archive", "gene", ts, set(out.split()), truth)
    # exact k-mer scan: GNU grep -F matches all 2x154 patterns in one pass (Aho-Corasick style)
    grep = f"paste - - - - | cut -f1,2 | (grep -F -f {work}/kmers.txt || true)"
    ts, out = repeat(f"cat {' '.join(a.fastq_gz)} | pigz -dc -p {t} | {grep}")
    rec("fastq.gz (as archived)", "gene", ts, names_from_hits(out), truth)
    ts, out = repeat(f"zstd -dc -q {' '.join(a.fastq_zst)} | {grep}")
    rec("fastq.zst (-19)", "gene", ts, names_from_hits(out), truth)
    sp = work / "sp"
    sp.mkdir()
    ts, out = repeat(f"cd {sp} && spring -d -i {a.spring} -o {sp}/r1.fq {sp}/r2.fq -t {t} -w {sp} >/dev/null "
                     f"&& cat {sp}/r1.fq {sp}/r2.fq | {grep}; rm -f {sp}/r1.fq {sp}/r2.fq")
    rec("spring", "gene", ts, names_from_hits(out), truth)
    ts, out = repeat(f"blastn -query {work}/gene.fa -db {a.blastdb} -outfmt '6 stitle' -evalue 1e-20 "
                     f"-max_target_seqs 1000000 -num_threads {t} | sort -u")
    # subject title is the FASTQ header ('SRR... 1/1'); both mates share the read name
    found = {line.split()[0] for line in out.splitlines() if line.strip()}
    rec("BLAST db", "gene", ts, found, truth)

    # ---------------- random region latency (indexed formats only)
    _, fai = sh(f"cut -f1,2 {a.ref}.fai")
    chrom, length = fai.split()[0], int(fai.split()[1])
    rng = random.Random(11)
    regions = [f"{chrom}:{s}-{s + a.region_len}" for s in
               (rng.randrange(1, length - a.region_len) for _ in range(a.random_regions))]
    for fmt, path, extra in [("BAM", a.bam, ""), ("CRAM 3.1", a.cram, f"-T {a.ref}"),
                             ("CRAM 3.1 archive", a.cram_archive, f"-T {a.ref}")]:
        lat = [sh(f"samtools view -c {extra} {path} {r}")[0] for r in regions]
        rows.append({"format": fmt, "query": f"random {a.region_len} bp region",
                     "median_s": statistics.median(lat), "min_s": min(lat), "n_runs": len(lat),
                     "p95_s": sorted(lat)[int(0.95 * len(lat))]})

    # ---------------- full scan
    for fmt, cmd in [
        ("BAM", f"samtools view -c -@ {t} {a.bam}"),
        ("CRAM 3.1", f"samtools view -c -@ {t} -T {a.ref} {a.cram}"),
        ("CRAM 3.1 archive", f"samtools view -c -@ {t} -T {a.ref} {a.cram_archive}"),
        ("fastq.gz (as archived)", f"cat {' '.join(a.fastq_gz)} | pigz -dc -p {t} | wc -l"),
        ("fastq.zst (-19)", f"zstd -dc -q {' '.join(a.fastq_zst)} | wc -l"),
        ("spring", f"cd {sp} && spring -d -i {a.spring} -o {sp}/r1.fq {sp}/r2.fq -t {t} -w {sp} >/dev/null "
                   f"&& cat {sp}/r1.fq {sp}/r2.fq | wc -l; rm -f {sp}/r1.fq {sp}/r2.fq"),
    ]:
        ts, _ = repeat(cmd)
        rows.append({"format": fmt, "query": "full scan", "median_s": statistics.median(ts),
                     "min_s": min(ts), "n_runs": len(ts)})

    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).round(4).to_csv(a.out, sep="\t", index=False)
    pd.DataFrame(sizes).round(4).to_csv(a.sizes, sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print(pd.DataFrame(sizes).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
