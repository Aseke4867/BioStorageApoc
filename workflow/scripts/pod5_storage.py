"""Raw Nanopore signal: how big is it, and is it safe to delete after basecalling? (task 5)

  inspect   --pod5 F --sample N --ids OUT --out TSV
            read count, signal samples, sample rate, bytes per sample (VBZ-compressed vs raw
            int16), and a random subset of N read ids for GPU basecalling
  archive   --bucket URL --prefix NAME=PREFIX ... --out TSV
            total bytes of whole public ONT datasets, raw signal vs their basecalls
  summarize --stats TSV --fastq NAME=FILE ... --identity NAME=TSV ... --bench NAME=TSV ... --out TSV
            per basecalling model: bytes of signal per base, accuracy, GPU time
"""
from __future__ import annotations

import argparse
import gzip
import random
import re
import sys
from pathlib import Path

import pandas as pd
import requests


def inspect(pod5_path: Path, n_sample: int, ids_out: Path, out: Path, seed: int) -> None:
    import pod5
    ids, samples, rates = [], 0, set()
    lens = []
    with pod5.Reader(pod5_path) as r:
        for rec in r.reads():
            ids.append(str(rec.read_id))
            samples += rec.num_samples
            lens.append(rec.num_samples)
            rates.add(rec.run_info.sample_rate)
    fbytes = pod5_path.stat().st_size
    rng = random.Random(seed)
    pick = ids if n_sample <= 0 or n_sample >= len(ids) else rng.sample(ids, n_sample)
    ids_out.write_text("\n".join(pick) + "\n")
    s = pd.Series(lens)
    pd.DataFrame([{
        "file": pod5_path.name, "file_bytes": fbytes, "reads": len(ids), "signal_samples": samples,
        "sample_rate_hz": ",".join(map(str, sorted(rates))),
        "median_samples_per_read": float(s.median()),
        "raw_int16_bytes": 2 * samples, "vbz_bits_per_sample": 8 * fbytes / samples,
        "vbz_ratio_vs_int16": 2 * samples / fbytes, "reads_sampled_for_basecalling": len(pick),
    }]).to_csv(out, sep="\t", index=False)


def s3_total(bucket: str, prefix: str) -> tuple[int, int]:
    total = n = 0
    token = None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        xml = requests.get(bucket + "/", params=params, timeout=60).text
        sizes = [int(x) for x in re.findall(r"<Size>(\d+)</Size>", xml)]
        total += sum(sizes)
        n += len(sizes)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            return total, n
        token = m.group(1)


def fastq_stats(path: Path) -> tuple[int, int]:
    op = gzip.open if str(path).endswith(".gz") else open
    reads = bases = 0
    with op(path, "rb") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                reads += 1
                bases += len(line) - 1
    return reads, bases


def kv(items: list[str]) -> dict[str, str]:
    return dict(x.split("=", 1) for x in items)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["inspect", "archive", "summarize"])
    ap.add_argument("--pod5", type=Path)
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--ids", type=Path)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--bucket")
    ap.add_argument("--prefix", nargs="*", default=[])
    ap.add_argument("--stats", type=Path)
    ap.add_argument("--fastq", nargs="*", default=[])
    ap.add_argument("--identity", nargs="*", default=[])
    ap.add_argument("--bench", nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    a.out.parent.mkdir(parents=True, exist_ok=True)

    if a.mode == "inspect":
        inspect(a.pod5, a.sample, a.ids, a.out, a.seed)
    elif a.mode == "archive":
        rows = []
        for name, prefix in kv(a.prefix).items():
            b, n = s3_total(a.bucket, prefix)
            rows.append({"dataset": name, "prefix": prefix, "objects": n, "bytes": b})
        df = pd.DataFrame(rows)
        df.to_csv(a.out, sep="\t", index=False)
        print(df.to_string(index=False))
    else:
        st = pd.read_csv(a.stats, sep="\t").iloc[0]
        # signal bytes attributable to the basecalled subset = per-read share of the file
        frac = st["reads_sampled_for_basecalling"] / st["reads"]
        signal_bytes = st["file_bytes"] * frac
        ident, bench = kv(a.identity), kv(a.bench)
        rows = []
        for model, fq in kv(a.fastq).items():
            reads, bases = fastq_stats(Path(fq))
            row = {"model": model, "reads_called": reads, "bases_called": bases,
                   "fastq_gz_bytes": Path(fq).stat().st_size,
                   "signal_bytes_subset": signal_bytes,
                   "signal_bytes_per_base": signal_bytes / bases,
                   "fastq_gz_bytes_per_base": Path(fq).stat().st_size / bases}
            row["signal_to_fastq_gz"] = row["signal_bytes_per_base"] / row["fastq_gz_bytes_per_base"]
            if model in ident:
                i = pd.read_csv(ident[model], sep="\t").iloc[0]
                row.update({"median_identity": i["median_identity"], "mean_identity": i["mean_identity"],
                            "phred_of_error": i["phred_of_error"], "reads_mapped": i["reads_mapped_total"]})
            if model in bench:
                b = pd.read_csv(bench[model], sep="\t").iloc[0]
                row.update({"gpu_wall_s": b["s"], "samples_per_s":
                            st["signal_samples"] * frac / b["s"]})
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(a.out, sep="\t", index=False)
        print(df.round(4).T.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
