"""Instrument throughput by generation, and what the archive says platforms actually produce (task 3).

1. config/instruments.tsv: vendor-stated maximum output per run and run time for 15 instruments
   (2006-2023), each with its source. Turned into Gb per instrument-day, and a log-linear fit of
   capacity against launch year gives the capacity doubling time per platform family.
2. ENA Portal API: number of public read runs per platform per release year, all organisms
   (one /count call per platform and year). Vendor capacity says what an instrument *can* do;
   run counts say which instruments the community actually used.

Outputs: instruments.tsv, instrument_capacity_fit.tsv, platform_runs_by_year.tsv
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ENA_COUNT = "https://www.ebi.ac.uk/ena/portal/api/count"
PLATFORMS = ["ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT", "ION_TORRENT", "LS454", "ABI_SOLID",
             "BGISEQ", "DNBSEQ", "ELEMENT", "ULTIMA", "CAPILLARY"]
UA = {"User-Agent": "aitu-bioinfo-course-project/1.0"}


def count(query: str, retries: int = 4) -> int:
    for k in range(retries):
        try:
            r = requests.get(ENA_COUNT, params={"result": "read_run", "query": query}, timeout=120, headers=UA)
            r.raise_for_status()
            return int(r.text.split()[-1])
        except (requests.RequestException, ValueError):
            time.sleep(2 ** k)
    raise RuntimeError(f"ENA count failed: {query}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instruments", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--first-year", type=int, default=2008)
    ap.add_argument("--last-year", type=int, default=2025, help="last complete year")
    ap.add_argument("--runs-table", type=Path, help="use this cached table instead of querying ENA")
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    ins = pd.read_csv(a.instruments, sep="\t", comment="#")
    ins["gb_per_day"] = ins["max_output_gb_per_run"] / (ins["run_time_h"] / 24)
    ins.to_csv(a.outdir / "instruments.tsv", sep="\t", index=False)
    fits = []
    for fam, g in [("Illumina", ins[ins.platform == "Illumina"]), ("all platforms", ins)]:
        b, c = np.polyfit(g["launch_year"], np.log2(g["gb_per_day"]), 1)
        fits.append({"family": fam, "instruments": len(g), "log2_gb_day_per_year": round(b, 4),
                     "capacity_doubling_years": round(1 / b, 2),
                     "first": int(g["launch_year"].min()), "last": int(g["launch_year"].max())})
    pd.DataFrame(fits).to_csv(a.outdir / "instrument_capacity_fit.tsv", sep="\t", index=False)

    if a.runs_table:
        runs = pd.read_csv(a.runs_table, sep="\t")
    else:
        jobs = [(p, y) for y in range(a.first_year, a.last_year + 1) for p in PLATFORMS + ["ALL"]]

        def one(job):
            p, y = job
            q = f"first_public>={y}-01-01 AND first_public<={y}-12-31"
            if p != "ALL":
                q = f'instrument_platform="{p}" AND ' + q
            return {"year": y, "platform": p, "runs": count(q)}

        with ThreadPoolExecutor(8) as ex:
            runs = pd.DataFrame(list(ex.map(one, jobs)))
    runs.to_csv(a.outdir / "platform_runs_by_year.tsv", sep="\t", index=False)
    piv = runs.pivot(index="year", columns="platform", values="runs")
    print(pd.DataFrame(fits).to_string(index=False))
    print((piv.div(piv["ALL"], axis=0) * 100).round(1).drop(columns="ALL").tail(5).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
