"""Tidy the two price histories the crossover model needs (task 8).

Inputs are downloaded by the workflow (fetch_data.py records URL, size and checksums):
  NHGRI 'DNA Sequencing Costs: Data' (.xls): cost per raw megabase and per human-sized genome
      at NHGRI-funded sequencing centres, 2001-2022. Nominal USD, includes labour, reagents,
      instrument depreciation and data processing up to raw reads (not storage).
  Our World in Data 'historical cost of computer memory and storage' (.csv): disk-drive and SSD
      price per TB by year, constant 2020 USD.
  AWS S3 price list (.csv): today's list price of S3 Standard, the anchor for cloud storage.

Outputs
  sequencing_cost.csv  date, year (decimal), usd_per_mb, usd_per_genome
  storage_cost.csv     year, medium, usd2020_per_tb
  storage_anchor.tsv   today's S3 Standard price per GB-month with the price-list row it came from
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def decimal_year(d: pd.Series) -> pd.Series:
    start = pd.to_datetime(d.dt.year.astype(str) + "-01-01")
    end = pd.to_datetime((d.dt.year + 1).astype(str) + "-01-01")
    return d.dt.year + (d - start) / (end - start)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nhgri", type=Path, required=True)
    ap.add_argument("--owid", type=Path, required=True)
    ap.add_argument("--aws", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    x = pd.read_excel(a.nhgri)
    x.columns = ["date", "usd_per_mb", "usd_per_genome"]
    x["date"] = pd.to_datetime(x["date"])
    x = x.dropna().sort_values("date")
    x["year"] = decimal_year(x["date"]).round(4)
    if (x["usd_per_mb"] <= 0).any():
        raise SystemExit("non-positive sequencing cost in the NHGRI table")
    x.to_csv(a.outdir / "sequencing_cost.csv", index=False)

    o = pd.read_csv(a.owid)
    o = o[o["entity"] == "World"].rename(columns={"ddrives": "disk drives", "ssd": "SSD", "flash": "flash", "memory": "DRAM"})
    long = o.melt(id_vars="year", value_vars=["disk drives", "SSD", "flash", "DRAM"], var_name="medium",
                  value_name="usd2020_per_tb").dropna()
    long.sort_values(["medium", "year"]).to_csv(a.outdir / "storage_cost.csv", index=False)

    pl = pd.read_csv(a.aws, skiprows=5, dtype=str)
    row = pl[pl["PriceDescription"].fillna("").str.endswith("per GB - first 50 TB / month of storage used")].iloc[0]
    pd.DataFrame([{"service": "AWS S3 Standard, us-east-1, first 50 TB", "usd_per_gb_month": float(row["PricePerUnit"]),
                   "effective_date": row["EffectiveDate"], "price_list_row": row["PriceDescription"]}]).to_csv(
        a.outdir / "storage_anchor.tsv", sep="\t", index=False)
    print(x.tail(3).to_string(index=False))
    print(long[long.medium == "disk drives"].tail(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
