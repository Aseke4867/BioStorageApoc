"""Tidy NCBI's SRA growth statistics (task 2).

The input is NCBI's own daily series of the SRA database size (bases and bytes, all and
open-access), downloaded by the workflow. It is reduced to one row per month (the last day
of each month that has data) and checked for the problems a time series from an archive can
have: unparsable dates, non-monotone totals, and gaps.

Output: sra_growth.csv  date, year (decimal), bases, open_access_bases, bytes, open_access_bytes
        sra_growth_issues.tsv  what was found and how it was handled
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--issues", type=Path, required=True)
    a = ap.parse_args(argv)
    d = pd.read_csv(a.input)
    issues = []
    d["date"] = pd.to_datetime(d["date"], format="%m/%d/%Y", errors="coerce")
    issues.append({"issue": "unparsable dates (dropped)", "rows": int(d["date"].isna().sum())})
    d = d.dropna(subset=["date"]).sort_values("date")
    for col in ("bases", "bytes"):
        down = (d[col].diff() < 0).sum()
        issues.append({"issue": f"day-to-day decreases in {col} (kept: withdrawn or suppressed data)", "rows": int(down)})
    gaps = d["date"].diff().dt.days
    issues.append({"issue": "gaps longer than 31 days between reports", "rows": int((gaps > 31).sum())})
    m = d.groupby(d["date"].dt.to_period("M")).tail(1).copy()
    start = pd.to_datetime(m["date"].dt.year.astype(str) + "-01-01")
    end = pd.to_datetime((m["date"].dt.year + 1).astype(str) + "-01-01")
    m["year"] = (m["date"].dt.year + (m["date"] - start) / (end - start)).round(4)
    m = m[["date", "year", "bases", "open_access_bases", "bytes", "open_access_bytes"]]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(a.out, index=False)
    issues.append({"issue": "monthly rows kept", "rows": len(m)})
    pd.DataFrame(issues).to_csv(a.issues, sep="\t", index=False)
    last = m.iloc[-1]
    print(f"SRA: {len(m)} months, last {last['date'].date()}: {last['bases'] / 1e15:.1f} Pbases, "
          f"{last['bytes'] / 1e15:.1f} PB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
