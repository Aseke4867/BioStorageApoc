"""When does storing sequence data cost more than producing it? (task 8, the headline projection)

The comparison is made per sequenced megabase, which is the unit both prices share:
  sequencing   NHGRI cost per raw Mb in year t (fitted trend)
  storage      cost of keeping that Mb for R years from year t:
               bytes per base x copies x sum over the R years of the price per GB-month,
               where the price starts at today's AWS S3 Standard list price and moves with the
               long-run decline of disk prices (Our World in Data, fitted trend).
Crossover = first year in which the storage cost of a newly sequenced Mb reaches its sequencing
cost: from then on, re-sequencing a sample (if it still exists) is cheaper than archiving its reads.

Two readings of the NHGRI series are both defensible and give different futures, so both are
kept: the whole second-generation era (2008-2022) and the recent plateau (2015-2022). Monte Carlo
draws pick a reading at random and sample its fitted parameters and the storage decline rate by
residual bootstrap, so the reported range covers model choice as well as fit uncertainty.

Scenarios change only what part B measured: bits per stored base (codec, quality binning) and the
storage price (tiering). Output tables: price_fits.tsv, crossover_summary.tsv,
crossover_curves.tsv, crossover_sensitivity.tsv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

T0 = 2000.0


def ols_boot(t, y, draws, rng):
    """Slope/intercept of y ~ t plus residual-bootstrap draws of both."""
    X = np.c_[np.ones_like(t), t - T0]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ beta
    out = np.empty((draws, 2))
    for i in range(draws):
        yb = X @ beta + rng.choice(res, len(res), replace=True)
        out[i] = np.linalg.lstsq(X, yb, rcond=None)[0]
    return beta, out


def bits_from(spec: str, tables: Path) -> float:
    if spec == "sra_archive":
        # what the SRA itself stores per base: bytes on disk / bases, latest month of NCBI's statistics
        last = pd.read_csv(tables / "sra_growth.csv").iloc[-1]
        return 8 * float(last["bytes"]) / float(last["bases"])
    file, codec = spec.split(":")
    d = pd.read_csv(tables / file, sep="\t")
    return float(d.loc[d["codec"] == codec, "bits_per_base"].iloc[0])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", type=Path, required=True, help="results/tables (inputs from both parts)")
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args(argv)
    cfg = json.loads(a.config_json)
    rng = np.random.default_rng(7)
    D = int(cfg["draws"])
    grid = np.arange(cfg["first_year"], cfg["last_year"] + 1e-9, 1 / 12)

    seq = pd.read_csv(a.tables / "sequencing_cost.csv")
    sto = pd.read_csv(a.tables / "storage_cost.csv")
    anchor = pd.read_csv(a.tables / "storage_anchor.tsv", sep="\t").iloc[0]
    anchor_price = float(anchor["usd_per_gb_month"])
    anchor_year = pd.to_datetime(anchor["effective_date"]).year + pd.to_datetime(anchor["effective_date"]).dayofyear / 365.25

    # ---- price trends
    fits, seq_draws = [], {}
    for name, (y0, y1) in cfg["sequencing_fits"].items():
        s = seq[(seq.year >= y0) & (seq.year <= y1 + 1)]
        beta, bd = ols_boot(s.year.to_numpy(), np.log10(s.usd_per_mb.to_numpy()), D, rng)
        seq_draws[name] = bd
        fits.append({"series": "sequencing cost per Mb (NHGRI)", "fit": name, "years": f"{y0}-{y1}", "points": len(s),
                     "annual_change_pct": 100 * (10 ** beta[1] - 1),
                     "annual_change_p05": 100 * (10 ** np.quantile(bd[:, 1], 0.05) - 1),
                     "annual_change_p95": 100 * (10 ** np.quantile(bd[:, 1], 0.95) - 1),
                     "halving_years": -np.log10(2) / beta[1]})
    y0, y1 = cfg["storage_fit_years"]
    dd = sto[(sto.medium == "disk drives") & sto.year.between(y0, y1)]
    sbeta, sbd = ols_boot(dd.year.to_numpy(float), np.log10(dd.usd2020_per_tb.to_numpy()), D, rng)
    fits.append({"series": "disk price per TB (OWID, 2020 USD)", "fit": "storage", "years": f"{y0}-{y1}", "points": len(dd),
                 "annual_change_pct": 100 * (10 ** sbeta[1] - 1),
                 "annual_change_p05": 100 * (10 ** np.quantile(sbd[:, 1], 0.05) - 1),
                 "annual_change_p95": 100 * (10 ** np.quantile(sbd[:, 1], 0.95) - 1),
                 "halving_years": -np.log10(2) / sbeta[1]})
    pd.DataFrame(fits).to_csv(a.outdir / "price_fits.tsv", sep="\t", index=False)

    # ---- scenarios
    scen = {}
    for key, sc in cfg["scenarios"].items():
        bits = bits_from(sc["bits_per_base_from"], a.tables)
        saving = sc.get("tier_saving", 0.0)
        if "tier_saving_from" in sc:
            saving = float(pd.read_csv(a.tables / sc["tier_saving_from"], sep="\t")["storage_cost_saving_fraction"].iloc[0])
        scen[key] = {"label": sc["label"], "bytes_per_base": bits / 8, "price_factor": 1 - saving}

    names = list(cfg["sequencing_fits"])
    months = np.arange(int(12 * cfg["retention_years"])) / 12

    def run(sc, draw_window, retention_months=months, copies=cfg["copies"]):
        """Return crossover years (nan = none before last_year) and the ratio curves for D draws."""
        cross = np.full(D, np.nan)
        ratios = np.empty((D, len(grid)))
        for i in range(D):
            w = draw_window[i]
            a_s, b_s = seq_draws[w][i]
            seq_mb = 10 ** (a_s + b_s * (grid - T0))
            slope = sbd[i, 1]
            # price per GB-month at time u: anchor x decline since the anchor date
            start = anchor_price * sc["price_factor"] * 10 ** (slope * (grid - anchor_year))
            # sum over retention months of price(t + m) = start(t) * sum 10^(slope*m)
            factor = np.sum(10 ** (slope * retention_months))
            store_mb = sc["bytes_per_base"] * 1e6 / 1e9 * copies * start * factor
            r = store_mb / seq_mb
            ratios[i] = r
            hit = np.nonzero(r >= 1)[0]
            if hit.size:
                cross[i] = grid[hit[0]]
        return cross, ratios

    def summarise(cross, label, **extra):
        ok = cross[~np.isnan(cross)]
        q = lambda p: float(np.quantile(ok, p)) if ok.size else np.nan  # noqa: E731
        return {"scenario": label, **extra,
                "p_cross_by_2030": float(np.mean(cross <= 2030)), "p_cross_by_2045": float(np.mean(cross <= 2045)),
                f"p_no_cross_by_{cfg['last_year']}": float(np.mean(np.isnan(cross))),
                "median_year_if_crossing": q(0.5), "p05_year": q(0.05), "p95_year": q(0.95),
                # unconditional median: counts 'never' as later than everything
                "median_year_all_draws": float(np.nanquantile(np.where(np.isnan(cross), np.inf, cross), 0.5))}

    mix = rng.choice(names, size=D)
    summary, curves = [], []
    for key, sc in scen.items():
        for wlabel, wdraw in [("both readings", mix)] + [(n, np.full(D, n)) for n in names]:
            cross, ratios = run(sc, wdraw)
            summary.append(summarise(cross, sc["label"], key=key, sequencing_trend=wlabel,
                                     bytes_per_base=round(sc["bytes_per_base"], 4), price_factor=round(sc["price_factor"], 4),
                                     ratio_in_first_year_median=float(np.median(ratios[:, 0]))))
            if wlabel == "both readings":
                for j in range(0, len(grid), 12):
                    curves.append({"scenario": sc["label"], "key": key, "year": round(grid[j], 2),
                                   "ratio_median": float(np.median(ratios[:, j])),
                                   "ratio_p05": float(np.quantile(ratios[:, j], 0.05)),
                                   "ratio_p95": float(np.quantile(ratios[:, j], 0.95))})
    pd.DataFrame(summary).to_csv(a.outdir / "crossover_summary.tsv", sep="\t", index=False)
    pd.DataFrame(curves).to_csv(a.outdir / "crossover_curves.tsv", sep="\t", index=False)

    # ---- sensitivity of the current-practice date to retention and number of copies
    sens = []
    base = scen["current"]
    for R in (5, 10, 20):
        for copies in (1, 3):
            cross, _ = run(base, mix, retention_months=np.arange(12 * R) / 12, copies=copies)
            sens.append(summarise(cross, base["label"], retention_years=R, copies=copies))
    pd.DataFrame(sens).to_csv(a.outdir / "crossover_sensitivity.tsv", sep="\t", index=False)

    print(pd.DataFrame(fits).round(2).to_string(index=False))
    show = ["scenario", "sequencing_trend", "bytes_per_base", "ratio_in_first_year_median", "p_cross_by_2030",
            "p_cross_by_2045", "median_year_all_draws", "p05_year", "p95_year"]
    print(pd.DataFrame(summary)[show].round(3).to_string(index=False))
    print(pd.DataFrame(sens).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
