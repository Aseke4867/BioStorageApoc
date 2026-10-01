"""Model the growth of the archives and test how much the forecast depends on the model (task 2).

Series (all real, from the archives' own statistics):
  GenBank traditional bases and WGS bases, release by release (gbrel.txt, 1982-2026)
  SRA bases and bytes, monthly (NCBI SRA statistics, 2007-2024)

Three functional forms, all fitted to log10(size) by least squares over a recent window
(config growth.fit_from; the archives' start-up years grew 5x a year and would dominate):
  exponential     log y = a + b t                    constant doubling time
  linear          y = a + b t                        constant yearly intake (relative growth slows)
  logistic        y = K / (1 + exp(-r (t - t0)))     growth saturates at K
(A quadratic in log y was tried first and rejected: it fits well in-sample but turns downwards,
i.e. predicts archives shrinking, as soon as it is extrapolated.)

Each model is judged two ways: in-sample (AIC on log10 residuals) and by a backtest that hides
the last N years, refits, and compares the forecast with what actually happened. Projections to
the forecast year come with a 5-95% band from a residual bootstrap (parameter uncertainty only).

Also: annual SRA growth divided by the capacity of the newest Illumina flagship gives the
number of flagship instruments, running all year, that the archive's intake corresponds to.

Outputs: growth_fits.tsv, growth_projections.tsv, growth_backtest.tsv, sra_intake_vs_capacity.tsv
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit

warnings.simplefilter("ignore", OptimizeWarning)
T0 = 2000.0


def f_exp(t, a, b):
    return a + b * (t - T0)


def f_linear(t, a, b):
    """log10 of a straight line in linear space; a = size at T0 (in units of 1e12), b per year."""
    return np.log10(np.maximum(a + b * (t - T0), 1e-6)) + 12


def f_logistic(t, logK, r, t0):
    return logK - np.log10(1 + np.exp(-r * (t - t0)))


MODELS = {"exponential": f_exp, "linear": f_linear, "logistic": f_logistic}


def p0(name, t, y):
    b = np.polyfit(t - T0, y, 1)
    if name == "exponential":
        return list(b[::-1])
    if name == "linear":
        yl = 10 ** (y - 12)
        bl = np.polyfit(t - T0, yl, 1)
        return [bl[1], max(bl[0], 1e-9)]
    return [y.max() + 1.0, b[0] * np.log(10), t.max()]


def bounds(name, y):
    if name == "logistic":
        return ([y.max() - 0.5, 1e-3, 1980], [y.max() + 6, 5.0, 2100])
    if name == "linear":
        return ([-np.inf, 0.0], [np.inf, np.inf])        # an archive does not shrink
    return (-np.inf, np.inf)


def fit(name, t, y):
    if name == "linear":
        # constant intake is estimated from absolute sizes (ordinary least squares in linear space);
        # a log-space fit lets the early, small values pull the line below the recent data
        b, a = np.polyfit(t - T0, 10 ** (y - 12), 1)
        return np.array([a, max(b, 0.0)])
    f = MODELS[name]
    p, _ = curve_fit(f, t, y, p0=p0(name, t, y), bounds=bounds(name, y), maxfev=50000)
    return p


def doubling_years(name, p, t):
    """Doubling time at year t from the local slope of log10 y."""
    h = 1e-3
    s = (MODELS[name](t + h, *p) - MODELS[name](t - h, *p)) / (2 * h)   # log10 per year
    return np.log10(2) / s if s > 0 else np.inf


def load_series(genbank: Path, sra: Path, fit_from: dict) -> dict[str, pd.DataFrame]:
    g = pd.read_csv(genbank).dropna(subset=["bases", "date"])
    g["year"] = pd.to_datetime(g["date"]).dt.year + (pd.to_datetime(g["date"]).dt.dayofyear - 1) / 365.25
    out = {}
    for comp, name in [("GenBank", "GenBank traditional (bases)"), ("WGS", "GenBank WGS (bases)")]:
        x = g[(g.component == comp) & (g.bases > 0) & (g.year >= fit_from.get("genbank", 0))]
        out[name] = pd.DataFrame({"year": x["year"], "value": x["bases"].astype(float)})
    s = pd.read_csv(sra)
    s = s[s["year"] >= fit_from.get("sra", 0)]
    out["SRA (bases)"] = pd.DataFrame({"year": s["year"], "value": s["bases"].astype(float)})
    out["SRA (bytes on disk)"] = pd.DataFrame({"year": s["year"], "value": s["bytes"].astype(float)})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genbank", type=Path, required=True,
                    help="parsed GenBank CSV, or the folder holding genbank_growth_gbrel_<release>.csv (newest is used)")
    ap.add_argument("--sra", type=Path, required=True)
    ap.add_argument("--instruments", type=Path, required=True)
    ap.add_argument("--end-year", type=int, default=2045)
    ap.add_argument("--backtest-years", type=int, default=5)
    ap.add_argument("--bootstrap", type=int, default=300)
    ap.add_argument("--fit-from", default='{"genbank": 2012, "sra": 2014}', help="JSON: first year per archive")
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    if a.genbank.is_dir():
        a.genbank = sorted(a.genbank.glob("genbank_growth_gbrel_*.csv"),
                           key=lambda p: float(p.stem.rsplit("_", 1)[1]))[-1]

    fits, proj, back = [], [], []
    import json
    for series, df in load_series(a.genbank, a.sra, json.loads(a.fit_from)).items():
        t, y = df["year"].to_numpy(), np.log10(df["value"].to_numpy())
        tmax = t.max()
        years = np.arange(np.floor(tmax) + 1, a.end_year + 1)
        for name, f in MODELS.items():
            try:
                p = fit(name, t, y)
            except RuntimeError:
                fits.append({"series": series, "model": name, "note": "did not converge"})
                continue
            res = y - f(t, *p)
            n, k = len(y), len(p)
            rss = float(np.sum(res ** 2))
            aic = n * np.log(rss / n) + 2 * k
            # backtest: hide the last N years
            m = t <= tmax - a.backtest_years
            try:
                pb = fit(name, t[m], y[m])
                pred = f(t[~m], *pb)
                err = pred - y[~m]
                bt_last = float(10 ** err[-1])
                bt_mae = float(np.mean(np.abs(err)))
            except (RuntimeError, ValueError):
                bt_last = bt_mae = np.nan
            back.append({"series": series, "model": name, "train_until": round(tmax - a.backtest_years, 2),
                         "forecast_over_actual_at_end": bt_last, "mean_abs_log10_error": bt_mae})
            # residual bootstrap for the projection band
            boots = []
            lin_fit = 10 ** (f(t, *p) - 12)
            lin_res = 10 ** (y - 12) - lin_fit
            for _ in range(a.bootstrap):
                if name == "linear":      # resample in the space the model is fitted in
                    yb = np.log10(np.maximum(lin_fit + rng.choice(lin_res, size=len(t), replace=True), 1e-6)) + 12
                else:
                    yb = f(t, *p) + rng.choice(res, size=len(res), replace=True)
                try:
                    boots.append(f(years, *fit(name, t, yb)))
                except RuntimeError:
                    continue
            boots = np.array(boots)
            central = f(years, *p)
            for i, yr in enumerate(years):
                proj.append({"series": series, "model": name, "year": int(yr), "value": 10 ** central[i],
                             "p05": 10 ** np.quantile(boots[:, i], 0.05), "p95": 10 ** np.quantile(boots[:, i], 0.95)})
            v_end = 10 ** central[-1]
            fits.append({"series": series, "model": name, "n_points": n, "first": round(t.min(), 2),
                         "last": round(tmax, 2), "params": ";".join(f"{x:.6g}" for x in p),
                         "rmse_log10": float(np.sqrt(rss / n)), "aic": aic,
                         "doubling_years_at_last_point": doubling_years(name, p, tmax),
                         f"projected_{a.end_year}": v_end,
                         f"projected_{a.end_year}_over_last_observed": v_end / 10 ** y[-1]})
    fits = pd.DataFrame(fits)
    fits.to_csv(a.outdir / "growth_fits.tsv", sep="\t", index=False)
    pd.DataFrame(proj).to_csv(a.outdir / "growth_projections.tsv", sep="\t", index=False)
    pd.DataFrame(back).to_csv(a.outdir / "growth_backtest.tsv", sep="\t", index=False)

    # archive intake vs instrument capacity
    s = pd.read_csv(a.sra)
    s["y"] = s["year"].astype(float).astype(int)
    annual = s.groupby("y")["bases"].max().diff().dropna()
    ins = pd.read_csv(a.instruments, sep="\t")
    il = ins[ins.platform == "Illumina"].sort_values("launch_year")
    rows = []
    for yr, added in annual.items():
        avail = il[il.launch_year <= yr]
        if avail.empty or yr >= int(s["year"].max()):     # the last year is incomplete
            continue
        top = avail.iloc[-1]
        cap = top["gb_per_day"] * 1e9 * 365
        rows.append({"year": int(yr), "sra_bases_added": added, "flagship": top["instrument"],
                     "flagship_bases_per_year": cap, "flagship_instrument_years": added / cap})
    pd.DataFrame(rows).to_csv(a.outdir / "sra_intake_vs_capacity.tsv", sep="\t", index=False)
    print(fits[["series", "model", "aic", "rmse_log10", "doubling_years_at_last_point",
                f"projected_{a.end_year}"]].round(3).to_string(index=False))
    print(pd.DataFrame(back).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
