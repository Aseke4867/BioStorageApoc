"""Figures for sections 2, 3 and the headline projection (student A). Reads results/tables.

A1 archive growth with the three model projections; A2 instrument capacity and platform mix;
A3 the two price histories with their fitted trends; A4 storage-to-sequencing cost ratio and the
crossover year per scenario. Axis labels carry units; captions are in report/report.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 10, "axes.titleweight": "bold",
    "legend.frameon": False,
})
MODEL_COLOR = {"exponential": C[1], "linear": C[2], "logistic": C[6]}


def save(fig, out, name):
    fig.tight_layout()
    fig.savefig(out / name, facecolor="white")
    plt.close(fig)
    print("wrote", name)


def fig_growth(t: Path, out: Path):
    proj = pd.read_csv(t / "growth_projections.tsv", sep="\t")
    sra = pd.read_csv(t / "sra_growth.csv")
    gb = sorted(t.glob("genbank_growth_gbrel_*.csv"))[-1]
    g = pd.read_csv(gb).dropna(subset=["bases"])
    g["year"] = pd.to_datetime(g["date"]).dt.year + (pd.to_datetime(g["date"]).dt.dayofyear - 1) / 365.25
    panels = [("GenBank traditional (bases)", g[g.component == "GenBank"][["year", "bases"]].values),
              ("GenBank WGS (bases)", g[g.component == "WGS"][["year", "bases"]].values),
              ("SRA (bases)", sra[["year", "bases"]].values)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, (series, obs) in zip(axes, panels):
        ax.plot(obs[:, 0], obs[:, 1], color=C[0], lw=2, label="observed")
        for m, col in MODEL_COLOR.items():
            p = proj[(proj.series == series) & (proj.model == m)]
            if p.empty:
                continue
            ax.plot(p.year, p.value, color=col, lw=1.8, ls="--", label=f"{m} fit")
            ax.fill_between(p.year, p.p05, p.p95, color=col, alpha=0.12, lw=0)
        ax.set_yscale("log")
        ax.set_xlabel("year")
        ax.set_ylabel("bases in the archive (log)")
        ax.set_title(series.replace(" (bases)", ""))
        ax.set_xlim(max(obs[:, 0].min(), 2000), 2046)
    axes[0].legend(fontsize=7.5, loc="upper left")
    save(fig, out, "fig_a1_archive_growth.png")


def fig_instruments(t: Path, out: Path):
    ins = pd.read_csv(t / "instruments.tsv", sep="\t")
    runs = pd.read_csv(t / "platform_runs_by_year.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    markers = {"Illumina": "o", "PacBio": "s", "Oxford Nanopore": "^"}
    for k, (plat, mk) in enumerate(markers.items()):
        d = ins[ins.platform == plat]
        ax.scatter(d.launch_year, d.gb_per_day, s=40, color=C[k], marker=mk, label=plat, zorder=3, edgecolor="white")
    other = ins[~ins.platform.isin(markers)]
    ax.scatter(other.launch_year, other.gb_per_day, s=40, color=MUTED, marker="D", label="454 / Ion Torrent", zorder=3)
    nudge = {"Genome Analyzer IIx": (4, 6), "Ion PGM": (-8, -14), "PromethION 48": (-30, 8), "NovaSeq 6000": (5, -11),
             "HiSeq X": (4, -11), "HiSeq 4000": (4, 5)}
    for r in ins.itertuples():
        name = r.instrument.split(" (")[0]
        ax.annotate(name, (r.launch_year, r.gb_per_day), xytext=nudge.get(name, (4, -10)),
                    textcoords="offset points", fontsize=6.5, color=INK)
    il = ins[ins.platform == "Illumina"]
    b, c = np.polyfit(il.launch_year, np.log2(il.gb_per_day), 1)
    xx = np.array([2006, 2024])
    ax.plot(xx, 2 ** (c + b * xx), color=C[0], lw=1, ls=":", label=f"Illumina trend: x2 every {1 / b:.1f} years")
    ax.set_yscale("log")
    ax.set_xlabel("launch year")
    ax.set_ylabel("maximum output per instrument (Gb/day, log)")
    ax.set_title("Instrument capacity by generation")
    ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    p = runs.pivot(index="year", columns="platform", values="runs")
    p = p[p["ALL"] > 0]
    share = p.div(p["ALL"], axis=0) * 100
    groups = {"Illumina": ["ILLUMINA"], "Nanopore": ["OXFORD_NANOPORE"], "PacBio": ["PACBIO_SMRT"],
              "other": [c for c in share.columns if c not in ("ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT", "ALL")]}
    for k, (lab, cols) in enumerate(groups.items()):
        ax.plot(share.index, share[cols].sum(axis=1), lw=2, color=C[k] if lab != "other" else MUTED, label=lab, marker="o", ms=3)
    ax.set_ylabel("share of public read runs (%)")
    ax.set_xlabel("year the runs became public (ENA)")
    ax.set_title("What the community actually runs")
    ax.legend(fontsize=7.5, loc="center right")
    peak = p["ALL"].idxmax()
    ax.text(0.02, 0.55, f"public runs per year: {p['ALL'].iloc[0] / 1e3:.0f}k ({p.index[0]})\n"
                        f"→ {p['ALL'].max() / 1e6:.1f}M ({peak}) → {p['ALL'].iloc[-1] / 1e6:.1f}M ({p.index[-1]})",
            transform=ax.transAxes, fontsize=7.5, color=MUTED)
    save(fig, out, "fig_a2_instruments_platforms.png")


def fig_prices(t: Path, out: Path):
    seq = pd.read_csv(t / "sequencing_cost.csv")
    sto = pd.read_csv(t / "storage_cost.csv")
    fits = pd.read_csv(t / "price_fits.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax = axes[0]
    ax.plot(seq.year, seq.usd_per_mb, color=C[0], lw=2, marker="o", ms=2.5, label="NHGRI cost per Mb")
    for k, r in enumerate(fits[fits.fit != "storage"].itertuples()):
        y0, y1 = map(int, r.years.split("-"))
        d = seq[(seq.year >= y0) & (seq.year <= y1 + 1)]
        b, a = np.polyfit(d.year, np.log10(d.usd_per_mb), 1)
        xx = np.linspace(y0, 2045, 50)
        ax.plot(xx, 10 ** (a + b * xx), ls="--", lw=1.5, color=C[k + 1],
                label=f"{r.fit} trend {y0}-{y1}: {r.annual_change_pct:.0f}%/yr")
    ax.set_yscale("log")
    ax.set_xlabel("year")
    ax.set_ylabel("USD per raw megabase (log)")
    ax.set_title("Sequencing got cheaper fast, then slower")
    ax.legend(fontsize=7)
    ax = axes[1]
    dd = sto[sto.medium == "disk drives"]
    ax.plot(dd.year, dd.usd2020_per_tb, color=C[0], lw=2, marker="o", ms=2.5, label="hard disks, 2020 USD per TB")
    r = fits[fits.fit == "storage"].iloc[0]
    y0, y1 = map(int, r.years.split("-"))
    d = dd[dd.year.between(y0, y1)]
    b, a = np.polyfit(d.year, np.log10(d.usd2020_per_tb), 1)
    xx = np.linspace(y0, 2045, 50)
    ax.plot(xx, 10 ** (a + b * xx), ls="--", lw=1.5, color=C[1], label=f"trend {y0}-{y1}: {r.annual_change_pct:.0f}%/yr")
    ax.set_yscale("log")
    ax.set_xlim(1980, 2046)
    ax.set_xlabel("year")
    ax.set_ylabel("USD per TB (log)")
    ax.set_title("Disk storage keeps getting cheaper, slowly")
    ax.legend(fontsize=7)
    save(fig, out, "fig_a3_prices.png")


def fig_crossover(t: Path, out: Path):
    cur = pd.read_csv(t / "crossover_curves.tsv", sep="\t")
    summ = pd.read_csv(t / "crossover_summary.tsv", sep="\t")
    keys = ["current", "spring", "spring_bin2", "spring_bin2_tiered"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4))
    ax = axes[0]
    for k, key in enumerate(keys):
        d = cur[cur.key == key]
        if d.empty:
            continue
        ax.plot(d.year, d.ratio_median, lw=2, color=C[k], label=d.scenario.iloc[0])
        ax.fill_between(d.year, d.ratio_p05, d.ratio_p95, color=C[k], alpha=0.12, lw=0)
    ax.axhline(1, color=INK, lw=1)
    ax.text(2022.3, 1.35, "above: storing a new Mb for 10 years costs more than sequencing it", fontsize=7, color=INK)
    ax.set_yscale("log")
    ax.set_xlabel("year the data are sequenced")
    ax.set_ylabel("10-year storage cost / sequencing cost (log)")
    ax.set_title("Storage vs sequencing cost per megabase")
    ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    s = summ[summ.key.isin(keys)].copy()
    order = {k: i for i, k in enumerate(keys)}
    trends = ["both readings", "ngs_era", "recent"]
    for j, tr in enumerate(trends):
        d = s[s.sequencing_trend == tr].sort_values("key", key=lambda c: c.map(order))
        yy = np.array([order[k] for k in d.key]) + (j - 1) * 0.22
        med = d.median_year_all_draws.replace(np.inf, np.nan)
        lo, hi = d.p05_year, d.p95_year.fillna(2060)
        ax.hlines(yy, lo, hi, color=C[j] if tr != "both readings" else INK, lw=2)
        ax.scatter(med, yy, color=C[j] if tr != "both readings" else INK, s=30, zorder=3,
                   label={"both readings": "both sequencing trends (headline)", "ngs_era": "if 2008-2022 pace resumes",
                          "recent": "if 2015-2022 pace continues"}[tr])
    ax.set_yticks(range(len(keys)), [s[s.key == k].scenario.iloc[0] for k in keys], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(2021, 2061)
    ax.set_xlabel("crossover year (median, 5-95% range; bars reaching 2060 = later)")
    ax.set_title("When the curves cross, by storage practice")
    ax.legend(fontsize=7, loc="upper right")
    save(fig, out, "fig_a4_crossover.png")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--figures", type=Path, required=True)
    a = ap.parse_args(argv)
    a.figures.mkdir(parents=True, exist_ok=True)
    for f in (fig_growth, fig_instruments, fig_prices, fig_crossover):
        try:
            f(a.tables, a.figures)
        except Exception as exc:  # one broken figure must not hide the others
            print(f"skipped {f.__name__}: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
