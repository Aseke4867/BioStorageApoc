"""Figures for report sections 4-7 (student B). Reads results/tables, writes results/figures/*.png.

Every figure has axis labels with units; captions live in report/report.md. A figure is skipped
(with a message) when its table is missing, so a partial run still produces what it can.
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

# validated categorical palette (dataviz reference instance, light mode), fixed order
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
})


def tsv(d: Path, name: str):
    p = d / name
    return pd.read_csv(p, sep="\t") if p.exists() else None


def save(fig, out: Path, name: str):
    fig.tight_layout()
    fig.savefig(out / name, facecolor="white")
    plt.close(fig)
    print("wrote", name)


def fig_compression(t: Path, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    for ax, plat, title in zip(axes, ["illumina", "ont"], ["Illumina SRR25629153 (2x251)", "Nanopore MinION SRR25637822"]):
        df = tsv(t, f"compress_{plat}.tsv")
        if df is None:
            continue
        ax.scatter(df["comp_MBps"], df["bits_per_base"], s=40, color=C[0], zorder=3, edgecolor="white", linewidth=1)
        for k, (_, r) in enumerate(df.sort_values("bits_per_base").iterrows()):
            dy = 4 if k % 2 == 0 else -10          # alternate above/below so neighbours do not collide
            ax.annotate(r["codec"], (r["comp_MBps"], r["bits_per_base"]), xytext=(5, dy),
                        textcoords="offset points", fontsize=7.5, color=INK)
        ax.set_xscale("log")
        ax.set_xlabel("compression speed (MB of FASTQ per second, log scale)")
        ax.set_ylabel("archive size (bits per sequenced base)")
        ax.set_title(title)
        ax.set_ylim(0, None)
    save(fig, out, "fig_b1_compression.png")


def fig_streams(t: Path, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"width_ratios": [1, 1.4]})
    ax = axes[0]
    labels, base = [], None
    for i, (plat, name) in enumerate([("illumina", "Illumina"), ("ont", "Nanopore")]):
        df = tsv(t, f"qual_profile_{plat}.tsv")
        if df is None:
            continue
        s = df[df.kind == "stream"].set_index("key")["value_b"].astype(float)
        tot = s.sum()
        left = 0.0
        for j, k in enumerate(["names", "bases", "qualities"]):
            w = 100 * s[k] / tot
            ax.barh(i, w, left=left, color=C[j], edgecolor="white", linewidth=2, label=k if i == 0 else None)
            if w > 8:
                ax.text(left + w / 2, i, f"{w:.0f}%", ha="center", va="center", color="white", fontsize=8)
            left += w
        labels.append(name)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("share of compressed FASTQ (%, each stream zstd -19)")
    ax.set_title("What the archive is made of")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3)
    ax = axes[1]
    for i, (plat, name) in enumerate([("illumina", "Illumina"), ("ont", "Nanopore")]):
        df = tsv(t, f"qual_profile_{plat}.tsv")
        if df is None:
            continue
        q = df[df.kind == "qual"].copy()
        q["key"] = q["key"].astype(int)
        ax.bar(q["key"] + (i - 0.5) * 0.42, 100 * q["value_b"].astype(float), width=0.42, color=C[i], label=name)
    ax.set_xlabel("Phred base quality")
    ax.set_ylabel("share of bases (%)")
    ax.set_title("Distinct quality values in the raw data")
    ax.legend()
    save(fig, out, "fig_b2_streams_qualities.png")


def fig_binning(t: Path, out: Path):
    vb = tsv(t, "variant_benchmark.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    ax = axes[0]
    schemes = ["orig", "illumina8", "bin4", "bin2", "noqual"]
    for i, plat in enumerate(["illumina", "ont"]):
        vals = []
        for s in schemes:
            df = tsv(t / "binning", f"compress_{plat}_{s}.tsv")
            vals.append(np.nan if df is None else df.loc[df.codec.str.startswith("spring"), "bits_per_base"].iloc[0])
        x = np.arange(len(schemes))
        ax.bar(x + (i - 0.5) * 0.4, vals, width=0.4, color=C[i], label="Illumina" if i == 0 else "Nanopore")
    ax.set_xticks(np.arange(len(schemes)), schemes)
    ax.set_ylabel("Spring archive (bits per base)")
    ax.set_xlabel("quality binning scheme")
    ax.set_title("Archive size after binning")
    ax.legend()
    if vb is not None:
        vb = vb[vb.min_qual == 20]
        for ax, plat, title in [(axes[1], "illumina", "Illumina: SNP F1 vs depth"),
                                (axes[2], "ont", "Nanopore: SNP F1 vs depth")]:
            d = vb[(vb.type == "SNP") & vb.callset.str.startswith(plat)].copy()
            d["scheme"] = d.callset.str.split("_").str[1]
            d["depth"] = d.callset.str.split("_").str[2].str.replace("x", "")
            order = sorted(d.depth.unique(), key=lambda v: 999 if v == "full" else int(v))
            for k, s in enumerate([s for s in schemes if s in set(d.scheme)]):
                e = d[d.scheme == s].set_index("depth").reindex(order)
                ax.plot(range(len(order)), e["f1"], marker="o", ms=5, lw=2, color=C[k], label=s)
            ax.set_xticks(range(len(order)), [f"{o}x" if o != "full" else "full" for o in order])
            ax.set_xlabel("sequencing depth")
            ax.set_ylabel("SNP F1 vs assembly truth (QUAL>=20)")
            ax.set_ylim(0.90, 1.0)    # same scale on both panels: differences are honestly tiny
            ax.set_title(title)
            ax.legend(fontsize=7)
    save(fig, out, "fig_b3_binning_cost.png")


def fig_errors(t: Path, out: Path, data_dir: Path | None):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    ax = axes[0]
    rows = []
    for plat, name in [("illumina", "Illumina"), ("ont", "Nanopore (MinION)"),
                       ("signal_hac", "Nanopore R10.4.1 hac"), ("signal_sup", "Nanopore R10.4.1 sup")]:
        df = tsv(t, f"error_profile_{plat}.tsv")
        if df is not None:
            r = df.iloc[0]
            rows.append((name, r["mismatch_per_100"], r["ins_per_100"], r["del_per_100"]))
    for i, (name, mm, ins, de) in enumerate(rows):
        left = 0
        for j, (v, lab) in enumerate([(mm, "mismatch"), (ins, "insertion"), (de, "deletion")]):
            ax.barh(i, v, left=left, color=C[j], edgecolor="white", linewidth=2, label=lab if i == 0 else None)
            left += v
        ax.text(left, i, f" {left:.2f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel("errors per 100 aligned bases")
    ax.set_title("Error profile by platform (raw reads)")
    ax.set_xlim(0, max(r[1] + r[2] + r[3] for r in rows) * 1.25)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
    ax = axes[1]
    hp = []
    for plat, name in [("illumina", "Illumina"), ("ont", "Nanopore (MinION)")]:
        df = tsv(t, f"error_profile_{plat}.tsv")
        if df is not None:
            hp.append((name, 100 * df.iloc[0]["hp_share_of_indel_events"], 100 * df.iloc[0]["hp_share_of_genome"]))
    if hp:
        x = np.arange(len(hp))
        ax.bar(x, [h[1] for h in hp], width=0.5, color=C[0], label="indel events in homopolymers")
        ax.scatter(x, [h[2] for h in hp], color=INK, marker="_", s=600, linewidth=2, zorder=3,
                   label="share of genome in homopolymers (>=4)")
        ax.set_xticks(x, [h[0] for h in hp])
        ax.set_ylabel("%")
        ax.set_title("Indels concentrate in homopolymers")
        ax.legend(fontsize=7, loc="upper left")
    save(fig, out, "fig_b4_error_profile.png")


def fig_search(t: Path, out: Path):
    q, s = tsv(t, "search_latency.tsv"), tsv(t, "search_sizes.tsv")
    if q is None or s is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    g = q[q["query"] == "gene"].merge(s, on="format")
    ax = axes[0]
    ax.scatter(g["bits_per_base"], g["median_s"], s=45, color=C[0], edgecolor="white", zorder=3)
    for k, (_, r) in enumerate(g.sort_values("bits_per_base").iterrows()):
        dy = 4 if k % 2 == 0 else -10
        ax.annotate(f"{r['format']} (recall {r['recall_vs_bam']:.2f})", (r["bits_per_base"], r["median_s"]),
                    xytext=(5, dy), textcoords="offset points", fontsize=7)
    ax.set_yscale("log")
    ax.set_xlabel("storage (bits per sequenced base)")
    ax.set_ylabel("time to find reads of one gene (s, log)")
    ax.set_title("Size vs searchability (lacZ query)")
    ax = axes[1]
    sc = q[q["query"] == "full scan"]
    rr = q[q["query"].str.startswith("random")]
    ax.barh(range(len(sc)), sc["median_s"], color=C[1], label="read every record")
    ax.set_yticks(range(len(sc)), sc["format"])
    for i, v in enumerate(sc["median_s"]):
        ax.text(v, i, f" {v:.1f} s", va="center", fontsize=7.5)
    ax.set_xlabel("full-scan time (s)")
    if not rr.empty:
        txt = "\n".join(f"{r.format}: {1000 * r.median_s:.0f} ms" for r in rr.itertuples())
        ax.text(0.98, 0.04, "random 1 kb region (median)\n" + txt, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=7, color=MUTED)
    ax.set_title("Full decompression cost")
    save(fig, out, "fig_b5_search_latency.png")


def fig_signal(t: Path, out: Path, data_dir: Path | None):
    sv, arch = tsv(t, "signal_vs_basecalls.tsv"), tsv(t, "ont_archive_sizes.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    if arch is not None:
        ax = axes[0]
        a = arch.set_index("dataset")["bytes"] / 1e9
        pairs = [("plasmid", "plasmid_raw_pod5", ["plasmid_basecalls_hac", "plasmid_basecalls_sup"]),
                 ("pathogen", "pathogen_raw_pod5", ["pathogen_basecalls_fastq"])]
        labels, k = [], 0
        for name, raw, calls in pairs:
            if raw not in a:
                continue
            ax.bar(k, a[raw], color=C[0], width=0.6, label="raw signal (POD5)" if k == 0 else None)
            ax.text(k, a[raw], f"{a[raw]:.0f}", ha="center", va="bottom", fontsize=7.5)
            labels.append(f"{name}\nPOD5")
            k += 1
            for c in calls:
                if c in a:
                    ax.bar(k, a[c], color=C[1], width=0.6, label="basecalls (FASTQ)" if k == 1 else None)
                    ax.text(k, a[c], f"{a[c]:.0f}", ha="center", va="bottom", fontsize=7.5)
                    labels.append(c.split("_")[-1])
                    k += 1
        ax.set_xticks(range(len(labels)), labels, fontsize=7.5)
        ax.set_ylabel("GB in the public ONT archive")
        ax.set_title("Raw signal vs basecalls, whole datasets")
        ax.legend(fontsize=7, loc="center", bbox_to_anchor=(0.4, 0.62))
    if sv is not None and data_dir is not None:
        ax = axes[1]
        for i, m in enumerate(sv["model"]):
            p = data_dir / "signal" / f"{m}_reads.tsv"
            if p.exists():
                r = pd.read_csv(p, sep="\t")
                q = -10 * np.log10((1 - r["identity"]).clip(lower=1e-5))
                ax.hist(q, bins=np.arange(5, 50, 1), histtype="step", lw=2, color=C[i],
                        label=f"{m}: median Q{np.median(q):.1f}")
        ax.set_xlabel("per-read accuracy vs known plasmid (Phred Q)")
        ax.set_ylabel("reads")
        ax.set_title("Same signal, two basecalling models")
        ax.legend(fontsize=7.5)
    save(fig, out, "fig_b6_signal.png")


def fig_tiering(t: Path, out: Path):
    tr, costs, grp = tsv(t, "tiering_tradeoff.tsv"), tsv(t, "tiering_costs.tsv"), tsv(t, "tiering_reuse_by_group.tsv")
    preds = tsv(t, "tiering_test_predictions.tsv")
    if tr is None or costs is None:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    ax = axes[0]
    if grp is not None:
        y = grp[grp.grouping == "year"].copy()
        y["group"] = y["group"].astype(str)
        ax.plot(y["group"], 100 * y["reused_share"], marker="o", lw=2, color=C[0], label="reused (>=2 papers)")
        ax.plot(y["group"], 100 * y["cited_share"], marker="o", lw=2, color=C[1], label="cited (>=1 paper)")
        ax.set_xlabel("year the study became public")
        ax.set_ylabel("studies (%), papers within 4 years")
        ax.tick_params(axis="x", rotation=45)
        ax.set_title("How often E. coli studies come back")
        ax.legend(fontsize=7)
    ax = axes[1]
    if preds is not None:
        from sklearn.metrics import roc_curve
        f, tp, _ = roc_curve(preds["reused"], preds["p_reuse"])
        ax.plot(f, tp, lw=2, color=C[0], label="model, test years")
        ax.plot([0, 1], [0, 1], lw=1, color=MUTED, ls="--", label="random")
        ax.set_xlabel("false positive rate")
        ax.set_ylabel("true positive rate")
        ax.set_title("Predicting reuse at submission time")
        ax.legend(fontsize=7)
    ax = axes[2]
    ax.plot(100 * tr["delayed_share"], tr["savings_vs_hot_pct"], lw=2, color=C[0], label="model threshold sweep")
    for i, r in enumerate(costs.itertuples()):
        if r.policy.startswith("all hot"):
            continue
        ax.scatter(100 * r.delayed_share, r.savings_vs_hot_pct, s=45, color=C[(i % 7) + 1], zorder=3,
                   edgecolor="white", label=r.policy)
    ax.set_xlabel("reuse requests that wait for a restore (%)")
    ax.set_ylabel("storage+retrieval saving vs all-hot (%)")
    ax.set_title("Savings vs access delay (test years)")
    ax.legend(fontsize=6.5, loc="lower right")
    save(fig, out, "fig_b7_tiering.png")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--figures", type=Path, required=True)
    ap.add_argument("--data-dir", type=Path, default=None, help="for per-read tables kept out of git")
    a = ap.parse_args(argv)
    a.figures.mkdir(parents=True, exist_ok=True)
    for f in (fig_compression, fig_streams, fig_binning, fig_search, fig_tiering):
        try:
            f(a.tables, a.figures)
        except Exception as exc:  # one broken figure must not hide the others
            print(f"skipped {f.__name__}: {exc}", file=sys.stderr)
    for f in (fig_errors, fig_signal):
        try:
            f(a.tables, a.figures, a.data_dir)
        except Exception as exc:
            print(f"skipped {f.__name__}: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
