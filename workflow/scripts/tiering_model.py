"""Tiering policy: predict which sequencing studies will be reused, price the policy (task 7).

Label. A study is 'reused' when >=2 papers mention its accession within a fixed window after the
data went public (default: the year before release to 4 years after). One paper is usually the
submitters' own; a second one is somebody coming back for the data. The fixed window matters:
an older study has simply had more years to be cited, and counting all papers would teach the
model 'old = reused'.

Validation is temporal, like the real decision: fit on studies public <= train_until_year,
evaluate on later ones.

Cost model (per study, horizon H months, prices re-checked against the AWS price list CSV):
  storage   GB x price(tier) x months in that tier
  retrieval each reuse event after the first paper = one full read of the study from its tier
  latency   reuse events that hit Deep Archive wait hours for a restore -> reported as 'delayed'
Tiers: hot = S3 Standard; warm = Glacier Instant Retrieval (millisecond access, 5.75x cheaper
storage, pays per GB read); cold = Deep Archive (23x cheaper, restore takes hours).
Policies: all-hot, all-warm, all-cold, age rule (hot 24 months, then cold), model (cold when the
predicted reuse probability < threshold, warm otherwise; threshold chosen on the TRAINING years
as the cheapest that delays <= 10% of reuse), and an oracle that knows the future.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

WINDOW_BEFORE, WINDOW_AFTER = 1, 4
PRICE_ROWS = {  # config tier -> (storage description, retrieval description) in the AWS CSV
    "hot": ("per GB - first 50 TB / month of storage used", None),   # matched exactly, see below
    "warm": ("of storage used in Glacier Instant Retrieval",
             "flat fee for all bytes retrieved in Glacier Instant Retrieval"),
    "cold": ("IntelligentTieringDAAStorage", "Bulk retrieval fee from Glacier Deep Archive"),
}


def check_prices(csv_path: Path, tiers: dict) -> pd.DataFrame:
    """Every price in config must equal the number in the downloaded AWS price list."""
    pl = pd.read_csv(csv_path, skiprows=5, dtype=str)
    rows = []
    for tier, (sdesc, rdesc) in PRICE_ROWS.items():
        for kind, desc, key in (("storage", sdesc, "usd_gb_month"), ("retrieval", rdesc, "usd_gb_retrieval")):
            if desc is None:
                continue
            d = pl["PriceDescription"].fillna("")
            # S3 Standard's row is the only one that ENDS with this text (Intelligent-Tiering rows go on)
            hit = pl[d.str.endswith(desc)] if tier == "hot" else pl[d.str.contains(desc, regex=False)]
            if hit.empty:
                raise SystemExit(f"price row '{desc}' not found in {csv_path}")
            price = float(hit["PricePerUnit"].iloc[0])
            cfg = tiers[tier][key]
            if not math.isclose(price, cfg, rel_tol=1e-6):
                raise SystemExit(f"{tier} {kind}: config says {cfg}, AWS price list says {price}")
            rows.append({"tier": tier, "kind": kind, "usd": price, "aws_description": hit["PriceDescription"].iloc[0]})
    return pd.DataFrame(rows)


def load(studies: Path, cache: Path) -> pd.DataFrame:
    st = pd.read_parquet(studies)
    years = {}
    for line in cache.read_text().splitlines():
        r = json.loads(line)
        if r.get("hits") is not None:
            years[r["acc"]] = r["years"]
    y0 = st["first_public"].dt.year

    def window_papers(i):
        ys = max((years.get(st.at[i, "study_accession"], []), years.get(st.at[i, "secondary_study_accession"], [])), key=len)
        lo, hi = y0.at[i] - WINDOW_BEFORE, y0.at[i] + WINDOW_AFTER
        return [y for y in ys if lo <= y <= hi]

    st["paper_years_window"] = [window_papers(i) for i in st.index]
    st["papers_window"] = st["paper_years_window"].str.len()
    st["reused"] = (st["papers_window"] >= 2).astype(int)
    st["year"] = y0
    st["gb"] = st["bytes"] / 1e9
    return st.dropna(subset=["gb"]).query("gb > 0").reset_index(drop=True)


def top_k(s: pd.Series, k: int) -> pd.Series:
    keep = s.value_counts().index[:k]
    return s.where(s.isin(keep), "other").fillna("other")


def features(st: pd.DataFrame, train_mask) -> tuple[pd.DataFrame, list, list]:
    X = pd.DataFrame({
        "log10_gb": np.log10(st["gb"]),
        "log10_runs": np.log10(st["n_runs"]),
        "log10_samples": np.log10(st["n_samples"].clip(lower=1)),
        "log10_gb_per_run": np.log10(st["gb"] / st["n_runs"]),
        "share_illumina": st["share_illumina"], "share_ont": st["share_ont"],
        "share_pacbio": st["share_pacbio"], "share_paired": st["share_paired"],
        "year": st["year"],
    })
    cats = {"library_strategy": 8, "library_source": 5, "instrument_model": 15, "center_name": 25}
    for c, k in cats.items():   # category vocabularies are learned on training years only
        keep = st.loc[train_mask, c].value_counts().index[:k]
        X[c] = st[c].where(st[c].isin(keep), "other").fillna("other")
    num = [c for c in X.columns if c not in cats]
    return X, num, list(cats)


def make_models(num, cat):
    pre = ColumnTransformer([("num", StandardScaler(), num),
                             ("cat", OneHotEncoder(handle_unknown="ignore"), cat)])
    return {
        "logistic": Pipeline([("pre", pre), ("m", LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5))]),
        "gradient_boosting": Pipeline([("pre", ColumnTransformer(
            [("num", "passthrough", num), ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), cat)])),
            ("m", HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                                 class_weight="balanced", random_state=1))]),
        "size_only": Pipeline([("pre", ColumnTransformer([("num", StandardScaler(), ["log10_gb"])])),
                               ("m", LogisticRegression(class_weight="balanced"))]),
    }


def policy_cost(st: pd.DataFrame, tier_of: np.ndarray, tiers: dict, H: int, cold_after: int | None = None) -> dict:
    """tier_of[i] in {'hot','warm','cold'} for the whole horizon; or, with cold_after, every study
    sits in hot for that many months and then moves to cold (the age rule)."""
    store = retr = 0.0
    delayed = events = 0
    for i, row in enumerate(st.itertuples(index=False)):
        gb = row.gb
        ys = sorted(row.paper_years_window)[1:]            # reuse = papers after the first one
        months = [max(0, (y - row.year) * 12 + 6) for y in ys]   # placed mid-year
        events += len(months)
        if cold_after is not None:
            t = min(cold_after, H)
            store += gb * (tiers["hot"]["usd_gb_month"] * t + tiers["cold"]["usd_gb_month"] * (H - t))
            late = [m for m in months if m >= t]
            retr += gb * tiers["cold"]["usd_gb_retrieval"] * len(late)
            delayed += len(late)
            continue
        tr = tiers[tier_of[i]]
        store += gb * tr["usd_gb_month"] * max(H, tr["min_months"])
        retr += gb * tr["usd_gb_retrieval"] * len(months)
        if tier_of[i] == "cold":
            delayed += len(months)
    return {"storage_usd": store, "retrieval_usd": retr, "total_usd": store + retr,
            "reuse_events": events, "delayed_events": delayed,
            "delayed_share": delayed / events if events else float("nan")}


def assign(p: np.ndarray, thr: float) -> np.ndarray:
    return np.where(p < thr, "cold", "warm")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--prices", type=Path, required=True)
    ap.add_argument("--tiering-json", required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--max-delayed-share", type=float, default=0.10,
                    help="threshold choice: cheapest policy that delays at most this share of reuse events")
    a = ap.parse_args(argv)
    tc = json.loads(a.tiering_json)
    tiers, H, split = tc["tiers"], int(tc["horizon_months"]), int(tc["train_until_year"])
    a.outdir.mkdir(parents=True, exist_ok=True)
    prices = check_prices(a.prices, tiers)
    prices.to_csv(a.outdir / "tiering_prices_checked.tsv", sep="\t", index=False)

    st = load(a.studies, a.studies.parent / "epmc_cache.jsonl")
    train = (st["year"] <= split).to_numpy()
    test = ~train
    X, num, cat = features(st, train)
    y = st["reused"].to_numpy()

    # ---- who gets reused (descriptive, all labelled studies)
    groups = []
    for col in ("library_strategy", "year"):
        key = st[col].astype(str) if col == "year" else top_k(st[col].astype(str), 10)
        g = st.assign(k=key).groupby("k").agg(
            studies=("reused", "size"), reused_share=("reused", "mean"), cited_share=("papers_window", lambda s: (s >= 1).mean()),
            tb=("gb", lambda s: s.sum() / 1e3))
        g["byte_share"] = g["tb"] / g["tb"].sum()
        groups.append(g.reset_index().assign(grouping=col))
    st["platform"] = np.select([st["share_ont"] > 0.5, st["share_pacbio"] > 0.5, st["share_illumina"] > 0.5],
                               ["ONT", "PacBio", "Illumina"], "mixed/other")
    g = st.groupby("platform").agg(studies=("reused", "size"), reused_share=("reused", "mean"),
                                   cited_share=("papers_window", lambda s: (s >= 1).mean()), tb=("gb", lambda s: s.sum() / 1e3))
    g["byte_share"] = g["tb"] / g["tb"].sum()
    groups.append(g.reset_index().rename(columns={"platform": "k"}).assign(grouping="platform"))
    pd.concat(groups).rename(columns={"k": "group"}).round(4).to_csv(a.outdir / "tiering_reuse_by_group.tsv", sep="\t", index=False)

    # ---- models
    metrics, preds = [], {}
    for name, model in make_models(num, cat).items():
        model.fit(X[train], y[train])
        p_tr, p_te = model.predict_proba(X[train])[:, 1], model.predict_proba(X[test])[:, 1]
        preds[name] = (p_tr, p_te)
        metrics.append({"model": name, "train_studies": int(train.sum()), "test_studies": int(test.sum()),
                        "train_reuse_rate": y[train].mean(), "test_reuse_rate": y[test].mean(),
                        "test_roc_auc": roc_auc_score(y[test], p_te), "test_pr_auc": average_precision_score(y[test], p_te),
                        "test_brier": brier_score_loss(y[test], p_te)})
    rng = np.random.default_rng(1)
    rand = rng.random(test.sum())
    metrics.append({"model": "random", "train_studies": int(train.sum()), "test_studies": int(test.sum()),
                    "train_reuse_rate": y[train].mean(), "test_reuse_rate": y[test].mean(),
                    "test_roc_auc": roc_auc_score(y[test], rand), "test_pr_auc": average_precision_score(y[test], rand),
                    "test_brier": float("nan")})
    mdf = pd.DataFrame(metrics)
    mdf.round(4).to_csv(a.outdir / "tiering_model_metrics.tsv", sep="\t", index=False)
    best = mdf[mdf.model.isin(["logistic", "gradient_boosting"])].sort_values("test_roc_auc").model.iloc[-1]

    # logistic coefficients = interpretable feature effects
    lr = make_models(num, cat)["logistic"].fit(X[train], y[train])
    names = lr.named_steps["pre"].get_feature_names_out()
    coef = pd.DataFrame({"feature": names, "coef": lr.named_steps["m"].coef_[0]})
    coef["odds_ratio_per_unit"] = np.exp(coef["coef"])
    coef.reindex(coef["coef"].abs().sort_values(ascending=False).index).round(4).to_csv(
        a.outdir / "tiering_feature_effects.tsv", sep="\t", index=False)

    # ---- threshold chosen on training years, then priced on test years
    st_tr, st_te = st[train].reset_index(drop=True), st[test].reset_index(drop=True)
    p_tr, p_te = preds[best]
    grid = np.unique(np.quantile(p_tr, np.linspace(0, 1, 101)))
    curve, chosen = [], None
    for thr in grid:
        c = policy_cost(st_tr, assign(p_tr, thr), tiers, H)
        curve.append({"threshold": thr, **c})
    curve = pd.DataFrame(curve)
    ok = curve[curve["delayed_share"] <= a.max_delayed_share]
    chosen = float(ok.sort_values("total_usd").threshold.iloc[0]) if not ok.empty else 0.0

    n = len(st_te)
    rows = []
    for name, kw in [("all hot (S3 Standard)", dict(tier_of=np.full(n, "hot"))),
                     ("all warm (Glacier Instant Retrieval)", dict(tier_of=np.full(n, "warm"))),
                     ("all cold (Deep Archive)", dict(tier_of=np.full(n, "cold"))),
                     ("age rule: hot 24 months, then Deep Archive", dict(tier_of=None, cold_after=24)),
                     (f"model ({best}): Deep Archive if p<thr, else Instant Retrieval", dict(tier_of=assign(p_te, chosen))),
                     ("oracle: Deep Archive if never reused, else Instant Retrieval",
                      dict(tier_of=np.where(st_te["papers_window"].to_numpy() < 2, "cold", "warm")))]:
        rows.append({"policy": name, **policy_cost(st_te, H=H, tiers=tiers, **kw)})
    costs = pd.DataFrame(rows)
    hot_cost = costs.loc[0, "total_usd"]
    costs["savings_vs_hot_pct"] = 100 * (1 - costs["total_usd"] / hot_cost)
    costs["test_tb"] = st_te["gb"].sum() / 1e3
    costs["horizon_months"] = H
    costs.round(4).to_csv(a.outdir / "tiering_costs.tsv", sep="\t", index=False)

    # sweep on the TEST years too (for the savings vs delay trade-off figure)
    sweep = []
    for thr in np.unique(np.quantile(p_te, np.linspace(0, 1, 51))):
        c = policy_cost(st_te, assign(p_te, thr), tiers, H)
        sweep.append({"threshold": thr, "cold_share_bytes": st_te.loc[p_te < thr, "gb"].sum() / st_te["gb"].sum(),
                      "savings_vs_hot_pct": 100 * (1 - c["total_usd"] / hot_cost), "delayed_share": c["delayed_share"]})
    pd.DataFrame(sweep).round(5).to_csv(a.outdir / "tiering_tradeoff.tsv", sep="\t", index=False)
    pd.DataFrame({"study": st_te["study_accession"], "year": st_te["year"], "gb": st_te["gb"].round(3),
                  "papers_window": st_te["papers_window"], "reused": st_te["reused"],
                  "p_reuse": np.round(p_te, 5)}).to_csv(a.outdir / "tiering_test_predictions.tsv", sep="\t", index=False)

    costs["savings_vs_warm_pct"] = 100 * (1 - costs["total_usd"] / costs.loc[1, "total_usd"])
    costs.round(4).to_csv(a.outdir / "tiering_costs.tsv", sep="	", index=False)
    m = costs.iloc[4]
    pd.DataFrame([{
        "policy": m["policy"], "storage_cost_saving_fraction": round(m["savings_vs_hot_pct"] / 100, 4),
        "delayed_share_of_reuse": round(m["delayed_share"], 4), "threshold": round(chosen, 5),
        "basis": f"E. coli SRA/ENA studies public {int(st_te['year'].min())}-{int(st_te['year'].max())}, "
                 f"{H}-month horizon, AWS us-east-1 list prices",
        "how_to_use": "multiply projected archive storage cost by (1 - saving fraction) in the crossover model",
    }]).to_csv(a.outdir / "tiering_savings_for_forecast.tsv", sep="\t", index=False)
    print(mdf.round(3).to_string(index=False))
    print(costs.round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
