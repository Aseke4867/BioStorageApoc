# Student B: tiering policy validated against real reuse signals (task 7).

TC = config["tiering"]


rule reuse_signals:
    """ENA run metadata (~600k runs) + Europe PMC mentions per study. Network-heavy, resumable."""
    output:
        studies=f"{D}/tiering/studies.parquet", runs=f"{D}/tiering/ena_runs.parquet",
        issues=f"{RES}/tables/tiering_metadata_issues.tsv",
    params:
        q=TC["query"], a=TC["public_from"], b=TC["public_to"],
        extra=TC.get("fetch_extra", ""),
    threads: 10
    benchmark: f"{RES}/benchmarks/reuse_signals.tsv"
    shell:
        "python workflow/scripts/fetch_reuse_signals.py --query {params.q:q} --out $(dirname {output.studies}) "
        "--public-from {params.a} --public-to {params.b} --threads {threads} {params.extra} "
        "&& cp $(dirname {output.studies})/metadata_issues.tsv {output.issues}"


rule tiering_model:
    input:
        studies=f"{D}/tiering/studies.parquet", runs=f"{D}/tiering/ena_runs.parquet",
        prices=f"{D}/raw/aws_s3_prices.csv",
    output:
        metrics=f"{RES}/tables/tiering_model_metrics.tsv",
        costs=f"{RES}/tables/tiering_costs.tsv",
        features=f"{RES}/tables/tiering_feature_effects.tsv",
        reuse=f"{RES}/tables/tiering_reuse_by_group.tsv",
        savings=f"{RES}/tables/tiering_savings_for_forecast.tsv",
        tradeoff=f"{RES}/tables/tiering_tradeoff.tsv",
        prices=f"{RES}/tables/tiering_prices_checked.tsv",
        preds=f"{RES}/tables/tiering_test_predictions.tsv",
    params: tc=json.dumps(TC)
    benchmark: f"{RES}/benchmarks/tiering_model.tsv"
    shell:
        "python workflow/scripts/tiering_model.py --studies {input.studies} --runs {input.runs} "
        "--prices {input.prices} --tiering-json {params.tc:q} --outdir {RES}/tables"
