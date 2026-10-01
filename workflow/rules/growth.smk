# Student A: archive growth (task 2), instruments (task 3), prices and the crossover projection (task 8).
# Inputs are downloaded by rules/data.smk (fetch_input) and recorded in downloads_*.tsv.

import json

GC = config["growth"]
XC = config["crossover"]


rule genbank_growth:
    """Downloads gbrel.txt (or parses a local copy given as config growth.gbrel_local)."""
    output: touch(config.get("results_dir", "results") + "/tables/.genbank_growth.done")
    params:
        inp=lambda w: f"--input {config['growth']['gbrel_local']}" if config["growth"].get("gbrel_local") else "",
        raw=f"{config['data_dir']}/raw",
        out=config.get("results_dir", "results") + "/tables",
        acc=config["growth"].get("accessions", "config/accessions.tsv"),
    shell:
        "python workflow/scripts/fetch_growth_stats.py {params.inp} --raw-dir {params.raw} "
        "--out-dir {params.out} --accessions {params.acc}"


rule sra_growth:
    input: f"{D}/raw/sra_stats.csv"
    output: csv=f"{RES}/tables/sra_growth.csv", issues=f"{RES}/tables/sra_growth_issues.tsv"
    shell: "python workflow/scripts/fetch_sra_growth.py --input {input} --out {output.csv} --issues {output.issues}"


rule prices:
    input: nhgri=f"{D}/raw/nhgri_costs.xls", owid=f"{D}/raw/owid_storage.csv", aws=f"{D}/raw/aws_s3_prices.csv"
    output:
        seq=f"{RES}/tables/sequencing_cost.csv", sto=f"{RES}/tables/storage_cost.csv",
        anchor=f"{RES}/tables/storage_anchor.tsv",
    shell:
        "python workflow/scripts/fetch_prices.py --nhgri {input.nhgri} --owid {input.owid} --aws {input.aws} "
        "--outdir {RES}/tables"


rule instruments:
    """Vendor capacity table + ENA run counts per platform and year (network; ~200 count queries)."""
    input: "config/instruments.tsv"
    output:
        ins=f"{RES}/tables/instruments.tsv", fit=f"{RES}/tables/instrument_capacity_fit.tsv",
        runs=f"{RES}/tables/platform_runs_by_year.tsv",
    params: cached=lambda w: f"--runs-table {GC['runs_table_local']}" if GC.get("runs_table_local") else ""
    benchmark: f"{RES}/benchmarks/instruments.tsv"
    shell:
        "python workflow/scripts/fetch_instrument_specs.py --instruments {input} --outdir {RES}/tables {params.cached}"


rule growth_models:
    input:
        flag=f"{RES}/tables/.genbank_growth.done", sra=f"{RES}/tables/sra_growth.csv",
        ins=f"{RES}/tables/instruments.tsv",
    output:
        fits=f"{RES}/tables/growth_fits.tsv", proj=f"{RES}/tables/growth_projections.tsv",
        back=f"{RES}/tables/growth_backtest.tsv", intake=f"{RES}/tables/sra_intake_vs_capacity.tsv",
    params: ff=json.dumps(GC["fit_from"])
    benchmark: f"{RES}/benchmarks/growth_models.tsv"
    shell:
        "python workflow/scripts/fit_growth_models.py --genbank {RES}/tables --sra {input.sra} --instruments {input.ins} "
        "--end-year {GC[forecast_end_year]} --backtest-years {GC[backtest_years]} --bootstrap {GC[bootstrap]} "
        "--fit-from {params.ff:q} "
        "--outdir {RES}/tables"


rule crossover:
    """Needs part B's measured bits/base and tiering saving, so it runs after those."""
    input:
        seq=f"{RES}/tables/sequencing_cost.csv", sto=f"{RES}/tables/storage_cost.csv",
        anchor=f"{RES}/tables/storage_anchor.tsv", sra=f"{RES}/tables/sra_growth.csv",
        comp=f"{RES}/tables/compress_illumina.tsv",
        binned=f"{RES}/tables/binning/compress_illumina_bin2.tsv",
        tier=f"{RES}/tables/tiering_savings_for_forecast.tsv",
    output:
        f"{RES}/tables/crossover_summary.tsv", f"{RES}/tables/crossover_curves.tsv",
        f"{RES}/tables/crossover_sensitivity.tsv", f"{RES}/tables/price_fits.tsv",
    params: cfg=json.dumps(XC)
    benchmark: f"{RES}/benchmarks/crossover.tsv"
    shell:
        "python workflow/scripts/project_crossover.py --tables {RES}/tables --config-json {params.cfg:q} "
        "--outdir {RES}/tables"


rule student_a_figures:
    input:
        rules.growth_models.output, rules.crossover.output, rules.instruments.output, rules.prices.output,
    output: touch(f"{RES}/figures/.done_student_a")
    shell: "python workflow/scripts/plots.py --tables {RES}/tables --figures {RES}/figures"
