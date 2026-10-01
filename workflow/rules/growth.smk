# Student A: archive growth. (Wiring added by student B so the existing fetcher runs from the
# repository root with correct paths; the analysis logic is student A's.)

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
