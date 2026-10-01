# Contributions

## Student A (Aseke)
- GenBank release-notes fetcher and parser (`workflow/scripts/fetch_growth_stats.py`): table
  location without hard-coded line numbers, date reconciliation across tables, header
  cross-checks, and the parsed release 273.0 table.
- Block 1 scope and data choices (growth, instruments, prices, forecast). The remaining block 1
  analysis (SRA growth, instrument table, price series, growth models, crossover projection,
  report sections 1-3 and 8) was produced with the AI agent in student B's working session, at
  student A's request and with his approval (see `report/ai_disclosure.md`).
  <!-- Aseke: edit this paragraph so that it describes what you reviewed and can defend. -->

## Student B (Ravil)
- Tasks 4-7 and reproducibility: data acquisition with provenance (`fetch_data.py`, `rules/data.smk`);
  short- and long-read QC; compression benchmark (`bench_compress.py`); quality-score binning and its
  downstream cost, measured against an assembly-derived truth set (`bin_qualities.py`,
  `compare_variants.py`); platform error profiles (`read_error_profile.py`); search latency on
  compressed formats (`bench_search.py`); raw-signal sizing and GPU re-basecalling (`pod5_storage.py`,
  `rules/signal.smk`); reuse-signal collection and the tiering cost model (`fetch_reuse_signals.py`,
  `tiering_model.py`); figures for sections 4-7 (`figures_b.py`); Snakemake workflow, pinned
  environment, test dataset and `make test`; report sections 4-7.
- Ran the block 1 pipeline (growth, prices, instruments, crossover) in the same environment.
