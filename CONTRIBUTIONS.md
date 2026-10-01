# Contributions

## Student A (Aseke)
- Block 1: growth, instruments, prices, forecast; report sections 1-3

## Student B (Ravil)
- Tasks 4-7 and reproducibility: data acquisition with provenance (`fetch_data.py`, `rules/data.smk`);
  short- and long-read QC; compression benchmark (`bench_compress.py`); quality-score binning and its
  downstream cost, measured against an assembly-derived truth set (`bin_qualities.py`,
  `compare_variants.py`); platform error profiles (`read_error_profile.py`); search latency on
  compressed formats (`bench_search.py`); raw-signal sizing and GPU re-basecalling (`pod5_storage.py`,
  `rules/signal.smk`); reuse-signal collection and the tiering cost model (`fetch_reuse_signals.py`,
  `tiering_model.py`); figures for sections 4-7 (`figures_b.py`); Snakemake workflow, pinned
  environment, test dataset and `make test`; report sections 4-7.
- Wiring of student A's `fetch_growth_stats.py` into the workflow (paths only; logic unchanged).
