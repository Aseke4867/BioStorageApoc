# Contributions

Both partners worked with an AI coding agent (Claude Code); its use is disclosed in
`report/ai_disclosure.md`. Each partner is responsible for, has reviewed, and defends the
components listed under his name, and reviewed the other partner's components where they
connect.

## Student A (Asulkhan): block 1, tasks 1-3 and 8
- **Led:** project scope and data choices for block 1; the GenBank release-notes parser
  (`fetch_growth_stats.py`), written by hand: table location without fixed line numbers,
  date reconciliation between tables, header cross-checks.
- **Owns and defends:** SRA growth series (`fetch_sra_growth.py`), price series
  (`fetch_prices.py`), instrument table and platform mix (`config/instruments.tsv`,
  `fetch_instrument_specs.py`), growth models with backtests (`fit_growth_models.py`), the
  crossover projection (`project_crossover.py`), figures A1-A4 (`plots.py`), report summary,
  sections 1-3 and 8, and block 1 limitations. These were generated with the AI agent at his
  request; he checked the modelling decisions (fit windows, the three functional forms, the
  per-megabase crossover framing, the two NHGRI trend readings) and the numbers in the text
  against `results/tables/`.
- **Reviewed from block 2:** the bytes-per-base and tiering-saving values that feed the
  crossover scenarios.

## Student B (Ravil): blocks 2-3, tasks 4-7, and reproducibility
- **Led:** data acquisition with provenance (`fetch_data.py`, `rules/data.smk`); short- and
  long-read QC; compression benchmark (`bench_compress.py`); quality-score binning and its
  downstream cost against an assembly-derived truth set (`bin_qualities.py`,
  `compare_variants.py`); platform error profiles (`read_error_profile.py`); search latency on
  compressed formats (`bench_search.py`); raw-signal sizing and GPU re-basecalling
  (`pod5_storage.py`, `rules/signal.smk`); reuse signals and the tiering model
  (`fetch_reuse_signals.py`, `tiering_model.py`); figures B1-B7; report sections 4-7.
- **Reproducibility:** Snakemake workflow, pinned environment, test dataset and `make test`;
  ran the full pipeline for both blocks on his workstation (WSL2, 16 threads, RTX 4060 Ti).
- **Reviewed from block 1:** how the crossover model uses block 2 results, and the
  conclusion section.

## Defence
Asulkhan presents the summary, archive growth, platforms and the crossover projection. Ravil
presents compression, quality binning, long reads and raw signal, searchability, tiering and a
`make test` run. Both answer questions on any part.
