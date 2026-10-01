# BioStorageApoc — The Storage Apocalypse (AITU Bioinformatics, Project 09)

Sequencing output has grown faster than storage got cheaper. This repository models when the
curves cross (growth, prices, forecast — student A) and measures what compression, quality-score
binning, raw-signal policy and storage tiering can buy (tasks 4–7 — student B).

Everything runs from public data with one Snakemake workflow. The report is `report/report.md`.

## Quick start (under ten steps)

1. Linux or WSL2 (Ubuntu). Keep data on the Linux filesystem, not under `/mnt/c`.
2. Install Miniforge: <https://github.com/conda-forge/miniforge>
3. `mamba env create -f environment.yml && conda activate storage-apocalypse`
4. (Optional, NVIDIA GPU) `bash workflow/scripts/install_dorado.sh` — needed only for re-basecalling.
5. Check the pipeline on the bundled real-data slices (no network, no GPU, a few minutes):
   `make test` → `test_run/results/`
6. Full run: `make repro DATA=$HOME/biostorage_data CORES=16` (downloads ~3.5 GB of inputs,
   ~1.5 h on 16 cores; without a GPU set `dorado: {enabled: false}` in `config/config.yaml`).
7. Results: tables in `results/tables/`, figures in `results/figures/`, per-step runtime and
   memory in `results/benchmarks/` (Snakemake `benchmark:`), provenance of every download in
   `results/tables/downloads_student_b.tsv` and `config/accessions.tsv`.

The headline numbers of part B are produced by `make repro`: compression ratio per codec
(`compress_*.tsv`), cost of quality binning in variant-calling F1 (`variant_benchmark.tsv`),
signal-vs-basecalls bytes and re-basecalling accuracy gain (`signal_vs_basecalls.tsv`), query
latency per storage format (`search_latency.tsv`), and the tiering saving
(`tiering_costs.tsv`, `tiering_savings_for_forecast.tsv`).

## Data (all public; accessions and checksums in `config/config.yaml`)

| What | Accession / source | Used for |
|---|---|---|
| Illumina 2×251, E. coli K-12 MG1655 | SRR25629153 (SRA Normalized file from NCBI; ENA FASTQ = SRA Lite) | QC, compression, binning, alignment, search |
| Nanopore MinION, same BioSample | SRR25637822 (SAMN36977599, PRJNA1005239) | long-read QC, compression, error profile, variants |
| Raw signal POD5, R10.4.1 plasmids | ont-open-data `plasmid_2025.04` FBC24981 file 17 | signal size, re-basecalling hac vs sup |
| Plasmid reference sequences | ont-open-data `plasmid_2025.04/analysis/inputs/references` | basecalling accuracy |
| Reference genomes | NC_000913.3 (K-12 MG1655), NC_012967.1 (B REL606) | alignment; assembly-derived truth set |
| Run metadata, all E. coli runs | ENA Portal API `read_run`, `tax_tree(562)` (~607k runs) | tiering features |
| Literature mentions of each study | Europe PMC REST (`ACCESSION_ID` + free text) | reuse signal for tiering |
| Cloud storage prices | AWS Price List API, AmazonS3 us-east-1 | tiering cost model |
| GenBank release statistics | `ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt` | growth model (student A) |

## Layout

```
config/config.yaml        all inputs, parameters and prices (with sources)
config/test.yaml          overrides for the bundled test slices
workflow/Snakefile        entry point; rules/ holds one file per block
workflow/rules/           growth (A) | data, compress, align_check, signal, tiering (B)
workflow/scripts/         one script per analysis step, each with a docstring
test_data/                ~10 MB of real data slices for `make test`
results/                  small result tables and figures (committed)
report/                   report, AI-usage disclosure
```

## Test dataset

`test_data/` contains the first 20,000 Illumina read pairs (with original qualities), 150
Nanopore reads ≥1 kb, three POD5 reads, both reference genomes (gzipped), the plasmid
references, a snapshot of the AWS S3 price list, a 2,000-study slice of the tiering inputs with
its Europe PMC answers, and GenBank release notes 273.0. All are cut from the real inputs above.

## Contributions

See `CONTRIBUTIONS.md`; AI assistance is disclosed in `report/ai_disclosure.md`.
