#!/usr/bin/env bash
# bootstrap.sh: create the repo skeleton (safe to re-run)
set -e

mkdir -p config workflow/rules workflow/scripts data test_data \
         results/tables results/figures results/benchmarks report

# --- Block 1 scripts (student A) ---
for f in fetch_growth_stats fetch_instrument_specs fetch_prices \
         fit_growth_models project_crossover; do
  [ -f workflow/scripts/$f.py ] || printf '"""TODO: %s (block 1, student A)"""\n' "$f" > workflow/scripts/$f.py
done
# --- Blocks 2/3 placeholders (student B owns these) ---
for f in bench_compress bin_qualities compare_variants \
         fetch_reuse_signals tiering_model plots; do
  [ -f workflow/scripts/$f.py ] || printf '"""TODO: %s (student B)"""\n' "$f" > workflow/scripts/$f.py
done

touch workflow/rules/{growth,compress,align_check,tiering}.smk
touch data/.gitkeep test_data/.gitkeep

[ -f config/accessions.tsv ] || printf 'accession\tarchive\tdescription\tsize_bytes\tdownload_date\tverified_by\n' > config/accessions.tsv

[ -f config/config.yaml ] || cat > config/config.yaml <<'EOF'
# All paths, parameters, prices live here. Every price needs a source URL + retrieval date.
paths:
  data: data
  results: results
growth:
  forecast_end_year: 2045      # assumption, change if needed
prices: {}                     # filled by fetch_prices.py, with source + date
EOF

[ -f workflow/Snakefile ] || cat > workflow/Snakefile <<'EOF'
configfile: "config/config.yaml"
include: "rules/growth.smk"
# include: "rules/compress.smk"     # student B
# include: "rules/align_check.smk"  # student B
# include: "rules/tiering.smk"      # student B

rule all:
    input: rules.growth_all.input if False else []   # replace when first rule exists
EOF

[ -f Makefile ] || cat > Makefile <<'EOF'
repro:
	snakemake -s workflow/Snakefile --cores 4 --use-conda
test:
	snakemake -s workflow/Snakefile --cores 2 --configfile config/config.yaml -n
EOF

[ -f .gitignore ] || cat > .gitignore <<'EOF'
data/*
!data/.gitkeep
.snakemake/
__pycache__/
.ipynb_checkpoints/
*.fastq* *.bam *.cram *.pod5 *.fast5
EOF

[ -f environment.yml ] || cat > environment.yml <<'EOF'
name: storage-apocalypse
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - snakemake
  - pandas
  - numpy
  - scipy
  - matplotlib
  - requests
  # pin exact versions after first successful run: `conda env export --no-builds`
EOF

[ -f CONTRIBUTIONS.md ] || printf '# Contributions\n\n## Student A (Aseke)\n- Block 1: growth, instruments, prices, forecast; report sections 1-3\n\n## Student B\n- Blocks 2-3: compression, alignment check, tiering; environment/reproducibility\n' > CONTRIBUTIONS.md
[ -f report/ai_disclosure.md ] || printf '# AI usage disclosure\n\n| Date | Tool | Used for | What we changed/verified |\n|---|---|---|---|\n' > report/ai_disclosure.md
[ -f report/report.md ] || printf '# The Storage Apocalypse\n\n## 1. Biological setting\n\n## 2. Data sources\n\n## 3. Platforms and throughput\n\n## Limitations\n' > report/report.md
[ -f README.md ] || printf '# BioStorageApoc\n\nProject 09: The Storage Apocalypse.\n\n## Quick start\n(TODO: up to 10 steps)\n' > README.md
