# Run inside the conda env:  mamba env create -f environment.yml && conda activate storage-apocalypse
# Large data lives outside the repository (keep it on a Linux filesystem under WSL).
DATA  ?= $(HOME)/biostorage_data
CORES ?= 16

.PHONY: env repro test dry clean-test

env:
	mamba env create -f environment.yml

# Full analysis from raw public data (downloads ~3.5 GB; ~1.5 h on 16 cores + GPU for basecalling).
repro:
	snakemake -s workflow/Snakefile --cores $(CORES) --resources gpu=1 --config data_dir=$(DATA) --keep-going

# Same pipeline on the bundled real-data slices in test_data/ (no network, no GPU, a few minutes).
test:
	snakemake -s workflow/Snakefile --configfile config/config.yaml config/test.yaml --cores 4 --resources gpu=1 --rerun-incomplete
	@echo "test outputs: test_run/results/tables, test_run/results/figures"

dry:
	snakemake -s workflow/Snakefile -n --cores $(CORES) --config data_dir=$(DATA)

clean-test:
	rm -rf test_run /tmp/storage_apocalypse_test
