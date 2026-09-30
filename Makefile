repro:
	snakemake -s workflow/Snakefile --cores 4 --use-conda
test:
	snakemake -s workflow/Snakefile --cores 2 --configfile config/config.yaml -n
