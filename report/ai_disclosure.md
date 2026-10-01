# AI usage disclosure

Tool: Claude (Anthropic) through Claude Code, an agent that runs commands and edits files.

| Part | Used for | What we did ourselves / verified |
|---|---|---|
| Block 1 (tasks 1-3, 8) |SRA growth, instrument table (transcribed from vendor sheets), NHGRI and OWID price series, growth models, crossover projection, plotsn | The GenBank parser (`fetch_growth_stats.py`) was written by Aseke. Aseke reviewed the generated block 1 code and its modelling decisions: fitting from 2012/2014 because the archives' start-up years dominate an exponential fit, replacing a quadratic-in-log form that turned downwards when extrapolated, and comparing costs per sequenced megabase rather than per archive. He also checked the report's numbers against `results/tables/`. Ravil ran the block 1 pipeline and reviewed how it uses block 2 results. |
| Block 2-3 (tasks 4-7), workflow | Locating accessions; drafting scripts, Snakemake rules and figure code; | Ravil chose the scope (all of tasks 4-7, real data only), reviewed every script and rule, ran the pipeline end to end, and checked that each number in sections 4-7 comes from a file in `results/tables/` and that every accession resolves. Two problems found during the run were checked by hand before use: ENA's FASTQ for SRR25629153 has constant Q30 qualities (SRA Lite), and the run's flow-cell ID is a NovaSeq ID although the metadata says HiSeq 4000. |
