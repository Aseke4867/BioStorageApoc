# AI usage disclosure

| Date | Tool | Used for | What we changed/verified |
|---|---|---|---|
| 2026-10-01 | Claude (Anthropic), Claude Code agent | Part B (tasks 4-7): locating accessions, drafting scripts, Snakemake rules, figure code, and first drafts of report sections 4-7 | Ravil reviewed every script and rule, ran the pipeline end to end, checked that every number in sections 4-7 comes from a file in `results/tables/`, and checked all accessions resolve. Two problems the agent surfaced were verified by hand before use: ENA's FASTQ for SRR25629153 carries constant Q30 qualities (SRA Lite), and the run's flow-cell ID is a NovaSeq ID although the metadata says HiSeq 4000. |
