# Student B: QC (tasks 4, 5), compression benchmark (4, 5), quality binning (4).

PLAT_FILES = {  # QC'd reads per platform, by mate
    "illumina": ["R1", "R2"],
    "ont": ["SE"],
}


rule fastp_illumina:
    """Short-read QC. Thresholds come from config and are argued in the report (section 4)."""
    input: r1=f"{D}/raw/illumina_R1.fastq.gz", r2=f"{D}/raw/illumina_R2.fastq.gz"
    output:
        r1=f"{D}/qc/illumina_R1.fastq", r2=f"{D}/qc/illumina_R2.fastq",
        json=f"{RES}/qc/fastp_illumina.json", html=f"{RES}/qc/fastp_illumina.html",
    params: c=config["fastp"]
    threads: min(T, 16)
    benchmark: f"{RES}/benchmarks/fastp_illumina.tsv"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --dedup --dup_calc_accuracy 3 --overrepresentation_analysis "
        "--qualified_quality_phred {params.c[qualified_quality]} "
        "--unqualified_percent_limit {params.c[unqualified_percent]} "
        "--length_required {params.c[length_required]} "
        "--cut_tail --cut_tail_mean_quality {params.c[cut_tail_quality]} "
        "-w {threads} -j {output.json} -h {output.html} 2> /dev/null"


rule nanoplot_ont:
    input: f"{D}/raw/ont.fastq.gz"
    output: f"{RES}/qc/nanoplot_ont/NanoStats.txt"
    threads: T
    benchmark: f"{RES}/benchmarks/nanoplot_ont.tsv"
    shell: "NanoPlot --fastq {input} -t {threads} -o $(dirname {output}) --no_static --tsv_stats 2> /dev/null"


rule filter_ont:
    """Long-read QC: no adapter/duplicate logic from short reads; filter on length and mean Q."""
    input: f"{D}/raw/ont.fastq.gz"
    output: f"{D}/qc/ont_SE.fastq"
    params: c=config["ont_filter"]
    threads: 4
    shell: "seqkit seq -j {threads} -m {params.c[min_length]} -Q {params.c[min_mean_q]} {input} > {output}"


rule read_stats:
    input:
        raw=[f"{D}/raw/illumina_R1.fastq.gz", f"{D}/raw/illumina_R2.fastq.gz", f"{D}/raw/ont.fastq.gz"]
            + [INPUT_FILES[k] for k in ("illumina_lite_R1", "illumina_lite_R2") if k in INPUT_FILES],
        qc=[f"{D}/qc/illumina_R1.fastq", f"{D}/qc/illumina_R2.fastq", f"{D}/qc/ont_SE.fastq"],
    output: f"{RES}/tables/read_stats.tsv"
    threads: 6
    shell: "seqkit stats -a -T -j {threads} {input.raw} {input.qc} > {output}"


rule plain_raw:
    """Uncompressed copy of the raw archive files: the common starting point for every codec."""
    input: f"{D}/raw/{{name}}.fastq.gz"
    output: f"{D}/plain/{{name}}.fastq"
    wildcard_constraints: name="illumina_R1|illumina_R2|ont|illumina_lite_R1|illumina_lite_R2"
    shell: "pigz -dc {input} > {output}"


rule bench_compress_raw:
    input: lambda w: [f"{D}/plain/illumina_R1.fastq", f"{D}/plain/illumina_R2.fastq"]
                     if w.plat == "illumina" else [f"{D}/plain/ont.fastq"]
    output: f"{RES}/tables/compress_{{plat}}.tsv"
    params: long=lambda w: "--long-reads" if w.plat == "ont" else "", work=f"{D}/scratch"
    threads: workflow.cores   # runs alone, so timings are not disturbed by other jobs
    benchmark: f"{RES}/benchmarks/bench_compress_{{plat}}.tsv"
    shell:
        "mkdir -p {params.work} && python workflow/scripts/bench_compress.py --fastq {input} "
        "--dataset {wildcards.plat} {params.long} --threads {T} --workdir {params.work} --out {output}"


rule qual_profile:
    input: lambda w: f"{D}/qc/{w.plat}_{PLAT_FILES[w.plat][0]}.fastq"
    output: f"{RES}/tables/qual_profile_{{plat}}.tsv"
    threads: T
    shell: "python workflow/scripts/bin_qualities.py --in {input} --profile {output} --threads {threads}"


rule bin_qualities:
    input: f"{D}/qc/{{plat}}_{{mate}}.fastq"
    output: f"{D}/binned/{{plat}}_{{scheme}}_{{mate}}.fastq"
    wildcard_constraints: scheme="|".join(config["binning_schemes"])
    shell: "python workflow/scripts/bin_qualities.py --in {input} --scheme {wildcards.scheme} --out {output}"


rule bench_binned:
    """What binning buys: archive size of each binned file with the two best lossless codecs."""
    input: lambda w: [f"{D}/binned/{w.plat}_{w.scheme}_{m}.fastq" for m in PLAT_FILES[w.plat]]
    output: f"{RES}/tables/binning/compress_{{plat}}_{{scheme}}.tsv"
    params: long=lambda w: "--long-reads" if w.plat == "ont" else "", work=f"{D}/scratch"
    threads: workflow.cores   # runs alone, so timings are not disturbed by other jobs
    shell:
        "mkdir -p {params.work} && python workflow/scripts/bench_compress.py --fastq {input} "
        "--dataset {wildcards.plat}_{wildcards.scheme} --codecs zstd-19,spring {params.long} "
        "--threads {T} --workdir {params.work} --out {output}"


rule bench_compress_lite:
    """The archive's own lossy step: ENA's SRA Lite FASTQ (all qualities = Q30) of the same run."""
    input: [f"{D}/plain/illumina_lite_R1.fastq", f"{D}/plain/illumina_lite_R2.fastq"]
    output: f"{RES}/tables/compress_illumina_lite.tsv"
    params: work=f"{D}/scratch"
    threads: workflow.cores   # runs alone, so timings are not disturbed by other jobs
    shell:
        "mkdir -p {params.work} && python workflow/scripts/bench_compress.py --fastq {input} "
        "--dataset illumina_lite --codecs gzip-6,zstd-19,spring --threads {T} --workdir {params.work} --out {output}"
