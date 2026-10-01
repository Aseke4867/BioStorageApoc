# Student B: raw Nanopore signal (task 5). Basecalling needs an NVIDIA GPU; set dorado.enabled
# false (the test config does) to run everything else without one.

DOR = config["dorado"]


rule pod5_inspect:
    input: f"{D}/raw/signal.pod5"
    output: stats=f"{RES}/tables/pod5_stats.tsv", ids=f"{D}/signal/read_ids.txt"
    params: n=DOR["max_reads"], seed=SEED
    benchmark: f"{RES}/benchmarks/pod5_inspect.tsv"
    shell:
        "python workflow/scripts/pod5_storage.py inspect --pod5 {input} --sample {params.n} "
        "--seed {params.seed} --ids {output.ids} --out {output.stats}"


rule pod5_recompress:
    """POD5 is already VBZ (zstd + streamvbyte) compressed; does a general codec find anything left?"""
    input: f"{D}/raw/signal.pod5"
    output: f"{RES}/tables/pod5_recompress.tsv"
    threads: T
    shell:
        """
        s0=$(date +%s.%N); b=$(zstd -19 -T{threads} -q -c {input} | wc -c); s1=$(date +%s.%N)
        printf "file_bytes\\tzstd19_bytes\\tratio\\tseconds\\n%s\\t%s\\t%s\\t%s\\n" $(stat -c %s {input}) $b \
          $(echo "scale=4; $(stat -c %s {input}) / $b" | bc) $(echo "$s1 - $s0" | bc) > {output}
        """


rule dorado_basecall:
    input: pod5=f"{D}/raw/signal.pod5", ids=f"{D}/signal/read_ids.txt"
    output: f"{D}/signal/{{model}}.fastq.gz"
    params: bin=DOR["bin"], dev=DOR["device"], models=f"{D}/models"
    resources: gpu=1             # one basecaller at a time on an 8 GB card
    threads: 4
    benchmark: f"{RES}/benchmarks/dorado_{{model}}.tsv"
    shell:
        "mkdir -p {params.models} && {params.bin} basecaller {wildcards.model} {input.pod5} -l {input.ids} "
        "-x {params.dev} --models-directory {params.models} --emit-fastq 2> {D}/signal/{wildcards.model}.log "
        "| pigz -p 4 > {output}"


rule signal_identity:
    """Plasmid sequences are known exactly, so identity = basecaller accuracy (no biology in it)."""
    input: fq=f"{D}/signal/{{model}}.fastq.gz", ref=f"{D}/raw/plasmid_refs.fa"
    output: per_read=f"{D}/signal/{{model}}_reads.tsv", summary=f"{RES}/tables/error_profile_signal_{{model}}.tsv"
    threads: T
    shell:
        "minimap2 -ax map-ont -t {threads} {input.ref} {input.fq} 2>/dev/null | samtools view -F 0x904 "
        "| python workflow/scripts/read_error_profile.py --ref {input.ref} --platform ont_{wildcards.model} "
        "--max-reads 50000 --per-read {output.per_read} --summary {output.summary}"


rule archive_sizes:
    output: f"{RES}/tables/ont_archive_sizes.tsv"
    params:
        b=config["plasmid_refs_s3"]["bucket"],
        p=" ".join(f"{k}={v}" for k, v in config["archive_size_prefixes"].items()),
    shell: "python workflow/scripts/pod5_storage.py archive --bucket {params.b} --prefix {params.p} --out {output}"


rule signal_summary:
    input:
        stats=f"{RES}/tables/pod5_stats.tsv",
        fq=expand(f"{D}/signal/{{m}}.fastq.gz", m=DOR["models"]),
        ident=expand(f"{RES}/tables/error_profile_signal_{{m}}.tsv", m=DOR["models"]),
        bench=expand(f"{RES}/benchmarks/dorado_{{m}}.tsv", m=DOR["models"]),
    output: f"{RES}/tables/signal_vs_basecalls.tsv"
    params:
        fq=" ".join(f"{m}={D}/signal/{m}.fastq.gz" for m in DOR["models"]),
        ident=" ".join(f"{m}={RES}/tables/error_profile_signal_{m}.tsv" for m in DOR["models"]),
        bench=" ".join(f"{m}={RES}/benchmarks/dorado_{m}.tsv" for m in DOR["models"]),
    shell:
        "python workflow/scripts/pod5_storage.py summarize --stats {input.stats} --fastq {params.fq} "
        "--identity {params.ident} --bench {params.bench} --out {output}"
