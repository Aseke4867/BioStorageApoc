# Student B: alignment core (task 6), downstream cost of binning (task 4), short vs long (task 5),
# searchability of compressed archives (task 6).

import json

G = config["genome_size"]
SEED = config["seed"]


def qc_bases(plat):
    """Bases after QC, used to turn a target depth into a sampling fraction."""
    if plat == "illumina":
        with open(f"{RES}/qc/fastp_illumina.json") as fh:
            return json.load(fh)["summary"]["after_filtering"]["total_bases"]
    tot = 0
    with open(f"{RES}/tables/read_stats.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r["file"].endswith("ont_SE.fastq"):
                tot = int(r["sum_len"])
    return tot


rule subsample:
    """Same seed for both mates keeps pairs in step (seqkit samples by record index)."""
    input:
        fq=f"{D}/binned/{{plat}}_{{scheme}}_{{mate}}.fastq",
        qc=lambda w: f"{RES}/qc/fastp_illumina.json" if w.plat == "illumina" else f"{RES}/tables/read_stats.tsv",
    output: temp(f"{D}/depth/{{plat}}_{{scheme}}_{{depth}}x_{{mate}}.fastq")
    wildcard_constraints: depth=r"\d+|full", scheme="|".join(config["binning_schemes"])
    params:
        frac=lambda w: 1.0 if w.depth == "full" else min(1.0, int(w.depth) * G / qc_bases(w.plat))
    shell:
        "if [ {wildcards.depth} = full ]; then ln -sf \"$(realpath {input.fq})\" {output}; "
        "else seqkit sample -p {params.frac} -s {SEED} {input.fq} > {output} 2>/dev/null; fi"


rule align_illumina:
    input:
        r1=f"{D}/depth/illumina_{{scheme}}_{{depth}}x_R1.fastq",
        r2=f"{D}/depth/illumina_{{scheme}}_{{depth}}x_R2.fastq",
        ref=f"{D}/raw/ref_{{ref}}.fa", idx=f"{D}/raw/ref_{{ref}}.fa.bwt.2bit.64",
    output: bam=f"{D}/aln/{{ref}}/illumina_{{scheme}}_{{depth}}x.bam"
    wildcard_constraints: depth=r"\d+|full"
    threads: T
    benchmark: f"{RES}/benchmarks/align/illumina_{{ref}}_{{scheme}}_{{depth}}x.tsv"
    shell:
        "bwa-mem2 mem -t {threads} -R '@RG\\tID:ill\\tSM:K12\\tPL:ILLUMINA' {input.ref} {input.r1} {input.r2} 2>/dev/null"
        " | samtools sort -@4 -o {output.bam} - && samtools index {output.bam}"


rule align_ont:
    """map-ont preset: k=15 seeds and gap costs tuned for ~5-10% indel-dominated error."""
    input:
        fq=f"{D}/depth/ont_{{scheme}}_{{depth}}x_SE.fastq",
        ref=f"{D}/raw/ref_{{ref}}.fa",
    output: bam=f"{D}/aln/{{ref}}/ont_{{scheme}}_{{depth}}x.bam"
    wildcard_constraints: depth=r"\d+|full"
    threads: T
    benchmark: f"{RES}/benchmarks/align/ont_{{ref}}_{{scheme}}_{{depth}}x.tsv"
    shell:
        "minimap2 -ax map-ont -t {threads} -R '@RG\\tID:ont\\tSM:K12\\tPL:ONT' {input.ref} {input.fq} 2>/dev/null"
        " | samtools sort -@4 -o {output.bam} - && samtools index {output.bam}"


rule call_variants:
    """Haploid bacterium -> --ploidy 1. ONT uses bcftools' own ONT error-model preset."""
    input:
        bam=f"{D}/aln/REL606/{{plat}}_{{scheme}}_{{depth}}.bam", ref=f"{D}/raw/ref_REL606.fa",
    output: vcf=f"{D}/vcf/{{plat}}_{{scheme}}_{{depth}}.vcf.gz"
    wildcard_constraints: depth=r"(\d+|full)x"
    params: x=lambda w: "-X ont-sup-1.20" if w.plat == "ont" else "-q 20 -Q 13"
    threads: 2
    benchmark: f"{RES}/benchmarks/call/{{plat}}_{{scheme}}_{{depth}}.tsv"
    shell:
        "bcftools mpileup {params.x} -f {input.ref} -a AD,DP -d 10000 {input.bam} 2>/dev/null"
        " | bcftools call -mv --ploidy 1 2>/dev/null"
        " | bcftools norm -f {input.ref} -m -any -a 2>/dev/null | bcftools sort -Oz -o {output.vcf} 2>/dev/null"
        " && bcftools index -f {output.vcf}"


def callsets():
    out = []
    for plat, depths in (("illumina", config["depths_illumina"]), ("ont", config["depths_ont"])):
        for s in config["binning_schemes"]:
            if plat == "ont" and s == "illumina8":
                continue          # Illumina's bins are meaningless for a Nanopore quality model
            for d in depths:
                out.append(f"{plat}_{s}_{d}x")
    return out


rule variant_benchmark:
    input:
        truth=f"{D}/truth/truth.vcf.gz", bed=f"{RES}/tables/truth_confident.bed",
        calls=expand(f"{D}/vcf/{{cs}}.vcf.gz", cs=callsets()),
    output: f"{RES}/tables/variant_benchmark.tsv"
    params: names=" ".join(callsets())
    shell:
        "python workflow/scripts/compare_variants.py --truth {input.truth} --bed {input.bed} "
        "--calls {input.calls} --names {params.names} --out {output}"


# ---------------------------------------------------------------- raw reads vs the correct reference
RAW_ALIAS = {"illumina_raw_fullx_R1": "illumina_R1", "illumina_raw_fullx_R2": "illumina_R2",
             "ont_raw_fullx_SE": "ont"}


rule raw_link:
    """Raw (pre-QC) reads enter the alignment rules as scheme 'raw' at full depth."""
    input: lambda w: f"{D}/plain/{RAW_ALIAS[w.name]}.fastq"
    output: f"{D}/depth/{{name}}.fastq"
    wildcard_constraints: name="|".join(RAW_ALIAS)
    shell: "ln -sf \"$(realpath {input})\" {output}"


ruleorder: raw_link > subsample


rule error_profile:
    input:
        bam=f"{D}/aln/MG1655/{{plat}}_raw_fullx.bam", ref=f"{D}/raw/ref_MG1655.fa",
    output:
        per_read=f"{D}/profile/{{plat}}_reads.tsv", summary=f"{RES}/tables/error_profile_{{plat}}.tsv",
    params: n=lambda w: 200000 if w.plat == "illumina" else 20000
    shell:
        "samtools view -F 0x904 {input.bam} | python workflow/scripts/read_error_profile.py "
        "--ref {input.ref} --platform {wildcards.plat} --max-reads {params.n} "
        "--per-read {output.per_read} --summary {output.summary}"


rule mapping_stats:
    input: bam=f"{D}/aln/MG1655/{{plat}}_raw_fullx.bam"
    output: flag=f"{RES}/tables/flagstat_{{plat}}.tsv", cov=f"{RES}/tables/coverage_{{plat}}.tsv"
    shell: "samtools flagstat -O tsv {input.bam} > {output.flag} && samtools coverage {input.bam} > {output.cov}"


rule aligned_formats:
    """BAM -> CRAM 3.1 (default profile) and CRAM 3.1 'archive' profile (slower, smaller)."""
    input: bam=f"{D}/aln/MG1655/{{plat}}_raw_fullx.bam", ref=f"{D}/raw/ref_MG1655.fa"
    output:
        cram=f"{D}/formats/{{plat}}.cram", arch=f"{D}/formats/{{plat}}.archive.cram",
        tsv=f"{RES}/tables/aligned_formats_{{plat}}.tsv",
    threads: T
    shell:
        """
        s0=$(date +%s.%N); samtools view -@ {threads} -C -T {input.ref} --output-fmt-option version=3.1 -o {output.cram} {input.bam}; s1=$(date +%s.%N)
        samtools view -@ {threads} -T {input.ref} -O cram,version=3.1,archive -o {output.arch} {input.bam}; s2=$(date +%s.%N)
        samtools index {output.cram}; samtools index {output.arch}
        printf "format\\tbytes\\tencode_s\\n" > {output.tsv}
        printf "BAM\\t%s\\tNA\\n" $(stat -c %s {input.bam}) >> {output.tsv}
        printf "CRAM 3.1\\t%s\\t%s\\n" $(stat -c %s {output.cram}) $(echo "$s1 - $s0" | bc) >> {output.tsv}
        printf "CRAM 3.1 archive\\t%s\\t%s\\n" $(stat -c %s {output.arch}) $(echo "$s2 - $s1" | bc) >> {output.tsv}
        """


rule search_formats:
    """The non-aligned searchable forms of the raw Illumina run."""
    input:
        r1=f"{D}/plain/illumina_R1.fastq", r2=f"{D}/plain/illumina_R2.fastq",
    output:
        z1=f"{D}/formats/illumina_R1.fastq.zst", z2=f"{D}/formats/illumina_R2.fastq.zst",
        sp=f"{D}/formats/illumina.spring", db=f"{D}/formats/blastdb/reads.ndb",
    params: w=f"{D}/scratch"
    threads: T
    shell:
        """
        zstd -19 -T{threads} -q -f {input.r1} -o {output.z1}; zstd -19 -T{threads} -q -f {input.r2} -o {output.z2}
        mkdir -p {params.w}; spring -c -i {input.r1} {input.r2} -o {output.sp} -t {threads} -w {params.w} > /dev/null
        cat {input.r1} {input.r2} | seqkit fq2fa | makeblastdb -in - -dbtype nucl -out $(dirname {output.db})/reads \
          -title illumina_reads -max_file_sz 4GB > /dev/null
        """


rule bench_search:
    input:
        gz=[f"{D}/raw/illumina_R1.fastq.gz", f"{D}/raw/illumina_R2.fastq.gz"],
        zst=[f"{D}/formats/illumina_R1.fastq.zst", f"{D}/formats/illumina_R2.fastq.zst"],
        sp=f"{D}/formats/illumina.spring", db=f"{D}/formats/blastdb/reads.ndb",
        bam=f"{D}/aln/MG1655/illumina_raw_fullx.bam", cram=f"{D}/formats/illumina.cram",
        arch=f"{D}/formats/illumina.archive.cram", ref=f"{D}/raw/ref_MG1655.fa",
        fai=f"{D}/raw/ref_MG1655.fa.fai", stats=f"{RES}/tables/read_stats.tsv",
    output: q=f"{RES}/tables/search_latency.tsv", s=f"{RES}/tables/search_sizes.tsv"
    params:
        bases=lambda w, input: sum(int(l.split("\t")[4]) for l in open(input.stats).read().splitlines()[1:]
                                   if "raw/illumina_R" in l.split("\t")[0]),
        c=config["search"], w=f"{D}/scratch",
    threads: workflow.cores   # runs alone, so timings are not disturbed by other jobs
    benchmark: f"{RES}/benchmarks/bench_search.tsv"
    shell:
        "python workflow/scripts/bench_search.py --fastq-gz {input.gz} --fastq-zst {input.zst} "
        "--spring {input.sp} --bam {input.bam} --cram {input.cram} --cram-archive {input.arch} "
        "--blastdb $(dirname {input.db})/reads --ref {input.ref} --region {params.c[gene_region]} "
        "--bases {params.bases} --random-regions {params.c[random_regions]} --region-len {params.c[region_len]} "
        "--threads {T} --workdir {params.w} --out {output.q} --sizes {output.s}"
