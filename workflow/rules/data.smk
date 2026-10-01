# Student B: inputs, references, truth set.
# Every download goes through fetch_data.py (streamed, MD5-checked, provenance sidecar).
import re

D = config["data_dir"]
RES = config.get("results_dir", "results")
T = config["threads"]
INPUT_NAMES = {
    "illumina_sra": "illumina.sra",                 # main run: NCBI SRA Normalized (with qualities)
    "illumina_R1": "illumina_R1.fastq.gz",          # test run: FASTQ given directly
    "illumina_R2": "illumina_R2.fastq.gz",
    "illumina_lite_R1": "illumina_lite_R1.fastq.gz",
    "illumina_lite_R2": "illumina_lite_R2.fastq.gz",
    "ont": "ont.fastq.gz",
    "ref_MG1655": "ref_MG1655.fa",
    "ref_REL606": "ref_REL606.fa",
    "pod5": "signal.pod5",
    "aws_s3_prices": "aws_s3_prices.csv",
}
INPUT_FILES = {k: f"{D}/raw/{v}" for k, v in INPUT_NAMES.items() if config["inputs"].get(k)}


FILE2KEY = {v.split("/")[-1]: k for k, v in INPUT_FILES.items()}


rule fetch_input:
    output: f"{D}/raw/{{name}}"
    wildcard_constraints: name="|".join(re.escape(n) for n in FILE2KEY)
    params:
        key=lambda w: FILE2KEY[w.name],
        url=lambda w: config["inputs"][FILE2KEY[w.name]]["url"],
        md5=lambda w: config["inputs"][FILE2KEY[w.name]].get("md5", ""),
        acc=lambda w: config["inputs"][FILE2KEY[w.name]]["accession"],
        note=lambda w: config["inputs"][FILE2KEY[w.name]].get("note", ""),
    shell:
        "python workflow/scripts/fetch_data.py url --url {params.url:q} --md5={params.md5:q} "
        "--out {output:q} --key {params.key:q} --accession {params.acc:q} --note {params.note:q}"


if "illumina_sra" in INPUT_FILES:
    rule sra_to_fastq:
        """fasterq-dump keeps the original base qualities of the SRA Normalized file."""
        input: INPUT_FILES["illumina_sra"]
        output: r1=f"{D}/raw/illumina_R1.fastq.gz", r2=f"{D}/raw/illumina_R2.fastq.gz"
        threads: 8
        benchmark: f"{RES}/benchmarks/sra_to_fastq.tsv"
        shell:
            """
            d=$(dirname {output.r1})/sra_dump; mkdir -p $d
            fasterq-dump --split-files --threads {threads} --outdir $d --temp $d {input} > /dev/null
            pigz -p {threads} -c $d/*_1.fastq > {output.r1}; pigz -p {threads} -c $d/*_2.fastq > {output.r2}; rm -rf $d
            """


rule fetch_plasmid_refs:
    output: f"{D}/raw/plasmid_refs.fa"
    params: b=config["plasmid_refs_s3"]["bucket"], p=config["plasmid_refs_s3"]["prefix"],
            local=config["plasmid_refs_s3"].get("local", "")
    shell:
        "if [ -n '{params.local}' ]; then python workflow/scripts/fetch_data.py url --url local:{params.local} "
        "--out {output} --key plasmid_refs --accession 'test subset' --note 'local test copy'; else "
        "python workflow/scripts/fetch_data.py s3cat --bucket {params.b} --prefix {params.p} "
        "--out {output} --key plasmid_refs --accession 'ont-open-data plasmid_2025.04 references' "
        "--note 'one full plasmid sequence per barcode'; fi"


rule provenance:
    input: [v for v in INPUT_FILES.values()] + [f"{D}/raw/plasmid_refs.fa"]
    output: f"{RES}/tables/downloads_student_b.tsv"
    params: prov=lambda w, input: " ".join(f"'{p}.prov.tsv'" for p in input)
    shell: "python workflow/scripts/fetch_data.py merge --prov {params.prov} --out {output}"


rule index_ref:
    input: f"{D}/raw/ref_{{ref}}.fa"
    output:
        fai=f"{D}/raw/ref_{{ref}}.fa.fai",
        bwa=f"{D}/raw/ref_{{ref}}.fa.bwt.2bit.64",
    benchmark: f"{RES}/benchmarks/index_ref_{{ref}}.tsv"
    shell: "samtools faidx {input} && bwa-mem2 index {input} > /dev/null 2>&1"


rule truth_set:
    """Differences between the two finished assemblies = truth for every read-based callset.
    -l/-L 10000: only 1-to-1 alignment blocks >=10 kb count as 'confident' (repeats excluded)."""
    input:
        rel=f"{D}/raw/ref_REL606.fa", mg=f"{D}/raw/ref_MG1655.fa", fai=f"{D}/raw/ref_REL606.fa.fai",
    output:
        vcf=f"{D}/truth/truth.vcf.gz", bed=f"{RES}/tables/truth_confident.bed",
        stats=f"{RES}/tables/truth_summary.tsv",
    threads: T
    shell:
        """
        d=$(dirname {output.vcf})
        minimap2 -cx asm5 --cs -t {threads} {input.rel} {input.mg} 2>/dev/null | sort -k6,6 -k8,8n > $d/asm.paf
        paftools.js call -l 10000 -L 10000 $d/asm.paf 2>/dev/null | awk '$1=="R"' | cut -f2-4 > {output.bed}
        paftools.js call -l 10000 -L 10000 -f {input.rel} -s MG1655 $d/asm.paf 2>/dev/null \
          | bcftools norm -f {input.rel} -m -any -a 2>/dev/null | bcftools sort -Oz -o {output.vcf} 2>/dev/null
        bcftools index -f {output.vcf}
        printf "confident_bp\\tsnps\\tindels\\n" > {output.stats}
        paste <(awk '{{s+=$3-$2}} END {{print s}}' {output.bed}) \
              <(bcftools view -H {output.vcf} | awk 'length($4)==1&&length($5)==1' | wc -l) \
              <(bcftools view -H {output.vcf} | awk 'length($4)!=length($5)' | wc -l) >> {output.stats}
        """
