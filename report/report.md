# The Storage Apocalypse

*AITU, Introduction to Bioinformatics, Project 09. Student A (Aseke): sections 1-3;
student B (KeonMei): sections 4-7. Every number comes from a file in `results/`
that `make repro` regenerates.*

## 1. Biological setting

<!-- student A -->

## 2. Data sources

<!-- student A. Part-B inputs (accessions, checksums, retrieval dates): results/tables/downloads_student_b.tsv
     and the Data table in README.md -->

## 3. Platforms and throughput

<!-- student A -->

## 4. What a sequencing run stores, and how small it can get (short reads)

**Question.** A FASTQ record has three parts: a read name, the bases, and one quality score per
base. Only the bases carry the biology. The quality scores carry the instrument's confidence,
and variant callers use them to weigh evidence. If quality scores make up most of an archive's
bytes, the cheapest saving is also the one most likely to damage later analyses. We measured
both sides of that trade on one real run.

**Data, and two problems found in the archive.** Illumina paired-end run SRR25629153
(*E. coli* K-12 MG1655, BioSample SAMN36977599): 832,053 pairs of 251 bp, 400 Mb, ~86× of the
4.64 Mb genome. The archives hold two different copies of this run.

* NCBI's *SRA Normalized* file (141 MB, MD5 from the SRA run record). We decoded it with
  `fasterq-dump`. Its qualities take only four values (Q2, Q11, Q25, Q37; 90.7% of bases
  are Q37). The instrument itself had already binned them.
* ENA's FASTQ (148 MB). Every quality in it is `?` (Q30). This is *SRA Lite*: the archive
  replaced all qualities with one constant. Nothing in the ENA record says so. We found it
  from the quality histogram, after seqkit reported exactly 100% Q30 bases.

The metadata also contradicts the data. SRA records the instrument as *Illumina HiSeq 4000*,
but the read names carry flow-cell ID `HWFV2DRXY`, a NovaSeq flow-cell pattern. 2×251 reads
with 4-level qualities are also a NovaSeq signature. We interpret the run as NovaSeq and keep
the conflict on record (Limitations). Either way, an archive user asking "what qualities does
this run have?" gets three answers: none (ENA), four levels (NCBI), and the instrument's
original values, which no longer exist anywhere.

**QC (fastp 1.3.7, thresholds in `config.yaml`).**

* Adapters were detected per pair. Only 774 reads carried adapter, as expected for 251 bp reads
  from long fragments.
* 3' ends were trimmed where a sliding window falls below Q20, the quality decay of the last
  cycles.
* A read was discarded if more than 40% of its bases are below Q15 (344 reads) or if it is
  shorter than 50 bp after trimming (18,188 reads). Shorter fragments map ambiguously in a
  genome with ~5 kb rRNA operon repeats.
* Duplicates were removed (5.1%).

Of 1,664,106 reads, 1,560,708 remained (376 Mb, 90.1% of bases ≥Q30). The raw reads were then
aligned to MG1655. 99.90% mapped and 99.69% were properly paired, with 99.999% of the genome
covered at a mean depth of 86×.

**Where the bytes are.** We compressed the three FASTQ streams separately (zstd -19). For
Illumina, read names make up 10% of the compressed bytes, bases 56% and qualities 34%. For
Nanopore, qualities make up 78% (§5). On this run qualities are "only" a third, because the
instrument already reduced them to four symbols. Where qualities still take ~40 values, as in
the Nanopore run, they dominate the archive.

**Compression benchmark (Figure B1, `compress_illumina.tsv`).** Every codec ran on the same
uncompressed FASTQ (two files, 965 MB). We timed each run with GNU `time` (wall clock, peak
memory), with nothing else running on the machine. Every output was decompressed and compared
with the input by MD5, and all codecs reproduced the input exactly. The one exception is
Spring: `fasterq-dump` repeats the read name on each record's `+` line, and Spring restores
that line as a bare `+`. Spring was therefore checked with this redundant copy ignored, so its
names, bases and qualities are identical. General-purpose codecs
compress each mate file. Spring is a FASTQ-specific compressor: it reorders reads to exploit
their overlaps and models qualities separately, and it receives both mates together.

| codec | ratio | bits/base | compress (MB/s) | decompress (MB/s) | peak memory (MB) |
|---|---|---|---|---|---|
| gzip -6 (the format archives serve) | 4.98× | 4.13 | 16 | 506 | 4 |
| pigz -9, 8 threads | 5.21× | 3.94 | 32 | 444 | 7 |
| bzip2 -9 | 6.24× | 3.29 | 31 | 62 | 8 |
| xz -6, 8 threads | 7.24× | 2.84 | 14 | 1,351 | 991 |
| zstd -3 | 4.71× | 4.36 | **2,282** | 1,194 | 116 |
| zstd -19 | 7.28× | 2.82 | 10 | **1,556** | 1,047 |
| zstd -19 --long=27 (128 MB window) | 12.29× | 1.67 | 8 | 1,316 | 1,294 |
| **Spring** (both mates) | **25.66×** | **0.80** | 99 | 104 | 1,778 |

The general-purpose codecs fall into a clear hierarchy. zstd -3 is the fastest by two orders
of magnitude but no better than gzip. A long matching window (zstd `--long`) almost halves
the size, because at 86× coverage every genome position recurs in ~86 reads that a 128 MB
window can still see. Spring, which reorders reads so that overlapping ones sit next to each
other, is 5× smaller than the gzip files archives distribute (part of the gap, a few percent, is
the redundant read-name copy on `+` lines that Spring drops), and it compresses ten times
faster than zstd -19. Its costs are memory (1.8 GB), slow decompression (104 MB/s) and no
random access (§6).

Binning on top of Spring (`binning/compress_illumina_*.tsv`) takes the archive from 0.81
bits/base (instrument 4-level) to 0.49 (2 levels) and 0.29 (no qualities). Relabelling the
four levels (8-level and 4-level maps) changes nothing, because no information is removed.
ENA's SRA Lite copy, compressed with Spring, comes to 0.24 bits/base (76.7× smaller than its
own uncompressed FASTQ, `compress_illumina_lite.tsv`): a whole bacterial run in 11.8 MB.

For aligned data (`aligned_formats_illumina.tsv`), the same reads as a sorted BAM take
105.0 MB. As CRAM 3.1 they take 36.7 MB, encoded in 1.1 s. CRAM stores only differences from
the reference, so the bases almost vanish. The CRAM "archive" profile shaves off another 5%,
at 7× the encoding time.

![Figure B1](../results/figures/fig_b1_compression.png)
*Figure B1. Archive size (bits per sequenced base, lower is better) against compression speed (MB of uncompressed FASTQ per second, log scale) for each codec; left Illumina, right Nanopore. Every point passed a lossless round-trip check.*

![Figure B2](../results/figures/fig_b2_streams_qualities.png)
*Figure B2. Left: share of the compressed FASTQ taken by read names, bases and qualities (each stream zstd -19). Right: distribution of quality values in the raw reads; the Illumina run has four values only.*

**Quality binning and what it costs (Figure B3, `variant_benchmark.tsv`).** We wrote five
versions of the QC'd reads (`bin_qualities.py`):

* the original (instrument 4-level) qualities;
* Illumina's published 8-level map;
* a 4-level map;
* a 2-level good/bad flag;
* no qualities: every base Q30, which is what SRA Lite does.

We compressed each version with Spring and called variants from it at four depths: 5×, 15×,
30× and the full ~80×. Subsamples used a fixed seed, so every scheme sees the same reads.
Reads were aligned to *E. coli* B REL606 with bwa-mem2 2.3 and variants called with bcftools
1.24 (`mpileup -q20 -Q13 | call -mv --ploidy 1`, since a bacterium is haploid).

The truth set does not come from reads at all. We aligned the finished MG1655 genome to the
finished REL606 genome (minimap2 `asm5`). paftools.js then called the differences inside
one-to-one alignment blocks of at least 10 kb: 25,952 SNPs and 354 indels over 4.02 Mb of
"confident" sequence. Reads from a K-12 strain aligned to REL606 must rediscover exactly these
differences. This gives an answer key that is large and independent of the platform. Its
known weakness: the sequenced strain is a laboratory derivative of MG1655, so its own few
private mutations count as false positives, equally in every callset.

Result:

* At 15×, 30× and full depth, no binning scheme, including discarding the qualities
  altogether, changes SNP F1 by more than 0.001. At full depth all schemes find 25,702–25,703
  of 25,952 true SNPs with 237–240 false positives (F1 0.9906).
* Only at 5× does information loss show, and it shows in score calibration rather than in the
  calls themselves. Without a QUAL filter, the no-quality calls fall from F1 0.9665 to 0.9506,
  because bcftools gives Q30 to every error and false calls stop looking weak. At the usual
  QUAL ≥20 filter, the difference disappears (0.9220 vs 0.9265).
* Indel F1 (~0.72) is limited by the aligner and caller, not by the qualities.

This is a negative result, and an important one. For haploid SNP calling at normal depth, the
quality stream is close to dead weight: 34% of this archive (78% for Nanopore) buys almost
nothing. Where it would matter (low depth, mixed or somatic samples, diploid heterozygous
calls, base-quality recalibration), this experiment cannot speak (Limitations).

![Figure B3](../results/figures/fig_b3_binning_cost.png)
*Figure B3. Left: Spring archive size per sequenced base under each quality-binning scheme. Middle and right: SNP F1 against the assembly-derived truth set (QUAL ≥20) by sequencing depth; the y axis is the same in both panels.*

## 5. Long reads and the raw signal

**Long-read data.** Nanopore MinION run SRR25637822 comes from the same BioSample as the
Illumina run: 41,926 reads, 587 Mb (~127×), read N50 15.2 kb, median read quality Q15.9
(NanoPlot). Short-read QC does not transfer. There are no adapters at fixed positions to trim,
duplicate rate is meaningless for single-molecule reads, and an Illumina-style per-base Q20
cut would discard 26% of all bases (only 74% are ≥Q20). We filtered on read length (≥1 kb) and mean quality (≥Q10). No read was
removed, because the shortest read is 5,063 bp: the submitter had already filtered before
upload, which the metadata does not say.

**Error profiles, measured rather than assumed (Figure B4, `error_profile_*.tsv`).** We
aligned the raw reads of both platforms to the correct reference (MG1655) and walked every
alignment base by base (`read_error_profile.py`; 200,000 Illumina reads and 20,000 Nanopore
reads, reservoir-sampled).

| | Illumina | Nanopore MinION |
|---|---|---|
| identity (matches / aligned columns) | 99.64% (Q24.5) | 96.97% (Q15.2) |
| median per-read identity | 100% | 97.5% |
| mismatches per 100 aligned bases | 0.354 | 1.056 |
| inserted bases per 100 | 0.002 | 0.730 |
| deleted bases per 100 | 0.001 | 1.241 |
| indels as share of errors | 0.8% | 65.1% |
| indel events in homopolymers ≥4 | 29% (of very few) | 16.5% |
| genome inside homopolymers ≥4 | 5.9% | 5.9% |

The two platforms fail differently. Illumina errors are substitutions; it practically never
gains or loses a base. Nanopore makes three times more substitutions, but two thirds of its
errors are deletions and insertions. Indel events are enriched in homopolymers (16.5% of events
vs 5.9% of the genome): the pore cannot count identical bases passing through it. This is why
indel calling from these reads stays poor (F1 0.54) while SNP calling does not (below).

![Figure B4](../results/figures/fig_b4_error_profile.png)
*Figure B4. Left: errors per 100 aligned bases by type for raw reads aligned to MG1655 (and, for the plasmid POD5 run, after hac and sup basecalling). Right: share of indel events inside homopolymers of ≥4 bases (bars) against the share of the genome in such homopolymers (black line).*

**Variants from long reads.** Nanopore reads were aligned with minimap2 2.31 (`-x map-ont`:
k=15 seeds and gap costs for an indel-dominated error model). Variants were called with
bcftools' Nanopore preset (`-X ont-sup-1.20`, which raises indel thresholds in homopolymers),
against the same assembly truth.

* SNP F1 at full depth is 0.9905, the same as Illumina's 0.9906, reached the other way round.
* Recall is higher (99.60% vs 99.04%): 15 kb reads span repeats that 251 bp reads cannot place.
* Precision is lower (98.52% vs 99.08%): residual systematic errors.
* At 15× Nanopore still gives F1 0.991, so the long-read run could be subsampled by 8× with no
  loss for this purpose.
* Removing the Nanopore qualities does not move SNP F1 (0.9903 vs 0.9905), although they are
  78% of the compressed file.

**Compression of long reads (`compress_ont.tsv`).** On the uncompressed Nanopore FASTQ (1.18 GB), every codec reproduced the input exactly:

| codec | ratio | bits/base | compress (MB/s) | decompress (MB/s) | peak memory (MB) |
|---|---|---|---|---|---|
| gzip -6 (what archives serve) | 2.06× | 7.80 | 13 | 270 | 4 |
| pigz -9, 8 threads | 2.07× | 7.74 | 44 | 265 | 8 |
| bzip2 -9 | 2.43× | 6.59 | 22 | 33 | 8 |
| xz -6, 8 threads | 2.59× | 6.19 | 10 | 591 | 1,086 |
| zstd -3 | 2.04× | 7.86 | **1,239** | 885 | 136 |
| zstd -19 | 2.54× | 6.32 | 7 | 616 | 1,130 |
| zstd -19 --long=27 | 3.05× | 5.25 | 6 | **705** | 1,612 |
| Spring (long-read mode) | **3.18×** | **5.03** | 43 | 46 | 2,722 |

Long reads compress far worse than short reads (5.03 vs ~0.8 bits/base for Illumina). The
reason is the quality stream: Nanopore qualities take ~40 distinct values that vary from base
to base and are close to incompressible, so they make up 78% of the compressed bytes. (2.8% of the bases carry Q90, far outside any
basecaller's calibrated range: a placeholder from some upstream tool that the metadata does
not explain. We left it untouched.) Binning
them is therefore where the savings are (`binning/compress_ont_*.tsv`):

| Nanopore qualities | Spring bits/base | archive size vs original | SNP F1 at full depth |
|---|---|---|---|
| original (~40 levels) | 5.03 | 1.00× | 0.9905 |
| 4 levels | 1.88 | 0.37× | 0.9906 |
| 2 levels | 1.49 | 0.30× | 0.9905 |
| none | 1.13 | 0.22× | 0.9903 |

Reducing Nanopore qualities to four levels shrinks the archive 2.7× with no measurable change
in SNP calls. This is the largest saving per unit of biological risk in the whole project.

**The raw signal: can it be deleted after basecalling? (Figure B6).** A Nanopore run's primary
data is the ionic-current trace (POD5). A neural network infers the bases from it, and that
network improves every year. We took one real POD5 file from Oxford Nanopore's open plasmid
dataset (`plasmid_2025.04`, flow cell FBC24981, R10.4.1, 5 kHz sampling, 54,946 reads,
2.06 GB) and measured the following.

* *Size.* The file holds 2.45 billion current samples. Stored as raw 16-bit integers they would
  take 4.90 GB; POD5's VBZ codec stores them at 6.74 bits per sample (2.37×). zstd -19 on top
  gains 0.16% in 84 s, so the signal is already at its practical compression limit.
* *Archive scale.* Across the whole public datasets, raw signal is 272.5 GB against 47.9 GB of
  uncompressed hac basecalls for the plasmid run (5.7×), and 273.1 GB against 39.5 GB for the
  pathogen-surveillance run (6.9×).
* *Per base.* For the 20,000 basecalled reads (65.2 Mb), the signal takes 11.5 bytes per
  base. As gzipped FASTQ the basecalls take 0.95 bytes per base (hac) or 0.77 (sup, whose
  cleaner qualities compress better). Keeping the signal therefore multiplies storage by
  12–15× over keeping the reads.
* *What the signal is worth: re-basecalling.* We basecalled the same 20,000 reads twice with
  dorado 1.4.0 on a consumer GPU (RTX 4060 Ti): once with the fast "hac" model and once with the
  slower "sup" model. Accuracy was measured against the known plasmid sequences, so no biology
  enters the number.

| | hac | sup |
|---|---|---|
| median per-read identity | 98.73% | **99.47%** |
| mismatches per 100 aligned bases | 0.82 | 0.45 |
| deleted bases per 100 | 1.92 | 1.67 |
| GPU time for 20,000 reads | 57 s | 308 s (5.4×) |

  Re-reading the *same* signal with a better model removes 58% of the errors of a typical
  read (median error 1.27% → 0.53%) and halves substitutions. Only a lab that kept the POD5
  files can do this when a better model ships. FASTQ is a frozen snapshot of what the
  basecaller knew on the day.

**Recommendation.** Signal costs 12–15× the reads and cannot be compressed further, yet it
can still yield a large accuracy gain. Neither "always delete" nor "always keep" is right.
Keep signal for reference samples, for clinical or legal cases and for novel organisms. For
routine resequencing, keep it until the next major basecaller release has been applied
(typically about a year), in deep archive, where at $0.00099/GB-month a year of a 2 GB file
costs $0.02. Then delete it.

![Figure B6](../results/figures/fig_b6_signal.png)
*Figure B6. Left: total size of raw signal (POD5) and basecalls (FASTQ) in two complete public ONT datasets. Right: per-read accuracy (Phred-scaled identity to the known plasmid sequence) of the same 20,000 signal traces basecalled with the hac and sup models.*

## 6. Alignment as an index: is a compressed archive still searchable?

**Question.** Compression only helps if people can still find things in the archive. We stored
the same Illumina run in six forms and asked each one the same questions
(`bench_search.py`, three repeats, median wall time, nothing else running):

* find every read that comes from one gene (*lacZ*, NC_000913.3:363,231–366,305);
* fetch reads from 200 random 1 kb windows;
* read every record.

For the FASTQ forms there are no coordinates, so we scanned for the gene's 31-mers on both
strands (154 k-mers × 2, GNU `grep -F`, one pass). BLAST searched a database built from the
reads (`blastn`, E < 1e-20). Recall is measured against the reads that the aligned BAM places
on the gene.

| storage form | bits/base | gene query | recall | random 1 kb region | full scan |
|---|---|---|---|---|---|
| FASTQ.gz (fasterq-dump, pigz) | 4.13 | 5.27 s | 0.982 | not possible | 2.50 s |
| FASTQ.zst (-19) | 2.82 | 5.20 s | 0.982 | not possible | 0.83 s |
| Spring | 0.80 | 15.4 s | 0.982 | not possible | 10.2 s |
| BAM (sorted, indexed) | 2.10 | **0.006 s** | 1 (reference) | 4.8 ms | 0.26 s |
| CRAM 3.1 | **0.73** | 0.015 s | 1.000 | 6.3 ms | 0.10 s |
| CRAM 3.1 archive profile | 0.70 | 0.34 s | 1.000 | 32 ms | 0.15 s |
| BLAST database of reads | 7.28 | 0.15 s | 0.972 | not possible | – |

![Figure B5](../results/figures/fig_b5_search_latency.png)
*Figure B5. Left: storage per sequenced base against the time to retrieve the reads of one gene (log scale), with recall against the BAM answer. Right: time to read every record in each form.*

**Interpretation.**

* *An alignment is the index.* CRAM is the smallest form we measured, smaller even than Spring
  (0.73 vs 0.80 bits/base), because it stores only differences from the reference. It also
  answers a gene query in 15 ms, 360× faster than scanning any compressed FASTQ.
* *Index lookup vs scan.* CRAM, like BAM, keeps an index of genome coordinates (`.crai`) and
  decodes only the blocks it needs. Reading a random 1 kb window costs 6 ms in CRAM against
  5 ms in BAM: random access survives compression.
* *Where the costs hide.* The "archive" CRAM profile buys 5% fewer bytes for 5–20× slower
  queries. Spring is excellent for cold storage but has no random access: every question costs
  a full 10 s decompression.
* *BLAST.* The heuristic seed-and-extend search over an indexed word list is 35× faster than
  scanning FASTQ, but the database is the largest form here: it keeps sequence headers and
  lookup tables, and no qualities. It is a search accelerator, not an archive.
* *Recall.* The k-mer scan misses 1.8% of the gene's reads: those overlapping its ends by fewer
  than 31 bases, which an exact k-mer cannot see. BLAST misses 2.8%: short overlaps fail the
  E-value cutoff. Both algorithms trade sensitivity at the edges for speed, and the alignment-
  based answer is the only exact one.

The trade-off for an archive follows directly. Data that will be queried should be kept as
CRAM against a stable reference. FASTQ in a FASTQ-specific compressor suits data that will be
restored whole if at all. The catch with CRAM is that the reference becomes part of the
archive, and losing it makes the data unreadable.

## 7. Tiering: which data can go cold?

**Question.** Most archived sequencing data is never read again, but nobody knows *which* data
in advance. A tiering policy has to decide at submission time, from metadata alone, which
studies can go to deep storage (cheap, but hours to restore) and which must stay
instantly readable.

**Reuse signal.** Download counts are not published by SRA or ENA, so we used the next best
public signal: papers that mention the study. For every *E. coli* study we queried Europe PMC
with both the text-mined accession field and the literal accession string
(`(ACCESSION_ID:"PRJNA…" OR "PRJNA…")`), for the BioProject and its SRP/ERP alias. One paper
is usually the submitters' own description of the data; a second paper means somebody came
back. A study counts as **reused** if at least two papers mention it between the year before
and four years after its release. The fixed window matters: counting all papers to date would
teach the model that "old = reused" simply because old studies had more years to be cited.

**Data at scale.** We took the complete ENA run table for *E. coli* (`tax_tree(562)`).

* *Download.* A single request was silently cut off by the server at 383,280 of ~608,000 rows,
  with no error raised. The fetcher therefore downloads one release year at a time and checks
  each year's row count against ENA's own `/count` endpoint. The result is 607,211 runs, 0.27 PB
  of FASTQ, 9,875 studies, held as Parquet (101 MB TSV → 17 MB).
* *Metadata problems* (`tiering_metadata_issues.tsv`). 261 runs have no study, 2,638 have no
  file size, and 3,751 have zero bases because they were uploaded as native ONT/PacBio files
  that ENA never converted.
* *Labels.* The 5,748 studies released 2010–2021 (155 TB) needed 11,497 Europe PMC queries
  (27.8 min at 8.4 queries/s, resumable cache, 0.6 GB peak memory). 30% of the studies are
  mentioned in at least one paper and 8.0% are reused.

**Who gets reused (`tiering_reuse_by_group.tsv`, Figure B7).**

* Reuse rises steeply with release year: 0.7% for 2012 studies and 9–12% for 2019–2021 studies.
  Part of this is real, as surveillance-era data is analysed more. Part is a measurement
  artefact: citing BioProject accessions in papers became common only after ~2014.
* Long-read studies are cited more often than Illumina ones (45% and 39% vs 28% cited at
  least once) but are 7% of the bytes.
* Whole-genome shotgun data is 85% of the bytes.

**Model.** The model sees only features available at submission:

* data volume, number of runs and samples, bytes per run;
* platform shares and paired-end share;
* library strategy and source;
* instrument model and submitting centre (top categories, learned on training years only);
* release year.

Validation is temporal, as in real use: we fitted on studies released ≤2017 (3,021 studies,
3.4% reused) and tested on 2018–2021 (2,727 studies, 9.0% reused).

| model | test ROC-AUC | test PR-AUC |
|---|---|---|
| logistic regression (class-balanced) | **0.752** | 0.224 |
| gradient boosting | 0.709 | 0.204 |
| study size only | 0.711 | 0.226 |
| random | 0.493 | 0.093 |

The logistic model beats chance clearly, but much of its signal is size: big studies (large
consortia, surveillance programmes) get reused. Gradient boosting overfits the training years,
whose reuse rate is a third of the test years'. The strongest coefficients are submitting
centres, for example NISC (odds ratio 7.9) and Sanger (4.0), positive, against GEO-routed
submissions (0.14) and Illumina GA II era data (0.12). These are proxies for "this is part of
a public-health programme" rather than causes.

**Cost (`tiering_costs.tsv`).** We priced each policy on the 75 TB of test-year studies over 5
years. AWS us-east-1 list prices were read from the AWS Price List API, and the script stops
if `config.yaml` disagrees with the downloaded price list. Each reuse event after the first
paper is charged as one full read of the study.

| policy | 5-year cost (USD) | saving vs S3 Standard | reuse requests delayed |
|---|---|---|---|
| all S3 Standard | 103,395 | 0% | 0% |
| all Glacier Instant Retrieval | 19,131 | **81.5%** | 0% |
| all Deep Archive | 4,546 | 95.6% | 100% |
| age rule: Standard 24 months, then Deep Archive | 44,106 | 57.3% | 73% |
| model: Deep Archive if p < 0.61, else Instant Retrieval | 17,094 | 83.5% | 17% |
| oracle (knows future reuse) | 8,446 | 91.8% | 0% |

Three conclusions follow, and the first one was not what we expected.

1. *The first decision is not hot vs cold.* Moving everything from S3 Standard to an instant-
   access archive tier saves 81.5% and delays nobody. Retrieval fees are negligible at the
   observed reuse rate ($1,149 over 5 years).
2. *The common age rule is the worst option.* Keeping data hot for two years costs 2.3× more
   than the instant archive and still delays 73% of reuse, because reuse arrives late
   (papers come 1–4 years after release).
3. *Prediction buys real but limited extra savings.* The model's threshold was chosen on the
   training years as the cheapest one that delays ≤10% of reuse. On the test years it saves a
   further 10.7% over all-Instant-Retrieval, but delays 17% of reuse requests: the training
   years reused data three times less, so the threshold is miscalibrated for the future. The
   oracle shows the ceiling: 79% of the test bytes are never reused, and knowing which would
   save 56% more than all-Instant-Retrieval. Better reuse signals (real download logs) are
   worth more than a better model.

![Figure B7](../results/figures/fig_b7_tiering.png)
*Figure B7. Left: share of E. coli studies mentioned in ≥1 and ≥2 papers within the fixed window, by release year. Middle: ROC curve of the logistic reuse model on the held-out years 2018–2021. Right: total 5-year cost saving against the share of reuse requests that must wait for a Deep Archive restore; the line sweeps the model's threshold, points are the fixed policies.*

`tiering_savings_for_forecast.tsv` passes the model policy's saving (83.5% of archive storage
cost) to the crossover forecast in section 3. The crossover date moves later by however many
years of storage-price decline a factor of 1/(1−0.835) ≈ 6× represents.


## 8. Engineering and reproducibility (part B)

* **One workflow.** Snakemake 9.27 runs about 180 jobs from raw downloads to figures
  (`make repro`). The same rules run on ~14 MB of real data slices in `test_data/`
  (`make test`, no network, no GPU).
* **Pinned environment.** Versions are pinned in `environment.yml`, with the full transitive
  lock in `environment.lock.yml`. Dorado, which is not on conda, is installed by a versioned
  script.
* **Provenance.** Every input is streamed to disk and checked against the archive's MD5 where
  one exists (ENA FASTQ, NCBI SRA file). URL, size, MD5, SHA-256 and UTC retrieval time are
  recorded (`downloads_student_b.tsv`). Prices are re-checked against the downloaded AWS price
  list on every run.
* **Scale handling.**
  * Large data never enters git and lives on the Linux filesystem.
  * The 607k-row metadata table is downloaded in verified yearly chunks and stored as Parquet.
  * Europe PMC queries run concurrently (10 threads) through a resumable cache.
  * Read profiling uses one-pass reservoir sampling, so memory stays bounded.
  * Indexes are used instead of scans wherever the format allows (§6).
* **Measured cost.** Every rule has a Snakemake `benchmark:` (wall time, CPU, memory) in
  `results/benchmarks/`. Compression and search benchmarks reserve the whole machine, so
  their timings are not disturbed by parallel jobs.
* **Failure handling.** A truncated download is never accepted: `.part` files, MD5 checks and
  row-count checks against ENA's `/count`. Codec results are verified by a lossless round
  trip.


## Limitations

<!-- student A: growth model, instrument and price limitations -->

### Part B

* **One organism, one run per platform.** All compression, binning and error-profile numbers
  come from one *E. coli* Illumina run and one MinION run.
  * Compression ratios depend on genome repetitiveness and coverage. A 3 Gb human genome at
  30× compresses differently (Spring gains more from read overlap at high coverage).
  * The binning result covers haploid SNP calling only. Diploid heterozygous calls, low-
  frequency or somatic variants, and base-quality recalibration may depend on qualities in
  ways this design cannot detect.
* **The "original" qualities were already binned.** The best copy in any archive is the
  instrument's 4-level output, and the instrument's raw quality values no longer exist. We
  could therefore test binning only from 4 levels downwards.
* **Metadata conflict.** SRA calls the Illumina run HiSeq 4000. Read names and the quality
  structure indicate NovaSeq. We could not resolve this from public records.
* **Truth set.** It is made from two reference assemblies. The sequenced K-12 strain is a
  laboratory derivative of MG1655, so its private mutations appear as false positives in
  every callset (equally). Only 4.02 of 4.70 Mb lie in confident regions, so repeats are
  excluded from the evaluation.
* **Long-read chemistry unknown.** The MinION run's chemistry and basecaller are not stated in
  its metadata. The submitter's pre-filtering (no read <5 kb) is undocumented and biases read
  length upwards.
* **Raw signal comes from a different organism.** No public run gave us POD5 for the same
  *E. coli* sample. The POD5 experiment uses ONT's plasmid dataset, chosen because its true
  sequences are known exactly.
* **Reuse is proxied by citations.** SRA and ENA publish no download counts. Paper mentions
  miss reuse that is never published and depend on text-mining coverage, which is incomplete
  before ~2014 and limited to open-access full text. The train/test reuse rates differ 3×,
  which miscalibrated the model's threshold. Retrieval costs assume each reuse reads the
  whole study once.
* **Prices.** These are list prices for one provider and region (AWS us-east-1, Sept 2026).
  Archives such as SRA negotiate their own prices, and egress ($0.09/GB) is excluded because
  it is the same in every tier.
* **Timing.** All timings come from one workstation (16 threads, WSL2, consumer GPU). Absolute
  speeds will differ elsewhere. The ranking of codecs is the transferable part.


## Contribution statement

See `CONTRIBUTIONS.md`.

## Appendix: AI assistance

See `report/ai_disclosure.md`.
