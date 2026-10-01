"""Per-read error profile from alignments (tasks 4, 5): what kind of errors each platform makes.

Reads `samtools view` text on stdin (primary alignments only), walks each CIGAR against the
reference and counts mismatches, inserted and deleted bases. An indel is called homopolymer-
associated when it sits in (or extends) a run of >=4 identical reference bases; Nanopore's signal
cannot tell how many identical bases passed through the pore, so it should be enriched there.

  samtools view -F 0x904 aln.bam | read_error_profile.py --ref ref.fa --platform ONT \
      --per-read reads.tsv --summary summary.tsv

identity = matches / (matches + mismatches + inserted + deleted bases)  (BLAST-style identity)
"""
from __future__ import annotations

import argparse
import math
import random
import re
import sys
from pathlib import Path

import pandas as pd

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
HP_MIN = 4


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def run_length(ref: str, pos: int) -> int:
    """Length of the homopolymer run that contains ref[pos]."""
    if pos < 0 or pos >= len(ref):
        return 0
    b, l, r = ref[pos], pos, pos
    while l > 0 and ref[l - 1] == b:
        l -= 1
    while r + 1 < len(ref) and ref[r + 1] == b:
        r += 1
    return r - l + 1


def profile_read(seq: str, ref: str, pos0: int, cigar: str) -> dict:
    mm = ins = dele = match = 0
    ins_ev = del_ev = hp_ins = hp_del = 0
    q, r = 0, pos0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                if seq[q + k] == ref[r + k]:
                    match += 1
                else:
                    mm += 1
            q += n
            r += n
        elif op == "I":
            ins += n
            ins_ev += 1
            base = seq[q]
            if (r < len(ref) and ref[r] == base and run_length(ref, r) >= HP_MIN - 1) or \
               (r > 0 and ref[r - 1] == base and run_length(ref, r - 1) >= HP_MIN - 1):
                hp_ins += 1
            q += n
        elif op in "DN":
            dele += n
            del_ev += 1
            if run_length(ref, r) >= HP_MIN:
                hp_del += 1
            r += n
        elif op == "S":
            q += n
    aligned = match + mm + ins + dele
    return {"aligned_cols": aligned, "matches": match, "mismatches": mm, "ins_bases": ins, "del_bases": dele,
            "ins_events": ins_ev, "del_events": del_ev, "hp_ins": hp_ins, "hp_del": hp_del,
            "identity": match / aligned if aligned else float("nan")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", type=Path, required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--max-reads", type=int, default=50000, help="reservoir sample size")
    ap.add_argument("--per-read", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args(argv)
    refs = read_fasta(a.ref)
    rng = random.Random(a.seed)
    sample, seen = [], 0
    for line in sys.stdin:                       # reservoir sampling: one pass, bounded memory
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        if f[5] == "*" or f[9] == "*":
            continue
        seen += 1
        rec = (f[0], f[2], int(f[3]) - 1, f[5], f[9], len(f[9]))
        if len(sample) < a.max_reads:
            sample.append(rec)
        else:
            j = rng.randrange(seen)
            if j < a.max_reads:
                sample[j] = rec
    rows = []
    for name, chrom, pos0, cigar, seq, rlen in sample:
        d = profile_read(seq, refs[chrom], pos0, cigar)
        d.update(read=name, read_len=rlen)
        rows.append(d)
    df = pd.DataFrame(rows)
    a.per_read.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.per_read, sep="\t", index=False)
    tot = df[["aligned_cols", "matches", "mismatches", "ins_bases", "del_bases", "ins_events",
              "del_events", "hp_ins", "hp_del"]].sum()
    err = 1 - tot["matches"] / tot["aligned_cols"]
    ref_all = "".join(refs.values())
    hp_share_genome = sum(len(m.group(0)) for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", ref_all)) / len(ref_all)
    summ = {
        "platform": a.platform, "reads_mapped_total": seen, "reads_profiled": len(df),
        "median_read_len": float(df["read_len"].median()),
        "median_identity": float(df["identity"].median()),
        "mean_identity": float(tot["matches"] / tot["aligned_cols"]),
        "phred_of_error": -10 * math.log10(err) if err > 0 else float("inf"),
        "mismatch_per_100": 100 * tot["mismatches"] / tot["aligned_cols"],
        "ins_per_100": 100 * tot["ins_bases"] / tot["aligned_cols"],
        "del_per_100": 100 * tot["del_bases"] / tot["aligned_cols"],
        "indel_share_of_errors": (tot["ins_bases"] + tot["del_bases"]) /
                                 max(1, tot["mismatches"] + tot["ins_bases"] + tot["del_bases"]),
        "hp_share_of_indel_events": (tot["hp_ins"] + tot["hp_del"]) / max(1, tot["ins_events"] + tot["del_events"]),
        "hp_share_of_genome": hp_share_genome,
    }
    pd.DataFrame([summ]).to_csv(a.summary, sep="\t", index=False)
    print(pd.Series(summ).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
