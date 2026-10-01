"""Compare variant call sets against a truth set inside confident regions (tasks 4, 5, 6).

Truth = differences between two finished genome assemblies (K-12 MG1655 vs B REL606), called by
paftools.js from a whole-genome alignment. That truth does not depend on reads, base qualities
or platform, so it can judge Illumina, binned-quality Illumina and Nanopore calls equally.

All VCFs must already be normalised (bcftools norm -f REF -m -any -a): left-aligned, split,
atomised. A variant matches when CHROM, POS, REF and ALT are identical.

Output TSV (one row per callset x type x QUAL cutoff):
  callset type min_qual tp fp fn precision recall f1
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import sys
from pathlib import Path

import pandas as pd


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    regs: dict[str, list[tuple[int, int]]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            c, s, e = line.split("\t")[:3]
            regs.setdefault(c, []).append((int(s), int(e)))
    for c in regs:
        regs[c].sort()
    return regs


def inside(regs, chrom: str, pos1: int) -> bool:
    iv = regs.get(chrom)
    if not iv:
        return False
    k = bisect.bisect_right(iv, (pos1 - 1, float("inf"))) - 1
    return k >= 0 and iv[k][0] <= pos1 - 1 < iv[k][1]


def read_vcf(path: Path, regs, pass_only: bool = True) -> dict[tuple, float]:
    """{(chrom,pos,ref,alt): qual} for ALT calls inside confident regions."""
    out = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual, filt = f[0], int(f[1]), f[3], f[4], f[5], f[6]
            if alt in (".", "*") or (pass_only and filt not in ("PASS", ".")):
                continue
            if len(f) > 9:                     # genotype present: skip reference / no-call genotypes
                gt = f[9].split(":")[0].replace("|", "/")
                if all(x in ("0", ".") for x in gt.split("/")):
                    continue
            if not inside(regs, chrom, pos):
                continue
            out[(chrom, pos, ref, alt)] = float(qual) if qual != "." else float("inf")
    return out


def vtype(key) -> str:
    return "SNP" if len(key[2]) == 1 and len(key[3]) == 1 else "INDEL"


def score(truth: dict, calls: dict, min_qual: float) -> list[dict]:
    rows = []
    q = {k for k, v in calls.items() if v >= min_qual}
    for t in ("SNP", "INDEL", "ALL"):
        T = {k for k in truth if t == "ALL" or vtype(k) == t}
        C = {k for k in q if t == "ALL" or vtype(k) == t}
        tp, fp, fn = len(T & C), len(C - T), len(T - C)
        p = tp / (tp + fp) if tp + fp else float("nan")
        r = tp / (tp + fn) if tp + fn else float("nan")
        f1 = 2 * p * r / (p + r) if p + r else float("nan")
        rows.append({"type": t, "min_qual": min_qual, "tp": tp, "fp": fp, "fn": fn,
                     "precision": p, "recall": r, "f1": f1})
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", type=Path, required=True)
    ap.add_argument("--bed", type=Path, required=True, help="confident regions of the truth set")
    ap.add_argument("--calls", type=Path, nargs="+", required=True)
    ap.add_argument("--names", nargs="+", help="callset names (default: file stem)")
    ap.add_argument("--qual-steps", default="0,10,20,30,50,100")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    regs = read_bed(a.bed)
    truth = read_vcf(a.truth, regs, pass_only=False)
    names = a.names or [p.name.split(".")[0] for p in a.calls]
    rows = []
    for name, path in zip(names, a.calls):
        calls = read_vcf(path, regs)
        for mq in (float(x) for x in a.qual_steps.split(",")):
            for r in score(truth, calls, mq):
                rows.append({"callset": name, **r})
    df = pd.DataFrame(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.round(5).to_csv(a.out, sep="\t", index=False)
    print(df[(df.min_qual == 0) & (df.type != "ALL")].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
