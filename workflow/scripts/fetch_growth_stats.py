"""Download GenBank release notes (gbrel.txt), cut out the growth tables, write one tidy CSV.

Steps
  1. download gbrel.txt (streamed) -> data/raw/gbrel_<release>.txt   (release read from file header)
  2. record URL, release, bytes, SHA-256, UTC date in config/accessions.tsv
  3. locate every growth table, tag it with a component (GenBank / WGS / TSA / TLS), parse rows
  4. write results/tables/genbank_growth_gbrel_<release>.csv
     columns: component, release, date, bases, entries, note, date_file, date_precision
       date            canonical date: the date the GenBank (traditional) table gives for that release number
       date_file       the date exactly as printed in this row's own table (differs from `date` only where
                       NCBI's tables disagree; such rows also say so in `note`)
       date_precision  'month' or 'day' (gbrel.txt gives month+year only, so its dates are the 1st of the month)

How tables are located (no hard-coded line numbers, the file changes every release):
  * a table HEADER is a line containing 'release', 'date', 'bases'/'base pairs' and one of
    'entries'/'sequences'/'records'/'loci'
  * data ROWS follow: <release> <date> <bases> <entries> [note]
  * the table ENDS at the first line that is not a row (after at least one row)
  * the COMPONENT is taken from the nearest preceding title lines (see COMPONENT_PATTERNS)

A release that was not delivered (e.g. Feb 2024) is kept as a row with empty release/bases/entries
and the file's explanation in `note`, so gaps stay visible in the time series.

Use --input FILE --inspect to print what the locator sees without writing anything.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import re
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import requests

URL = "https://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt"
EXPECTED_COMPONENTS = ("GenBank", "WGS", "TSA", "TLS")
OUT_COLUMNS = ["component", "release", "date", "bases", "entries", "note", "date_file", "date_precision"]
ACC_COLUMNS = ["accession_or_url", "source", "release", "bytes", "sha256",
               "download_date", "verified_by", "note"]

# checked nearest-title-line first; first match wins
COMPONENT_PATTERNS = [
    ("WGS", re.compile(r"\bWGS\b|whole\s+genome\s+shotgun", re.I)),
    ("TSA", re.compile(r"\bTSA\b|transcriptome\s+shotgun", re.I)),
    ("TLS", re.compile(r"\bTLS\b|targeted\s+locus", re.I)),
    ("GenBank", re.compile(r"genbank|traditional", re.I)),
]

ROW_RE = re.compile(
    r"^\s*(?P<rel>\d+(?:\.\d+)?|n/a)\s+"
    r"(?P<date>[A-Za-z]{3,9}\.?\s+\d{1,2},?\s+\d{4}|[A-Za-z]{3,9}\.?\s+\d{4}"
    r"|\d{1,2}/\d{1,2}/\d{2,4}|\d{1,2}/\d{4}|\d{4}-\d{2}(?:-\d{2})?)\s+"
    r"(?P<bases>\d[\d,]*|n/a)\s+(?P<entries>\d[\d,]*|n/a)"
    r"(?:\s+(?P<note>\S.*?))?\s*$"
)
RELEASE_RE = re.compile(r"Release\s+(\d+\.\d+)")
SEPARATOR_RE = re.compile(r"^\s*[-=_*~ ]{3,}\s*$")
log = logging.getLogger("fetch_growth_stats")


@dataclass
class Table:
    header_line: int                      # 1-based
    component: str
    title: str
    rows: list[dict] = field(default_factory=list)


# ---------------------------------------------------------------- download
def download(url: str, raw_dir: Path, retries: int = 3) -> tuple[Path, str, int, str]:
    """Stream to disk, hash on the fly, name the file by release. Returns (path, release, bytes, sha256)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    tmp = raw_dir / "gbrel.txt.part"
    for attempt in range(1, retries + 1):
        try:
            sha, n = hashlib.sha256(), 0
            with requests.get(url, stream=True, timeout=60,
                              headers={"User-Agent": "aitu-bioinfo-course-project/1.0"}) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        fh.write(chunk)
                        sha.update(chunk)
                        n += len(chunk)
            break
        except requests.RequestException as exc:
            log.warning("download attempt %d/%d failed: %s", attempt, retries, exc)
            if attempt == retries:
                raise
            time.sleep(2 ** attempt)
    release = detect_release(tmp.read_text(errors="replace").splitlines())
    final = raw_dir / f"gbrel_{release}.txt"
    tmp.replace(final)
    log.info("downloaded %s (%d bytes, release %s, sha256 %s)", final, n, release, sha.hexdigest())
    return final, release, n, sha.hexdigest()


def detect_release(lines: list[str]) -> str:
    for line in lines[:80]:
        m = RELEASE_RE.search(line)
        if m:
            return m.group(1)
    raise ValueError("no 'Release NNN.N' found in the first 80 lines; the file layout has changed")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- parsing
def is_header(line: str) -> bool:
    t = re.sub(r"\s+", " ", line.lower())
    return ("release" in t and "date" in t and ("bases" in t or "base pairs" in t)
            and any(w in t for w in ("entries", "sequences", "records", "loci")))


def parse_date(raw: str) -> tuple[str, str, str]:
    """Return (ISO date, precision, note). Month+year input gives the 1st of the month with
    precision 'month'. note is only set for unparseable dates."""
    s = raw.replace(".", "").replace(",", "").strip()
    for fmt, precision in (("%b %d %Y", "day"), ("%B %d %Y", "day"), ("%b %Y", "month"), ("%B %Y", "month"),
                           ("%m/%d/%Y", "day"), ("%m/%d/%y", "day"), ("%m/%Y", "month"),
                           ("%Y-%m-%d", "day"), ("%Y-%m", "month")):
        try:
            return datetime.strptime(s, fmt).date().isoformat(), precision, ""
        except ValueError:
            continue
    return "", "", f"unparsed date '{raw}'"


def classify(lines: list[str], header_idx: int, prev_end: int) -> tuple[str, str]:
    """Component = the keyword that appears FIRST in the intro text between the previous
    table (or 40 lines back) and this header. Nearest-line matching fails on the real file,
    where the lines right above a header are footnotes that mention other components."""
    start = max(prev_end, header_idx - 15)        # fallback
    for k in range(header_idx - 1, max(prev_end, header_idx - 80) - 1, -1):
        if "following table" in lines[k].lower():
            start = k
            break
    text = "\n".join(lines[start:header_idx])
    best = None                                   # (position, name, matched text)
    for name, pat in COMPONENT_PATTERNS:
        m = pat.search(text)
        if m and (best is None or m.start() < best[0]):
            best = (m.start(), name, m.group(0))
    if best is None:
        return "UNCLASSIFIED", ""
    ctx = text[max(0, best[0] - 30): best[0] + 30].replace("\n", " ").strip()
    return best[1], ctx


def find_tables(lines: list[str]) -> list[Table]:
    tables, i, prev_end = [], 0, 0
    while i < len(lines):
        if not is_header(lines[i]):
            i += 1
            continue
        comp, title = classify(lines, i, prev_end)
        tbl = Table(header_line=i + 1, component=comp, title=title)
        j, skipped = i + 1, 0
        while j < len(lines):
            m = ROW_RE.match(lines[j])
            if m:
                iso, precision, dnote = parse_date(m["date"])
                notes = "; ".join(x for x in (dnote, (m["note"] or "").strip()) if x)
                tbl.rows.append({
                    "component": comp, "release": "" if m["rel"] == "n/a" else m["rel"], "date": iso,
                    "date_file": iso, "date_precision": precision,
                    "bases": None if m["bases"] == "n/a" else int(m["bases"].replace(",", "")),
                    "entries": None if m["entries"] == "n/a" else int(m["entries"].replace(",", "")),
                    "note": notes,
                })
            elif not tbl.rows and skipped < 3 and (not lines[j].strip() or SEPARATOR_RE.match(lines[j])):
                skipped += 1                    # blank / rule lines between header and first row
            else:
                break
            j += 1
        if tbl.rows:
            tables.append(tbl)
            prev_end = j
        i = max(j, i + 1)
    return tables


def reconcile_dates(tables: list[Table], reference: str = "GenBank") -> list[str]:
    """Give every release ONE date. NCBI's tables sometimes disagree (release 272 is Jun 2026 in the
    GenBank and WGS tables but Apr 2026 in TSA and TLS). The GenBank table is the reference: releases
    are bimonthly and that table is the one the release notes are about. The original value stays in
    `date_file` and the correction is written to `note`. Returns a list of human-readable corrections."""
    ref = {r["release"]: r["date"] for t in tables if t.component == reference
           for r in t.rows if r["release"]}
    fixes = []
    for t in tables:
        if t.component == reference:
            continue
        for r in t.rows:
            if r["release"] in ref and r["date"] != ref[r["release"]]:
                msg = (f"date corrected {r['date']} -> {ref[r['release']]} "
                       f"(release {r['release']} in the {reference} table)")
                r["note"] = "; ".join(x for x in (r["note"], msg) if x)
                r["date"] = ref[r["release"]]
                fixes.append(f"{t.component} {msg}")
            elif r["release"] and r["release"] not in ref:
                log.info("%s release %s is not in the %s table; its date is kept as in the file",
                         t.component, r["release"], reference)
    return fixes


def find_decreases(tables: list[Table]) -> dict[tuple[str, str], list[str]]:
    """{(component, release): [columns that went down vs the previous row]} for bases and entries.
    Rows with n/a values are skipped, so a gap does not hide or fake a decrease."""
    out: dict[tuple[str, str], list[str]] = {}
    for t in tables:
        prev = {"bases": None, "entries": None}
        for r in t.rows:
            for col in ("bases", "entries"):
                v = r[col]
                if v is None:
                    continue
                if prev[col] is not None and v < prev[col]:
                    out.setdefault((t.component, r["release"]), []).append(f"{col} {prev[col]} -> {v}")
                prev[col] = v
    return out


def check_decreases(tables: list[Table]) -> list[str]:
    """Compare observed decreases with the rows the file itself flags as decreases.
    Every observed decrease is logged at WARNING. Returns problems: flagged rows we did NOT observe
    (parser or file problem)."""
    observed = find_decreases(tables)
    flagged = {(t.component, r["release"]) for t in tables for r in t.rows if "decrease" in r["note"].lower()}
    for key, cols in sorted(observed.items(), key=lambda kv: (kv[0][0], float(kv[0][1] or 0))):
        log.warning("%s release %s: decrease in %s (%s in the file)", key[0], key[1], "; ".join(cols),
                    "flagged" if key in flagged else "NOT flagged")
    return [f"{c} release {r} is flagged as a decrease in the file but no decrease was found"
            for c, r in sorted(flagged - set(observed))]


def validate(tables: list[Table], expected: tuple[str, ...]) -> list[str]:
    problems = []
    found = [t.component for t in tables]
    if sorted(found) != sorted(expected):
        problems.append(f"expected components {sorted(expected)}, found {sorted(found)}")
    seen: dict[tuple[str, str], str] = {}
    for t in tables:
        for r in t.rows:
            if not r["release"]:
                continue                        # 'not delivered' rows have no release
            if not r["date"]:
                problems.append(f"{t.component} release {r['release']}: date could not be parsed")
            key = (t.component, r["date"])
            if key in seen:
                problems.append(f"{t.component}: releases {seen[key]} and {r['release']} share date {r['date']}")
            seen[key] = r["release"]
    problems += check_decreases(tables)
    return problems


HEADER_TOTALS_RE = re.compile(r"(\d+)\s+sequences,\s+(\d+)\s+bases,\s+for\s+(traditional|set-based)", re.I)


def cross_check_header(lines: list[str], tables: list[Table]) -> list[str]:
    """Last rows must agree with the two totals printed on the first page of gbrel.txt."""
    tot = {m.group(3).lower(): (int(m.group(1)), int(m.group(2)))
           for l in lines[:40] for m in [HEADER_TOTALS_RE.search(l)] if m}
    last = {t.component: t.rows[-1] for t in tables}
    issues = []
    if "traditional" in tot and "GenBank" in last:
        e, b = tot["traditional"]
        if (last["GenBank"]["entries"], last["GenBank"]["bases"]) != (e, b):
            issues.append(f"GenBank last row != header total ({last['GenBank']['entries']},{last['GenBank']['bases']} vs {e},{b})")
    if "set-based" in tot and all(c in last for c in ("WGS", "TSA", "TLS")):
        e, b = tot["set-based"]
        se = sum(last[c]["entries"] for c in ("WGS", "TSA", "TLS"))
        sb = sum(last[c]["bases"] for c in ("WGS", "TSA", "TLS"))
        if (se, sb) != (e, b):
            issues.append(f"WGS+TSA+TLS last rows != header set-based total (entries {se} vs {e}; bases {sb} vs {b}, diff {b - sb})")
    return issues


# ---------------------------------------------------------------- outputs
def write_csv(tables: list[Table], out: Path) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = [r for t in tables for r in t.rows]
    rows.sort(key=lambda r: (r["component"], r["date"], float(r["release"] or 0)))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def upsert_accession(path: Path, row: dict) -> None:
    """Insert or replace the row keyed by (accession_or_url, release)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    existing = []
    if path.exists():
        with open(path, newline="") as fh:
            existing = list(csv.DictReader(fh, delimiter="\t"))
    existing = [r for r in existing
                if not (r.get("accession_or_url") == row["accession_or_url"] and r.get("release") == row["release"])]
    existing.append(row)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ACC_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(existing)


def inspect(lines: list[str]) -> None:
    for t in find_tables(lines):
        r0, r1 = t.rows[0], t.rows[-1]
        print(f"header at line {t.header_line}: component={t.component!r} title={t.title!r}\n"
              f"    {len(t.rows)} rows, release {r0['release']}..{r1['release']}, {r0['date']}..{r1['date']}")
    print("\nall lines that look like a header:")
    for n, l in enumerate(lines, 1):
        if is_header(l):
            print(f"  {n}: {l.strip()}")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", default=URL)
    ap.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    ap.add_argument("--out-dir", type=Path, default=Path("results/tables"))
    ap.add_argument("--accessions", type=Path, default=Path("config/accessions.tsv"))
    ap.add_argument("--input", type=Path, help="parse this local file instead of downloading (testing)")
    ap.add_argument("--verified-by", default="UNVERIFIED",
                    help="name of the person who checked the download; written to accessions.tsv")
    ap.add_argument("--expect", default=",".join(EXPECTED_COMPONENTS),
                    help="comma-separated components that must all be found")
    ap.add_argument("--inspect", action="store_true", help="print located tables and exit, write nothing")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO, format="%(levelname)s %(message)s")

    if a.input:
        path = a.input
        lines = path.read_text(errors="replace").splitlines()
        release, nbytes, sha = detect_release(lines), path.stat().st_size, sha256_of(path)
    else:
        if a.inspect:
            log.error("--inspect needs --input (it does not download)")
            return 2
        try:
            path, release, nbytes, sha = download(a.url, a.raw_dir)
        except requests.RequestException as exc:
            log.error("download failed: %s", exc)
            log.error("if your network blocks NCBI, download gbrel.txt elsewhere and use --input")
            return 1
        lines = path.read_text(errors="replace").splitlines()

    if a.inspect:
        inspect(lines)
        return 0

    tables = find_tables(lines)
    for fix in reconcile_dates(tables):
        log.warning("%s", fix)
    problems = validate(tables, tuple(x.strip() for x in a.expect.split(",")))
    for issue in cross_check_header(lines, tables):
        log.warning("header cross-check: %s", issue)
    if problems:
        for p in problems:
            log.error(p)
        log.error("run with --input %s --inspect to see what the locator found", path)
        return 3

    out = a.out_dir / f"genbank_growth_gbrel_{release}.csv"
    n = write_csv(tables, out)
    log.info("wrote %d rows from %d tables -> %s", n, len(tables), out)

    if a.verified_by == "UNVERIFIED":
        log.warning("verified_by is UNVERIFIED: re-run with --verified-by <your name> after checking the row")
    upsert_accession(a.accessions, {
        "accession_or_url": a.url, "source": "NCBI GenBank release notes", "release": release,
        "bytes": nbytes, "sha256": sha,
        "download_date": datetime.now(timezone.utc).date().isoformat(),
        "verified_by": a.verified_by,
        "note": f"{len(tables)} tables, {n} rows; parsed CSV committed to git",
    })
    return 0


if __name__ == "__main__":
    sys.exit(main())