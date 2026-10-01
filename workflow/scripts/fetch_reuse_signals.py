"""Collect reuse signals for the tiering model (task 7).

Two archives, joined on the BioProject accession:
  1. ENA Portal API  -> every read run for a taxon (default: Escherichia coli, tax_tree(562)).
     Streamed to disk, then stored as Parquet (runs) and aggregated to one row per study.
  2. Europe PMC REST -> number of papers that mention each study accession, and the years
     they were published. Text-mined accessions (ACCESSION_ID) and the literal accession
     string are both searched, because each misses cases the other catches.

Outputs
  <out>/ena_runs.parquet        one row per run (raw metadata, typed)
  <out>/studies.parquet         one row per study (features + reuse label)
  <out>/epmc_cache.jsonl        one line per queried accession (resumable cache)
  <out>/metadata_issues.tsv     metadata problems found while joining (documented, not hidden)

Everything is resumable: the run table is only downloaded if missing, and studies already in
the Europe PMC cache are not queried again.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_COUNT = "https://www.ebi.ac.uk/ena/portal/api/count"
EPMC_SEARCH ="https://www.ebi.ac.uk/europepmc/webservices/rest/search"
RUN_FIELDS = ["run_accession", "study_accession", "secondary_study_accession", "sample_accession",
              "first_public", "instrument_platform", "instrument_model", "library_strategy",
              "library_source", "library_selection", "library_layout", "base_count", "read_count",
              "fastq_bytes", "submitted_bytes", "center_name"]
UA = {"User-Agent": "aitu-bioinfo-course-project/1.0 (storage apocalypse; tiering)"}
log = logging.getLogger("fetch_reuse_signals")


# ---------------------------------------------------------------- ENA
def ena_count(query: str) -> int:
    r = requests.get(ENA_COUNT, params={"result": "read_run", "query": query}, timeout=120, headers=UA)
    r.raise_for_status()
    return int(r.text.split()[-1])


def download_chunk(query: str, dest: Path, limit: int = 0, retries: int = 4) -> int:
    """Stream one query to dest; return data rows written. Retried, and checked by the caller."""
    params = {"result": "read_run", "query": query, "fields": ",".join(RUN_FIELDS),
              "format": "tsv", "limit": limit}
    for attempt in range(1, retries + 1):
        try:
            with requests.post(ENA_SEARCH, data=params, stream=True, timeout=600, headers=UA) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            with open(dest, "rb") as fh:
                return sum(1 for _ in fh) - 1
        except requests.RequestException as exc:
            log.warning("chunk %s attempt %d failed: %s", query, attempt, exc)
            time.sleep(5 * attempt)
    raise RuntimeError(f"ENA download failed for {query}")


def download_runs(query: str, dest: Path, limit: int = 0) -> None:
    """Download the ENA run table one first_public year at a time.

    A single 600k-row stream was silently cut short by the server on the first attempt (383k rows
    arrived, no error raised). So each year is fetched separately and its row count is compared
    with ENA's own /count for the same query; a short chunk is retried, never accepted."""
    t0 = time.time()
    parts_dir = dest.parent / "ena_chunks"
    parts_dir.mkdir(exist_ok=True)
    if limit:
        download_chunk(query, dest, limit)
        return
    chunks = [(f"{query} AND first_public<2008-01-01", "pre2008")] + [
        (f"{query} AND first_public>={y}-01-01 AND first_public<={y}-12-31", str(y))
        for y in range(2008, time.gmtime().tm_year + 1)]
    chunks.append((f"{query} AND first_public>{time.gmtime().tm_year}-12-31", "future"))
    files = []
    for q, tag in chunks:
        part = parts_dir / f"runs_{tag}.tsv"
        expected = ena_count(q)
        got = -1
        for attempt in range(3):
            if part.exists():
                with open(part, "rb") as fh:
                    got = sum(1 for _ in fh) - 1
            if got == expected:
                break
            got = download_chunk(q, part)
        if got != expected:
            raise RuntimeError(f"{tag}: ENA says {expected} runs, received {got}")
        log.info("  %s: %d runs", tag, got)
        files.append(part)
    total = ena_count(query)
    n_chunks = sum(sum(1 for _ in open(f, "rb")) - 1 for f in files)
    if n_chunks != total:
        log.warning("year chunks hold %d runs but the whole query counts %d (runs without a usable "
                    "first_public date); difference recorded, not fatal", n_chunks, total)
    with open(dest, "wb") as out:
        for k, f in enumerate(files):
            with open(f, "rb") as fh:
                header = fh.readline()
                if k == 0:
                    out.write(header)
                for line in fh:
                    out.write(line)
    log.info("ENA run table: %d chunks, %.1f MB in %.0f s", len(files), dest.stat().st_size / 1e6,
             time.time() - t0)


def sum_semicolon(s: pd.Series) -> pd.Series:
    """fastq_bytes is '123;456' for paired runs -> 579. Empty -> NaN."""
    parts = s.fillna("").astype(str).str.split(";", expand=True)
    return parts.apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=1)


def load_runs(tsv: Path) -> tuple[pd.DataFrame, list[dict]]:
    df = pd.read_csv(tsv, sep="\t", dtype=str, low_memory=False)
    issues = []
    df["first_public"] = pd.to_datetime(df["first_public"], errors="coerce")
    for col in ("base_count", "read_count"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["fastq_bytes"] = sum_semicolon(df["fastq_bytes"])
    df["submitted_bytes"] = sum_semicolon(df["submitted_bytes"])
    # storage footprint: archived FASTQ if ENA made one, otherwise what the submitter uploaded
    df["bytes"] = df["fastq_bytes"].fillna(df["submitted_bytes"])
    for what, mask in [("no study_accession", df["study_accession"].isna()),
                       ("no first_public date", df["first_public"].isna()),
                       ("no byte size (neither fastq nor submitted)", df["bytes"].isna()),
                       ("base_count = 0 (no reads extracted, e.g. native ONT/PacBio upload)", df["base_count"].eq(0))]:
        issues.append({"issue": what, "count": int(mask.sum()), "unit": "runs", "share": round(float(mask.mean()), 5)})
    return df, issues


def aggregate_studies(runs: pd.DataFrame) -> pd.DataFrame:
    r = runs.dropna(subset=["study_accession", "first_public"]).copy()
    r["is_illumina"] = r["instrument_platform"].eq("ILLUMINA")
    r["is_ont"] = r["instrument_platform"].eq("OXFORD_NANOPORE")
    r["is_pacbio"] = r["instrument_platform"].eq("PACBIO_SMRT")
    r["is_paired"] = r["library_layout"].eq("PAIRED")
    mode = lambda s: s.mode().iat[0] if s.notna().any() else "unknown"  # noqa: E731
    g = r.groupby("study_accession")
    st = pd.DataFrame({
        "first_public": g["first_public"].min(),
        "last_public": g["first_public"].max(),
        "n_runs": g.size(),
        "n_samples": g["sample_accession"].nunique(),
        "bytes": g["bytes"].sum(min_count=1),
        "bases": g["base_count"].sum(min_count=1),
        "share_illumina": g["is_illumina"].mean(),
        "share_ont": g["is_ont"].mean(),
        "share_pacbio": g["is_pacbio"].mean(),
        "share_paired": g["is_paired"].mean(),
        "library_strategy": g["library_strategy"].agg(mode),
        "library_source": g["library_source"].agg(mode),
        "instrument_model": g["instrument_model"].agg(mode),
        "center_name": g["center_name"].agg(mode),
        "secondary_study_accession": g["secondary_study_accession"].first(),
    })
    st.index.name = "study_accession"
    return st.reset_index()


# ---------------------------------------------------------------- Europe PMC
def epmc_query(acc: str, session: requests.Session, retries: int = 4) -> dict:
    q = f'(ACCESSION_ID:"{acc}" OR "{acc}")'
    params = {"query": q, "format": "json", "pageSize": 1000, "resultType": "lite"}
    for attempt in range(retries):
        try:
            r = session.get(EPMC_SEARCH, params=params, timeout=60)
            if r.status_code == 429 or r.status_code >= 500:
                raise requests.HTTPError(f"status {r.status_code}")
            r.raise_for_status()
            d = r.json()
            years = sorted(int(x["pubYear"]) for x in d["resultList"]["result"] if x.get("pubYear"))
            return {"acc": acc, "hits": int(d.get("hitCount", 0)), "years": years}
        except (requests.RequestException, ValueError, KeyError) as exc:
            if attempt == retries - 1:
                return {"acc": acc, "hits": None, "years": [], "error": str(exc)}
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def query_all(accs: list[str], cache: Path, threads: int) -> dict[str, dict]:
    done: dict[str, dict] = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            rec = json.loads(line)
            if rec.get("hits") is not None:
                done[rec["acc"]] = rec
    todo = [a for a in accs if a not in done]
    log.info("Europe PMC: %d cached, %d to query with %d threads", len(done), len(todo), threads)
    lock, t0 = threading.Lock(), time.time()
    local = threading.local()

    def worker(acc):
        if not hasattr(local, "s"):
            local.s = requests.Session()
            local.s.headers.update(UA)
        return epmc_query(acc, local.s)

    with open(cache, "a") as fh, ThreadPoolExecutor(threads) as ex:
        futs = [ex.submit(worker, a) for a in todo]
        for i, f in enumerate(as_completed(futs), 1):
            rec = f.result()
            with lock:
                fh.write(json.dumps(rec) + "\n")
                if rec.get("hits") is not None:
                    done[rec["acc"]] = rec
            if i % 1000 == 0:
                fh.flush()
                log.info("  %d/%d (%.1f q/s)", i, len(todo), i / (time.time() - t0))
    failed = [a for a in accs if a not in done]
    if failed:
        log.warning("%d accessions failed after retries; they are dropped, not counted as zero", len(failed))
    return done


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", default="tax_tree(562)", help="ENA portal query for read runs")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--public-from", default="2010-01-01")
    ap.add_argument("--public-to", default="2021-12-31",
                    help="only studies first public by this date get queried: they need years to be reused")
    ap.add_argument("--max-studies", type=int, default=0, help="0 = all (use a small number for tests)")
    ap.add_argument("--run-limit", type=int, default=0, help="ENA row limit, 0 = all (tests only)")
    ap.add_argument("--runs-tsv", type=Path, help="use this local run table instead of downloading")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed-cache", type=Path, help="pre-filled Europe PMC cache (offline test run)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    a.out.mkdir(parents=True, exist_ok=True)

    cache = a.out / "epmc_cache.jsonl"
    if a.seed_cache and not cache.exists():
        cache.write_bytes(a.seed_cache.read_bytes())
    tsv = a.runs_tsv or a.out / "ena_runs.tsv"
    if not tsv.exists():
        download_runs(a.query, tsv, a.run_limit)
    runs, issues = load_runs(tsv)
    runs.to_parquet(a.out / "ena_runs.parquet", index=False)
    log.info("runs: %d rows, %.2f PB of FASTQ", len(runs), runs["bytes"].sum() / 1e15)

    st = aggregate_studies(runs)
    window = st["first_public"].between(pd.Timestamp(a.public_from), pd.Timestamp(a.public_to))
    sel = st[window].copy()
    if sel.empty:
        log.error("no study first public between %s and %s", a.public_from, a.public_to)
        return 1
    if a.max_studies:
        sel = sel.sample(n=min(a.max_studies, len(sel)), random_state=1)
    issues.append({"issue": "studies outside the date window (not labelled)", "count": int((~window).sum()), "unit": "studies",
                   "share": round(float((~window).mean()), 5)})
    # BioProjects are what papers cite; ERP/SRP secondary accessions are also searched below
    acc_list = sorted(set(sel["study_accession"]) | set(sel["secondary_study_accession"].dropna()))
    res = query_all(acc_list, a.out / "epmc_cache.jsonl", a.threads)

    def combine(row):
        recs = [res.get(row["study_accession"]), res.get(row["secondary_study_accession"])]
        recs = [r for r in recs if r]
        if not recs:
            return pd.Series({"n_papers": None, "first_paper_year": None, "n_papers_later": None})
        years = max((r["years"] for r in recs), key=len)  # the two searches overlap; take the richer one
        hits = max(r["hits"] for r in recs)
        y0 = years[0] if years else None
        later = sum(1 for y in years if y0 is not None and y > y0)
        return pd.Series({"n_papers": hits, "first_paper_year": y0, "n_papers_later": later})

    sel = pd.concat([sel.reset_index(drop=True), sel.reset_index(drop=True).apply(combine, axis=1)], axis=1)
    sel = sel.dropna(subset=["n_papers"])
    sel["reused"] = (sel["n_papers"] >= 2).astype(int)
    sel.to_parquet(a.out / "studies.parquet", index=False)
    pd.DataFrame(issues).to_csv(a.out / "metadata_issues.tsv", sep="\t", index=False)
    log.info("studies labelled: %d; cited>=1: %.1f%%; reused (>=2 papers): %.1f%%", len(sel),
             100 * (sel["n_papers"] >= 1).mean(), 100 * sel["reused"].mean())
    return 0


if __name__ == "__main__":
    sys.exit(main())
