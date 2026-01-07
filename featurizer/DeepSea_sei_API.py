#!/usr/bin/env python3
"""
DeepSea/Sei batch client (with sequence_id merge + sanity checks + resubmit):
- Chunk FASTA and submit to HumanBase Sei API
- Poll to completion with windowed inflight jobs
- Download per-chunk feature_predictions.tsv.gz
- Sanity-check: verify TSV sequence IDs match chunk FASTA headers; resubmit on mismatch (optional)
- Merge into a single TSV.GZ in the working directory, with a 'sequence_id' column
- Resume-safe via a state JSON

API endpoints:
  POST https://hb.flatironinstitute.org/api/deepsea/jobs/
  GET  https://hb.flatironinstitute.org/api/deepsea/jobs/{job_id}/

USAGE: 
python3 DeepSea_sei_API.py locs.fasta  \
    --outdir deepsea_sei   --chunk 1000 \
    --assembly hg38   --submit-interval 5 \
    --max-inflight 20   --poll-every 15 \
    --timeout 21600   --queue-timeout 7200 \
    --on-queue-timeout wait \
    --resubmit-on-error \
    --resubmit-on-mismatch \ 
    --max-retries 3 
"""

import argparse
import os
import sys
import time
import json
import gzip
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import requests
from requests.adapters import HTTPAdapter, Retry

# =========================== CONFIG DEFAULTS ===========================

BASE_URL = "https://hb.flatironinstitute.org"
JOBS_ENDPOINT = "/api/deepsea/jobs/"       # POST submit; GET {job_id}/ for status
USER_AGENT = "sei-batch-client/1.5"

DEFAULT_CHUNK = 5000           # sequences per job
DEFAULT_ASSEMBLY = "hg19"      # or "hg38"
DEFAULT_SUBMIT_INTERVAL = 2.0  # seconds between submissions
DEFAULT_TIMEOUT = 60 * 60      # per-job timeout (sec) used only for "running" state
DEFAULT_POLL = 10              # seconds between poll cycles

DEFAULT_MAX_INFLIGHT = 20      # max concurrently submitted (not yet downloaded) jobs
DEFAULT_QUEUE_TIMEOUT = 0      # if >0, max seconds to tolerate 'queued' before action
DEFAULT_ON_QUEUE_TIMEOUT = "wait"  # wait | skip | requeue

BACKOFF_FACTOR = 1.5
BACKOFF_CAP = 60               # cap per-job poll backoff (seconds)

# ======================================================================

def mk_session(token=None):
    s = requests.Session()
    s.headers.update({"User-Agent": USER_AGENT, "Accept": "*/*"})
    if token:
        s.headers.update({"Authorization": f"Bearer {token}"})
    retry = Retry(
        total=8,
        connect=5,
        read=5,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
        raise_on_status=False,
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s

# ---------------- FASTA utilities ----------------

def read_fasta(path):
    """Yield (header_without_>, sequence_no_newlines)."""
    hdr, seq = None, []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)

def read_fasta_headers(path):
    return [h for h, _ in read_fasta(path)]

def chunk_records(records, chunk_size):
    batch = []
    for rec in records:
        batch.append(rec)
        if len(batch) >= chunk_size:
            yield batch
            batch = []
    if batch:
        yield batch

def fasta_text(records):
    return "\n".join(f">{h}\n{s}" for h, s in records)

# -------------- TSV sequence-id detection ----------------

def detect_seqid_column(header_fields):
    """
    Heuristics to find the sequence identifier column from a TSV header.
    Returns column index (0-based) or None.
    """
    keys = [h.strip().lower() for h in header_fields]
    candidates = [
        "sequence", "sequence_id", "seqname", "seq_name",
        "name", "id", "header"
    ]
    for c in candidates:
        if c in keys:
            return keys.index(c)
    # fallback: prefer first non-empty column name
    return 0 if header_fields else None

def tsv_sequence_ids(tsv_gz_path):
    """Return ordered list of IDs from the detected sequence ID column."""
    ids = []
    with gzip.open(tsv_gz_path, "rt") as fh:
        header = fh.readline()
        if not header:
            return ids
        header_fields = header.rstrip("\n").split("\t")
        col = detect_seqid_column(header_fields)
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if col is not None and col < len(fields):
                ids.append(fields[col])
    return ids

def validate_tsv_matches_fasta(fasta_headers_for_chunk, tsv_gz_path):
    """
    Strict equality check of N and order between FASTA headers (for this chunk) and TSV sequence IDs.
    Returns (ok: bool, reason: str, details: dict)
    """
    t_ids = tsv_sequence_ids(tsv_gz_path)
    details = {"fasta_n": len(fasta_headers_for_chunk), "tsv_n": len(t_ids)}
    if len(fasta_headers_for_chunk) != len(t_ids):
        return (False, f"count_mismatch({len(fasta_headers_for_chunk)} fasta vs {len(t_ids)} tsv)", details)
    for i, (a, b) in enumerate(zip(fasta_headers_for_chunk, t_ids), start=1):
        if a != b:
            details.update({"first_mismatch_index": i, "fasta_id": a, "tsv_id": b})
            return (False, "order_or_name_mismatch", details)
    return (True, "ok", details)

# ---------------- HTTP helpers ----------------

def submit_job(session, upload_text, title, assembly):
    url = BASE_URL + JOBS_ENDPOINT
    payload = {
        "ds_model": "sei",
        "input_type": "fasta",
        "title": title,
        "upload_textarea": upload_text,
        "genome_assembly": assembly,
    }
    r = session.post(url, json=payload, timeout=120)
    if not r.ok:
        raise RuntimeError(f"Submit HTTP {r.status_code}:\n{r.text}")
    js = r.json()
    job_id = js.get("job_id")
    if not job_id:
        raise RuntimeError(f"Submit response missing job_id: {js}")
    return job_id

def get_job(session, job_id):
    url = f"{BASE_URL}{JOBS_ENDPOINT}{job_id}/"
    r = session.get(url, timeout=60)
    if not r.ok:
        raise RuntimeError(f"Status HTTP {r.status_code} for {job_id}:\n{r.text}")
    return r.json()

def parse_iso(ts):
    try:
        dt = datetime.fromisoformat(ts)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    except Exception:
        return None

def adaptive_sleep(base, factor, cap, step):
    return min(cap, base * (factor ** max(0, step - 1)))

def _wget(url, outfile):
    """Try wget with resume & retries. Return True on success."""
    if shutil.which("wget") is None:
        return False
    if url.startswith("http://s3-") or url.startswith("http://s3."):
        url = "https://" + url[len("http://"):]
    args = [
        "wget", "-q", "-c",
        "--retry-connrefused",
        "--waitretry=5",
        "--read-timeout=30",
        "--tries=10",
        "-O", outfile,
        url
    ]
    return subprocess.call(args) == 0 and os.path.exists(outfile) and os.path.getsize(outfile) > 0

def download_stream(url, outfile):
    if url.startswith("http://s3-") or url.startswith("http://s3."):
        url = "https://" + url[len("http://"):]
    tmp = outfile + ".part"
    with requests.get(url, stream=True, timeout=600, headers={"User-Agent": USER_AGENT}) as r:
        r.raise_for_status()
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                if chunk:
                    f.write(chunk)
    os.replace(tmp, outfile)

def fetch_predictions(js, outfile):
    files = js.get("files") or {}
    feat = files.get("feature_predictions")
    if not feat or "url" not in feat:
        raise RuntimeError(f"No feature_predictions URL in: {js}")
    url = feat["url"]
    if not _wget(url, outfile):
        download_stream(url, outfile)
    print(f"[download] {outfile}")

# ---------------- Merge with sequence_id ----------------

def count_data_rows_gz(tsv_gz_path):
    """Count data rows (excluding header) in a gzipped TSV."""
    with gzip.open(tsv_gz_path, "rt") as fin:
        header = fin.readline()
        if not header:
            return 0
        return sum(1 for _ in fin)

def merge_results_with_seqid(tsv_files, merged_outfile, fasta_path, chunk_size):
    """
    Concatenate gzipped TSVs in order, writing header once.
    Prepend a 'sequence_id' column using FASTA headers mapped by chunk and row index.
    """
    headers = [h for h, _ in read_fasta(fasta_path)]
    total_headers = len(headers)

    header_slices = []
    start = 0
    for f in sorted(tsv_files):
        n_rows = count_data_rows_gz(f)
        end = min(start + n_rows, total_headers)
        header_slices.append(headers[start:end])
        start = end

    with gzip.open(merged_outfile, "wt") as fout:
        wrote_header = False
        for idx, f in enumerate(sorted(tsv_files)):
            seqids = header_slices[idx] if idx < len(header_slices) else []
            with gzip.open(f, "rt") as fin:
                for line_no, line in enumerate(fin):
                    if line_no == 0:
                        if not wrote_header:
                            line = line.rstrip("\n")
                            fout.write("sequence_id\t" + line + "\n")
                            wrote_header = True
                        continue
                    row_idx = line_no - 1
                    try:
                        sid = seqids[row_idx]
                    except IndexError:
                        raise RuntimeError(
                            f"Row/header mismatch while merging {f}: data row {row_idx+1} "
                            "has no matching FASTA header."
                        )
                    fout.write(f"{sid}\t{line}")
    print(f"[merged] {merged_outfile} (with sequence_id)")

# ---------------- State helpers ----------------

def load_state(path):
    if path and os.path.exists(path):
        with open(path, "r") as f:
            return json.load(f)
    return {}

def save_state(path, state):
    if not path:
        return
    tmp = path + ".part"
    with open(tmp, "w") as f:
        json.dump(state, f, indent=2)
    os.replace(tmp, path)

# ---------------- Main runner ----------------

def run(fasta_path, outdir, chunk_size, assembly, token,
        submit_interval, poll_every, timeout, state_path,
        max_inflight, queue_timeout, on_queue_timeout,
        resubmit_on_error, resubmit_on_mismatch, max_retries,
        skip_non4096):

    session = mk_session(token)
    fasta_path = os.path.expanduser(fasta_path)
    outdir = os.path.expanduser(outdir)
    os.makedirs(outdir, exist_ok=True)
    stem = Path(fasta_path).stem

    state = load_state(state_path)  # str(chunk_idx) -> {job_id,status,outfile,retries,...}

    # Load and optionally filter sequences
    records_all = list(read_fasta(fasta_path))
    if len(records_all) == 0:
        print("No sequences in FASTA.", file=sys.stderr)
        return

    if skip_non4096:
        kept = [(h, s) for (h, s) in records_all if len(s) == 4096]
        dropped = len(records_all) - len(kept)
        if dropped > 0:
            print(f"[filter] skipping {dropped} sequences not length 4096")
        records = kept
    else:
        records = records_all

    total_seq = len(records)
    if total_seq == 0:
        print("No sequences left to process after filtering.", file=sys.stderr)
        return

    # Build chunk list
    chunks = []
    idx = 0
    while idx < total_seq:
        end = min(idx + chunk_size, total_seq)
        chunks.append((len(chunks) + 1, records[idx:end]))  # (chunk_idx, batch)
        idx = end

    print(f"[plan] {len(chunks)} chunk(s), {total_seq} sequences total")

    inflight = {}  # chunk_idx -> {job_id, submitted_at, backoff_step}
    completed_outfiles = []
    next_to_submit = 0

    backoff_base = max(2, int(poll_every))
    backoff_cap = BACKOFF_CAP

    def submit_for_chunk(chunk_idx, batch, title_suffix=""):
        outfile = os.path.join(outdir, f"{stem}.chunk{chunk_idx:04d}_feature_predictions.tsv.gz")

        # If already downloaded, mark complete
        if os.path.exists(outfile):
            print(f"[skip] chunk {chunk_idx} already exists")
            entry = state.get(str(chunk_idx), {})
            entry.update({"status": "completed", "outfile": outfile})
            state[str(chunk_idx)] = entry
            save_state(state_path, state)
            completed_outfiles.append(outfile)
            return False

        entry = state.get(str(chunk_idx), {})
        job_id = entry.get("job_id")
        if job_id:
            print(f"[resume] chunk {chunk_idx} -> job_id={job_id}")
        else:
            title = f"Sei batch {stem} [{chunk_idx}]"
            if title_suffix:
                title += f" {title_suffix}"
            job_id = submit_job(session, fasta_text(batch), title, assembly)
            retries = entry.get("retries", 0)
            print(f"[submit] chunk {chunk_idx} ({len(batch)} seq) -> job_id={job_id} (retries={retries})")
            entry.update({"job_id": job_id, "status": "submitted", "outfile": outfile, "retries": entry.get("retries", 0)})
            state[str(chunk_idx)] = entry
            save_state(state_path, state)
            time.sleep(submit_interval)

        inflight[chunk_idx] = {"job_id": job_id, "submitted_at": time.time(), "backoff_step": 1}
        return True

    def resubmit_chunk(chunk_idx, reason):
        entry = state.get(str(chunk_idx), {})
        retries = int(entry.get("retries", 0))
        if retries >= max_retries:
            print(f"[resubmit] chunk {chunk_idx} NOT resubmitted (retries={retries} >= max_retries={max_retries}); reason={reason}")
            return False
        # Prepare batch
        batch = chunks[chunk_idx - 1][1]
        # Remember previous job_id
        prev = entry.get("job_id")
        if prev:
            entry.setdefault("job_id_prev", []).append(prev)
        entry["retries"] = retries + 1
        state[str(chunk_idx)] = entry
        save_state(state_path, state)
        print(f"[resubmit] chunk {chunk_idx} (retry {retries+1}/{max_retries}) due to {reason}")
        # Submit replacement
        submit_for_chunk(chunk_idx, batch, title_suffix=f"(retry {retries+1})")
        return True

    # Prime window
    while next_to_submit < len(chunks) and len(inflight) < max_inflight:
        ci, batch = chunks[next_to_submit]
        submit_for_chunk(ci, batch)
        next_to_submit += 1

    # Process windows
    while inflight or next_to_submit < len(chunks):
        for ci in list(inflight.keys()):
            info = inflight.get(ci)
            if not info:
                continue
            job_id = info["job_id"]
            entry = state.get(str(ci), {})
            outfile = entry.get("outfile") or os.path.join(outdir, f"{stem}.chunk{ci:04d}_feature_predictions.tsv.gz")

            # Poll
            try:
                js = get_job(session, job_id)
                status = (js.get("status") or "").lower()
            except Exception as e:
                # transient error → backoff
                sleep_s = adaptive_sleep(backoff_base, BACKOFF_FACTOR, backoff_cap, info["backoff_step"])
                info["backoff_step"] += 1
                time.sleep(sleep_s)
                continue

            # Queue guard
            if status == "queued" and queue_timeout > 0:
                created = parse_iso(js.get("created", ""))
                waited = (datetime.now(timezone.utc) - created).total_seconds() if created else (time.time() - info["submitted_at"])
                if waited >= queue_timeout:
                    if on_queue_timeout == "skip":
                        print(f"[skip-queued] chunk {ci} job {job_id} queued {int(waited)}s >= {queue_timeout}s")
                        inflight.pop(ci, None)
                        if next_to_submit < len(chunks):
                            cnext, bnext = chunks[next_to_submit]; next_to_submit += 1
                            submit_for_chunk(cnext, bnext)
                        continue
                    elif on_queue_timeout == "requeue":
                        print(f"[requeue] chunk {ci} job {job_id} queued {int(waited)}s; submitting replacement")
                        resubmit_chunk(ci, reason="queue_timeout")
                        inflight.pop(ci, None)
                        continue
                    # else "wait": fall through

            if status == "completed":
                try:
                    fetch_predictions(js, outfile)
                except Exception as e:
                    print(f"[download-error] chunk {ci} job {job_id}: {e}")
                    # Let it repoll next loop (could be transient)
                    sleep_s = adaptive_sleep(backoff_base, BACKOFF_FACTOR, backoff_cap, info["backoff_step"])
                    info["backoff_step"] = min(info["backoff_step"] + 1, 10)
                    time.sleep(sleep_s)
                    continue

                # Post-download sanity check (TSV vs FASTA headers for this chunk)
                fasta_headers_for_chunk = [h for h, _ in chunks[ci - 1][1]]
                ok, reason, details = validate_tsv_matches_fasta(fasta_headers_for_chunk, outfile)
                if ok:
                    completed_outfiles.append(outfile)
                    entry.update({"status": "completed", "outfile": outfile})
                    state[str(ci)] = entry
                    save_state(state_path, state)
                    inflight.pop(ci, None)
                    # free a slot
                    if next_to_submit < len(chunks):
                        cnext, bnext = chunks[next_to_submit]; next_to_submit += 1
                        submit_for_chunk(cnext, bnext)
                    time.sleep(0.5)
                else:
                    print(f"[mismatch] chunk {ci} {reason} details={details}")
                    if resubmit_on_mismatch and resubmit_chunk(ci, reason=f"mismatch:{reason}"):
                        inflight.pop(ci, None)  # replacement will get submitted
                    else:
                        # keep as completed_mismatch and move on
                        entry.update({"status": f"completed_mismatch:{reason}", "outfile": outfile})
                        state[str(ci)] = entry
                        save_state(state_path, state)
                        inflight.pop(ci, None)
                        if next_to_submit < len(chunks):
                            cnext, bnext = chunks[next_to_submit]; next_to_submit += 1
                            submit_for_chunk(cnext, bnext)
                    time.sleep(0.5)

            elif status.startswith("error") or status == "failed":
                print(f"[failed] chunk {ci} job {job_id}: {js.get('status')}")
                if resubmit_on_error and resubmit_chunk(ci, reason="error"):
                    inflight.pop(ci, None)
                else:
                    inflight.pop(ci, None)
                    entry.update({"status": "failed", "error": js.get("status", "error")})
                    state[str(ci)] = entry
                    save_state(state_path, state)
                    if next_to_submit < len(chunks):
                        cnext, bnext = chunks[next_to_submit]; next_to_submit += 1
                        submit_for_chunk(cnext, bnext)
            else:
                # running or queued: adaptive sleep
                sleep_s = adaptive_sleep(backoff_base, BACKOFF_FACTOR, backoff_cap, info["backoff_step"])
                info["backoff_step"] = min(info["backoff_step"] + 1, 10)
                time.sleep(sleep_s)

        # If inflight drained but there are more chunks, refill window
        if not inflight and next_to_submit < len(chunks):
            while next_to_submit < len(chunks) and len(inflight) < max_inflight:
                ci, batch = chunks[next_to_submit]
                submit_for_chunk(ci, batch)
                next_to_submit += 1

    # Merge whatever we downloaded (with sequence_id)
    if completed_outfiles:
        merged = os.path.join(os.getcwd(), f"{stem}_sei_merged_feature_predictions.tsv.gz")
        merge_results_with_seqid(sorted(completed_outfiles), merged, fasta_path=fasta_path, chunk_size=chunk_size)
        print(f"[done] merged {len(completed_outfiles)} chunk(s) -> {merged}")
    else:
        print("[done] no completed outputs to merge")

# ---------------- CLI ----------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Submit FASTA to HumanBase Sei in batches, download & merge predictions with sequence_id (resume-safe)."
    )
    p.add_argument("fasta", help="Input FASTA file")
    p.add_argument("-o", "--outdir", default="sei_results", help="Directory for per-chunk outputs (default: sei_results)")
    p.add_argument("-c", "--chunk", type=int, default=DEFAULT_CHUNK, help="Sequences per job (default: 5000)")
    p.add_argument("-a", "--assembly", default=DEFAULT_ASSEMBLY, help='Genome assembly (default: "hg19")')
    p.add_argument("-t", "--token", default=None, help="Bearer token if required")

    p.add_argument("--submit-interval", type=float, default=DEFAULT_SUBMIT_INTERVAL,
                   help=f"Seconds to wait between submissions (default: {DEFAULT_SUBMIT_INTERVAL})")
    p.add_argument("--poll-every", type=int, default=DEFAULT_POLL,
                   help=f"Polling base interval seconds (default: {DEFAULT_POLL})")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                   help=f"Per-job running timeout seconds (default: {DEFAULT_TIMEOUT})")
    p.add_argument("--state", default="sei_batch_state.json", help="Path to state JSON for resume (default: sei_batch_state.json)")

    p.add_argument("--max-inflight", type=int, default=DEFAULT_MAX_INFLIGHT,
                   help=f"Max jobs submitted but not yet downloaded (default: {DEFAULT_MAX_INFLIGHT})")
    p.add_argument("--queue-timeout", type=int, default=DEFAULT_QUEUE_TIMEOUT,
                   help=f"If >0, max seconds to tolerate 'queued' before action (default: {DEFAULT_QUEUE_TIMEOUT}=no limit)")
    p.add_argument("--on-queue-timeout", choices=["wait", "skip", "requeue"], default=DEFAULT_ON_QUEUE_TIMEOUT,
                   help=f"Action when 'queued' exceeds queue-timeout (default: {DEFAULT_ON_QUEUE_TIMEOUT})")

    # New: resubmission and checks
    p.add_argument("--resubmit-on-error", action="store_true",
                   help="Resubmit a chunk if job status is error/failed")
    p.add_argument("--resubmit-on-mismatch", action="store_true",
                   help="Resubmit a chunk if downloaded TSV sequence IDs do not match chunk FASTA headers")
    p.add_argument("--max-retries", type=int, default=2,
                   help="Maximum resubmissions per chunk (default: 2)")
    p.add_argument("--skip-non4096", action="store_true",
                   help="Skip sequences whose length != 4096 before submission")

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    try:
        if args.chunk < 50:
            print(f"[note] chunk size {args.chunk} is very small; consider >= 500 to reduce server load.", file=sys.stderr)
        run(
            fasta_path=args.fasta,
            outdir=args.outdir,
            chunk_size=args.chunk,
            assembly=args.assembly,
            token=args.token,
            submit_interval=args.submit_interval,
            poll_every=args.poll_every,
            timeout=args.timeout,
            state_path=args.state,
            max_inflight=args.max_inflight,
            queue_timeout=args.queue_timeout,
            on_queue_timeout=args.on_queue_timeout,
            resubmit_on_error=args.resubmit_on_error,
            resubmit_on_mismatch=args.resubmit_on_mismatch,
            max_retries=args.max_retries,
            skip_non4096=args.skip_non4096,
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
