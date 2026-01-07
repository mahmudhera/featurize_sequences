import os
import shutil
import subprocess
import pandas as pd

from featurizer_util import make_dir


def run_sei(working_dir, args):
    """
    HumanBase Sei batch API featurization.

    Contract with featurizer.generate_features():
      - working_dir ends with '/'
      - must produce: working_dir + 'sei.csv'
      - should align rows to args['locs']['name']

    Inputs:
      - uses working_dir/locs_sei.fasta if present, else working_dir/locs.fasta
      - requires DeepSea_sei_API.py to exist in the same folder as featurizer.py

    Environment variables (optional):
      - SEI_TOKEN
      - SEI_ASSEMBLY (default hg38)
      - SEI_CHUNK (default 5000)
      - SEI_MAX_INFLIGHT (default 20)
      - SEI_POLL_EVERY (default 10)
      - SEI_TIMEOUT (default 3600)
      - SEI_MAX_RETRIES (default 2)
      - SEI_SKIP_NON4096 (0/1; default 0)
      - SEI_RESUBMIT_ON_ERROR (0/1; default 0)
      - SEI_RESUBMIT_ON_MISMATCH (0/1; default 0)
    """
    # --- working layout
    sei_dir = make_dir(working_dir + 'sei/')
    per_chunk_dir = make_dir(sei_dir + 'per_chunk/')

    # --- choose input fasta
    fasta_in = working_dir + 'locs_sei.fasta'
    if not os.path.exists(fasta_in):
        fasta_in = working_dir + 'locs.fasta'
    if not os.path.exists(fasta_in):
        raise RuntimeError(f"run_sei requires FASTA input: missing {working_dir}locs.fasta (and locs_sei.fasta).")

    # Normalize filename to locs.fasta inside sei_dir so the merged output name is predictable:
    #   locs_sei_merged_feature_predictions.tsv.gz
    fasta_norm = os.path.join(sei_dir, 'locs.fasta')
    shutil.copy2(fasta_in, fasta_norm)

    # --- locate API script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sei_api = os.path.join(script_dir, 'DeepSea_sei_API.py')
    if not os.path.exists(sei_api):
        raise RuntimeError(f"Missing DeepSea_sei_API.py at {sei_api}. "
                           f"Place it in /home/jl2791/scripts/featurizer/.")

    # --- config from env
    token = os.environ.get('SEI_TOKEN', None)
    assembly = os.environ.get('SEI_ASSEMBLY', 'hg38')
    chunk = int(os.environ.get('SEI_CHUNK', '5000'))
    max_inflight = int(os.environ.get('SEI_MAX_INFLIGHT', '20'))
    poll_every = int(os.environ.get('SEI_POLL_EVERY', '10'))
    timeout = int(os.environ.get('SEI_TIMEOUT', '3600'))
    max_retries = int(os.environ.get('SEI_MAX_RETRIES', '2'))

    skip_non4096 = os.environ.get('SEI_SKIP_NON4096', '0') == '1'
    resubmit_on_error = os.environ.get('SEI_RESUBMIT_ON_ERROR', '0') == '1'
    resubmit_on_mismatch = os.environ.get('SEI_RESUBMIT_ON_MISMATCH', '0') == '1'

    state_path = os.path.join(sei_dir, 'sei_batch_state.json')

    cmd = [
        'python3', sei_api,
        fasta_norm,
        '--outdir', per_chunk_dir,
        '--chunk', str(chunk),
        '--assembly', assembly,
        '--poll-every', str(poll_every),
        '--timeout', str(timeout),
        '--state', state_path,
        '--max-inflight', str(max_inflight),
        '--max-retries', str(max_retries),
    ]
    if token:
        cmd += ['--token', token]
    if skip_non4096:
        cmd += ['--skip-non4096']
    if resubmit_on_error:
        cmd += ['--resubmit-on-error']
    if resubmit_on_mismatch:
        cmd += ['--resubmit-on-mismatch']

    # Run in sei_dir so merged output lands there
    subprocess.run(cmd, cwd=sei_dir, check=True)

    merged = os.path.join(sei_dir, 'locs_sei_merged_feature_predictions.tsv.gz')
    if not os.path.exists(merged):
        raise RuntimeError(f"Sei merged output not found: {merged}")

    # TSV.GZ -> sei.csv, align to locs order
    df = pd.read_csv(merged, sep='\t', compression='gzip')

    # Prefer 'sequence_id' column; fall back to first column
    if 'sequence_id' in df.columns:
        df = df.set_index('sequence_id')
    else:
        df = df.set_index(df.columns[0])

    # Align to featurizer's canonical order
    df = df.reindex(index=args['locs']['name'])

    out_csv = working_dir + 'sei.csv'
    df.to_csv(out_csv, float_format='%.5g')
