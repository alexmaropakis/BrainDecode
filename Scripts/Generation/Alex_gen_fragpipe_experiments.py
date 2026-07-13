#!/usr/bin/env python3
"""
Alex_gen_fragpipe_experiments.py

One RUN, many EXPERIMENTS, one manifest. This is the experiment-level sibling of
Alex_gen_fragpipe.py (which is one-plex-per-call). Use this when several TMT plexes belong to
the same FragPipe run/manifest -- e.g. a study whose two plexes are two experiments that must be
searched and quantified together (TMT-Integrator groups by experiment across plexes).

The mental model:
  RUN         = the whole FragPipe job (one workflow, one manifest, one submit script, one FASTA).
  EXPERIMENT  = one TMT plex inside that run. You name it (cortex_1, cortex_2, hippocampus-1 ...),
                point it at the folder holding that plex's .raw/.mzML, and give the plex KEY that
                selects that plex's rows from the shared sample_map (its 'TMT plex' or 'sample_ID').

For each experiment this builds its own annotation.txt from the sample_map, converts/stages its
spectra into <spectra-root>/<run>/<experiment>/ (symlinked mzML + that experiment's annotation.txt
-- exactly the per-folder layout FragPipe expects for multi-experiment TMT). Then it writes ONE
manifest listing every file across all experiments (each tagged with its experiment name), patches
ONE workflow (db-path + tmtintegrator.channel_num), and writes ONE submit_<run>.sh.

Multi-plex sample_map layout (one .xlsx for the whole run), e.g. sample_map_tsumagari_cortex.xlsx:
  TMT plex | TMT channel | ParticipantID | Group | MQ | sample_name | sample_ID
     3     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S3        <- experiment cortex_1
     4     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S4        <- experiment cortex_2
The plex KEY you pass with each --experiment is matched against 'TMT plex' (3/4) OR 'sample_ID'
(S3/S4), case-insensitively.

Example (Tsumagari cortex: plex 3 -> cortex_1, plex 4 -> cortex_2, both in one run):
  python Alex_gen_fragpipe_experiments.py \
    --run        cortex_tsumagari \
    --species    mouse \
    --workflow   /home/maropakis.a/scripts/FragPipe/templates/TMT10_MS2_Val.workflow \
    --sample-map /scratch/maropakis.a/Dependencies/sample_map/sample_map_tsumagari_cortex.xlsx \
    --experiment cortex_1 /scratch/maropakis.a/MQ_raw/Tsumagari_2023/cortex_1 3 \
    --experiment cortex_2 /scratch/maropakis.a/MQ_raw/Tsumagari_2023/cortex_2 4 \
    --fasta      /scratch/maropakis.a/Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta \
    --out-dir    /scratch/maropakis.a/Frag_outputs \
    --spectra-root /scratch/maropakis.a/spectra

Then: sbatch /scratch/maropakis.a/Frag_outputs/submit/submit_cortex_tsumagari.sh

You point --fasta straight at the FragPipe FASTA to search against (built by 2_buildFragFASTA.py);
it is patched into the workflow's database.db-path. Run after the sample_map exists.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter

# canonical TMTpro/TMT channel order (covers TMT10 through TMTpro18)
CHANNEL_ORDER = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C',
                 '131', '131N', '131C', '132N', '132C', '133N', '133C', '134N', '134C', '135N']
ORD = {c: i for i, c in enumerate(CHANNEL_ORDER)}

SPECIES_TAG = {
    'human': ('Homo sapiens', 9606),
    'mouse': ('Mus musculus', 10090),
}

# Helpers
def norm_col(col):
    """'Sample name'/'TMT channel'/'sample_ID' -> 'sample_name'/'tmt_channel'/'sample_id'."""
    return re.sub(r'\s+', '_', str(col).strip().lower())

def load_sample_map(path):
    """Read the (possibly multi-plex) sample_map .xlsx once; normalize columns + string cells."""
    import pandas as pd
    df = pd.read_excel(path)
    df.columns = [norm_col(c) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        sys.exit(f'{path}: need tmt_channel + sample_name columns, got {list(df.columns)}')
    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    df['tmt_channel'] = df['tmt_channel'].astype(str).str.strip()
    df['sample_name'] = df['sample_name'].astype(str).str.strip()
    return df

def select_plex(df, key):
    """Return the sample_map rows for one plex, matching KEY against tmt_plex or sample_id.

    KEY is matched case-insensitively against either column so both '3' and 'S3' work. A run with
    a single-plex sample_map (no plex column) may pass key=None to take the whole table.
    """
    if key is None:
        return df
    key = str(key).strip().lower()
    for col in ('tmt_plex', 'sample_id'):
        if col in df.columns:
            sub = df[df[col].astype(str).str.strip().str.lower() == key]
            if len(sub):
                return sub
    have = [c for c in ('tmt_plex', 'sample_id') if c in df.columns] or ['<none>']
    sys.exit(f'sample_map: no rows for plex key {key!r} (searched columns: {", ".join(have)})')

def write_annotation(df, out_path, expected=None):
    """One plex's rows -> FragPipe annotation.txt ('<channel> <sample>' per line). Return count.

    Repeated labels within the plex (e.g. a Bridge on two channels) are disambiguated with _<ch>,
    matching Alex_gen_fragpipe.py. Rows are ordered by canonical TMT channel.
    """
    counts = Counter(df['sample_name'])
    rows = [(ch, f'{name}_{ch}' if counts[name] > 1 else name)
            for ch, name in zip(df['tmt_channel'], df['sample_name'])]
    rows.sort(key=lambda x: ORD.get(x[0], 999))
    with open(out_path, 'w') as f:
        for ch, name in rows:
            f.write(f'{ch} {name}\n')
    if expected is not None and len(rows) != expected:
        print(f'  WARN: {os.path.basename(out_path)} has {len(rows)} channels but expected {expected}')
    return len(rows)

def convert_raws(raw_dir, trfp):
    """Convert every .raw in raw_dir to plain mzML in place; skip existing. Return mzML paths."""
    mzmls = []
    raws = sorted(glob.glob(os.path.join(raw_dir, '*.raw')))
    if not raws:
        # why: raw_dir may already hold mzML (pre-converted); fall through to collect those.
        mzmls = sorted(glob.glob(os.path.join(raw_dir, '*.mzML')))
        print(f'    no .raw in {raw_dir}; found {len(mzmls)} existing .mzML')
        return mzmls
    for raw in raws:
        base = os.path.splitext(os.path.basename(raw))[0]
        out = os.path.join(raw_dir, f'{base}.mzML')
        if os.path.getsize(out) if os.path.exists(out) else 0:
            print(f'    SKIP convert {base} (mzML exists)')
        else:
            # why: -f=2 = plain indexed mzML; -g would gzip and FragPipe silently skips .mzML.gz.
            subprocess.run([trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
            print(f'    converted {base}')
        if os.path.exists(out):
            mzmls.append(out)
    return mzmls

def stage_experiment(mzmls, annotation_path, spectra_root, run, experiment):
    """Symlink one experiment's mzML + copy its annotation.txt into <spectra-root>/<run>/<experiment>/."""
    import shutil
    dst = os.path.join(spectra_root, run, experiment)
    os.makedirs(dst, exist_ok=True)
    for f in mzmls:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    shutil.copy(annotation_path, os.path.join(dst, 'annotation.txt'))
    print(f'    staged {len(mzmls)} mzML + annotation -> {dst}')
    return dst

def patch_line(text, key, value):
    """Replace `key=...` in a workflow, or append it if absent."""
    pat = re.compile(rf'^{re.escape(key)}=.*$', re.MULTILINE)
    line = f'{key}={value}'
    out = pat.sub(line, text)
    return out if out != text else text.rstrip('\n') + f'\n{line}\n'

def write_workflow(template, fasta, channels, out_path):
    wf = open(template).read()
    wf = patch_line(wf, 'database.db-path', os.path.abspath(fasta))
    wf = patch_line(wf, 'tmtintegrator.channel_num', str(channels))
    open(out_path, 'w').write(wf)

def write_manifest(entries, out_path):
    """FragPipe manifest across all experiments: <mzML>\t<experiment>\t1\tDDA."""
    with open(out_path, 'w') as fh:
        for path, experiment in entries:
            fh.write(f'{os.path.abspath(path)}\t{experiment}\t1\tDDA\n')

def write_sample_table(entries, out_path):
    """Human-readable record of file -> experiment (not used by FragPipe)."""
    with open(out_path, 'w') as fh:
        fh.write('file\texperiment\n')
        for path, experiment in entries:
            fh.write(f'{os.path.basename(path)}\t{experiment}\n')

SUBMIT_TEMPLATE = """\
#!/usr/bin/env bash
#SBATCH --job-name=fp_{run}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={ram}G
#SBATCH --time={time}
#SBATCH --output={logdir}/fp_{run}_%j.out
#SBATCH --error={logdir}/fp_{run}_%j.err
set -euo pipefail
export JAVA_HOME={java_home}
export PATH=$JAVA_HOME/bin:$PATH

{fragpipe_bin} --headless \\
  --workflow {workflow} \\
  --manifest {manifest} \\
  --workdir  {workdir} \\
  --threads  {threads} \\
  --ram      {ram} \\
  --config-tools-folder {tools_folder}
"""

def write_submit(run, paths, opts, out_path):
    os.makedirs(paths['logdir'], exist_ok=True)
    os.makedirs(paths['workdir'], exist_ok=True)
    text = SUBMIT_TEMPLATE.format(run=run, partition=opts.partition, threads=opts.threads,
                                  ram=opts.ram, time=opts.time, logdir=paths['logdir'],
                                  java_home=opts.java_home, fragpipe_bin=opts.fragpipe_bin,
                                  workflow=paths['workflow'], manifest=paths['manifest'],
                                  workdir=paths['workdir'], tools_folder=opts.tools_folder)
    open(out_path, 'w').write(text)
    os.chmod(out_path, 0o755)


def main():
    ap = argparse.ArgumentParser(
        description='Prep one FragPipe RUN spanning multiple TMT experiments (plexes) -> one manifest + submit script.')
    ap.add_argument('--run', required=True, help='run token, e.g. cortex_tsumagari (names workflow/manifest/submit + spectra subdir)')
    ap.add_argument('--species', required=True, choices=sorted(SPECIES_TAG))
    ap.add_argument('--experiment', required=True, action='append', nargs=3,
                    metavar=('NAME', 'RAW_DIR', 'PLEX'),
                    help='one TMT experiment: its name, the dir with its .raw/.mzML, and the plex '
                         'key selecting its rows from the sample_map (its TMT plex or sample_ID, '
                         'e.g. 3 or S3). Repeat --experiment for each plex in the run.')
    ap.add_argument('--sample-map', required=True, help='the run\'s multi-plex sample_map .xlsx')
    ap.add_argument('--workflow', required=True, help='FragPipe .workflow template to patch')
    ap.add_argument('--fasta', required=True, help='FragPipe FASTA file to search against (patched into the workflow)')
    ap.add_argument('--out-dir', required=True, help='Frag_outputs root (workflows/manifests/...)')
    ap.add_argument('--spectra-root', required=True)
    ap.add_argument('--trfp', default=os.path.expanduser('~/thermoRawFileParser/ThermoRawFileParser'))
    ap.add_argument('--no-convert', action='store_true', help='skip raw->mzML (spectra already mzML)')
    # submit-script knobs
    ap.add_argument('--fragpipe-bin', default='/home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe')
    ap.add_argument('--tools-folder', default='/home/maropakis.a/fragpipe/fragpipe-24.0/tools')
    ap.add_argument('--java-home', default=os.path.expanduser('~/bin/jdk-17.0.18+8'))
    ap.add_argument('--partition', default='short')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--ram', type=int, default=64)
    ap.add_argument('--time', default='24:00:00')
    a = ap.parse_args()

    run = a.run.strip().lower()   # token may contain underscores (cortex_tsumagari); keep them
    if not os.path.isfile(a.fasta):
        sys.exit(f'--fasta not found: {a.fasta}')

    wf_dir = os.path.join(a.out_dir, 'workflows')
    mf_dir = os.path.join(a.out_dir, 'manifests')
    annot_dir = os.path.join(a.out_dir, 'annotations')
    submit_dir = os.path.join(a.out_dir, 'submit')
    for d in (wf_dir, mf_dir, annot_dir, submit_dir):
        os.makedirs(d, exist_ok=True)

    print(f'[{run}] species={a.species}  {len(a.experiment)} experiments')

    smap = load_sample_map(a.sample_map)

    # Build every experiment: annotation from its plex rows, then convert + stage its spectra.
    entries = []           # (mzML_path, experiment_name) across all experiments -> single manifest
    channel_counts = {}    # experiment -> channel count (must agree; one workflow, one channel_num)
    seen = set()
    for name, raw_dir, plex_key in a.experiment:
        name = name.strip()
        if name in seen:
            sys.exit(f'duplicate experiment name {name!r}')
        seen.add(name)
        if not os.path.isdir(raw_dir):
            sys.exit(f'[{name}] raw_dir not found: {raw_dir}')
        print(f'  experiment {name!r}  plex={plex_key}  raw_dir={raw_dir}')

        sub = select_plex(smap, plex_key)
        annot_path = os.path.join(annot_dir, f'{run}_{name}_annotation.txt')
        channel_counts[name] = write_annotation(sub, annot_path)
        print(f'    {channel_counts[name]} channels')

        if a.no_convert:
            mzmls = sorted(glob.glob(os.path.join(raw_dir, '*.mzML')))
            print(f'    --no-convert: using {len(mzmls)} existing .mzML')
        else:
            mzmls = convert_raws(raw_dir, a.trfp)
        if not mzmls:
            sys.exit(f'[{name}]: no mzML to stage')

        staged_dir = stage_experiment(mzmls, annot_path, a.spectra_root, run, name)
        staged = sorted(glob.glob(os.path.join(staged_dir, '*.mzML')))
        entries.extend((s, name) for s in staged)

    # One workflow means one TMT channel count -- all experiments must agree.
    uniq = set(channel_counts.values())
    if len(uniq) > 1:
        print(f'  WARN: experiments disagree on channel count {channel_counts}; using max')
    channels = max(uniq)

    # One FASTA, one workflow, one manifest, one submit for the whole run.
    fasta = a.fasta
    wf_path = os.path.join(wf_dir, f'{run}.workflow')
    mf_path = os.path.join(mf_dir, f'{run}.fp-manifest')
    tbl_path = os.path.join(annot_dir, f'{run}_samples.txt')
    write_workflow(a.workflow, fasta, channels, wf_path)
    write_manifest(entries, mf_path)
    write_sample_table(entries, tbl_path)
    print(f'  channel_num = {channels}\n  workflow -> {wf_path}'
          f'\n  manifest -> {mf_path} ({len(entries)} files, {len(seen)} experiments)'
          f'\n  samples  -> {tbl_path}\n  fasta = {os.path.basename(fasta)}')

    submit_path = os.path.join(submit_dir, f'submit_{run}.sh')
    paths = dict(workflow=wf_path, manifest=mf_path,
                 workdir=os.path.join(a.out_dir, 'results', run),
                 logdir=os.path.join(a.out_dir, 'logs'))
    write_submit(run, paths, a, submit_path)
    print(f'  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
