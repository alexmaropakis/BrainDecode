#!/usr/bin/env python3
"""
gen_fragpipe_plex.py

One plex, one call: prepare everything FragPipe needs for a single TMT plex and write a
standalone SLURM submit script you sbatch individually. Adding a new dataset/plex requires NO
edits to any script — you just call this with that plex's parameters on the command line.

For the given plex this does, in order:
  1. annotation.txt   from the plex's sample_map (.xlsx) -- '<channel> <sample>' per line,
                       ordered by TMT channel, repeated labels disambiguated with _<channel>.
  2. raw -> mzML       converts every .raw in RAW_DIR with ThermoRawFileParser (-f=2 plain mzML),
                       skipping any that already exist (idempotent). Skipped if --no-convert.
  3. spectra staging   builds <spectra-root>/<plex>/ with symlinks to the mzML + annotation.txt.
  4. workflow+manifest patches the given .workflow template (db-path, tmtintegrator.channel_num)
                       and writes the .fp-manifest pointing at the staged spectra.
  5. submit script     writes submit_<plex>.sh that runs FragPipe headless for this plex.

Everything plex-specific is a FLAG (mirrors the MQ-XML generator style):
  --plex, --species, --channels, --workflow, --sample-map, plus paths.
The plex token (e.g. 'acgb1', 'pooled', 'aorta') ties this run to its FASTA
(<fasta-dir>/*<plex>*_fragpipe.fasta) and its per-plex CSV from 1_PrepFASTA.py.

Example:
  python gen_fragpipe_plex.py /scratch/maropakis.a/MQ_raw/Ping_2018/ACG/b1 \
    --plex       acgb1 \
    --species    human \
    --channels   10 \
    --workflow   /home/maropakis.a/scripts/FragPipe/templates/TMT10_MS3_Val.workflow \
    --sample-map /scratch/maropakis.a/Dependencies/sample_map/acgb1.xlsx \
    --fasta-dir  /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
    --out-dir    /scratch/maropakis.a/Frag_outputs \
    --spectra-root /scratch/maropakis.a/spectra

Then: sbatch /scratch/maropakis.a/Frag_outputs/submit/submit_acgb1.sh

Run after 2_buildFragFASTA.py (needs the FragPipe FASTAs) and after the sample_map .xlsx exists.
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
    """'Sample name'/'TMT channel' -> 'sample_name'/'tmt_channel'."""
    return re.sub(r'\s+', '_', str(col).strip().lower())

def write_annotation(sample_map, out_path, channels=None):
    """sample_map .xlsx -> FragPipe annotation.txt; return channel count written.

    If `channels` is given, warn on mismatch; otherwise the count is whatever the sample_map holds.
    """
    import pandas as pd
    df = pd.read_excel(sample_map)
    df.columns = [norm_col(c) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        sys.exit(f'{sample_map}: need tmt_channel + sample_name columns, got {list(df.columns)}')
    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    df['tmt_channel'] = df['tmt_channel'].astype(str).str.strip()
    df['sample_name'] = df['sample_name'].astype(str).str.strip()

    # why: a label repeated within a plex (e.g. GIS on 126 and 131) must stay unique downstream.
    counts = Counter(df['sample_name'])
    rows = [(ch, f'{name}_{ch}' if counts[name] > 1 else name)
            for ch, name in zip(df['tmt_channel'], df['sample_name'])]
    rows.sort(key=lambda x: ORD.get(x[0], 999))

    with open(out_path, 'w') as f:
        for ch, name in rows:
            f.write(f'{ch} {name}\n')

    if channels is not None and len(rows) != channels:
        print(f'  WARN: annotation has {len(rows)} channels but --channels {channels}')
    return len(rows)

def convert_raws(raw_dir, trfp):
    """Convert every .raw in raw_dir to plain mzML in place; skip existing. Return mzML paths."""
    mzmls = []
    raws = sorted(glob.glob(os.path.join(raw_dir, '*.raw')))
    if not raws:
        # why: raw_dir may already hold mzML (pre-converted); fall through to collect those.
        mzmls = sorted(glob.glob(os.path.join(raw_dir, '*.mzML')))
        print(f'  no .raw in {raw_dir}; found {len(mzmls)} existing .mzML')
        return mzmls
    for raw in raws:
        base = os.path.splitext(os.path.basename(raw))[0]
        out = os.path.join(raw_dir, f'{base}.mzML')
        if os.path.getsize(out) if os.path.exists(out) else 0:
            print(f'  SKIP convert {base} (mzML exists)')
        else:
            # why: -f=2 = plain indexed mzML; -g would gzip and FragPipe silently skips .mzML.gz.
            subprocess.run([trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
            print(f'  converted {base}')
        if os.path.exists(out):
            mzmls.append(out)
    return mzmls

def stage_spectra(mzmls, annotation_path, spectra_root, plex):
    """Symlink mzML + copy annotation.txt into <spectra-root>/<plex>/. Return staged dir."""
    import shutil
    dst = os.path.join(spectra_root, plex)
    os.makedirs(dst, exist_ok=True)
    for f in mzmls:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    shutil.copy(annotation_path, os.path.join(dst, 'annotation.txt'))
    print(f'  staged {len(mzmls)} mzML + annotation -> {dst}')
    return dst

def find_fasta(fasta_dir, plex):
    """Find the FragPipe FASTA for this plex (built per-plex by 2_buildFragFASTA.py)."""
    # why: names look like S1_ACGB1_fragpipe.fasta or S9_cortex_keele_fragpipe.fasta; the label
    # between S#_ and _fragpipe is the plex token (underscores kept), matched case-insensitively.
    for p in sorted(glob.glob(os.path.join(fasta_dir, '*_fragpipe.fasta'))):
        label = os.path.basename(p)
        label = re.sub(r'^S\d+_', '', label, flags=re.I)
        label = re.sub(r'(?:_MTP)?_fragpipe\.fasta$', '', label, flags=re.I)
        if label.lower() == plex.lower():
            return p
    sys.exit(f'no *_fragpipe.fasta in {fasta_dir} for plex {plex!r}')

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

def write_manifest(mzmls, plex, out_path):
    with open(out_path, 'w') as fh:
        for s in mzmls:
            fh.write(f'{os.path.abspath(s)}\t{plex}\t1\tDDA\n')

SUBMIT_TEMPLATE = """\
#!/usr/bin/env bash
#SBATCH --job-name=fp_{plex}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={ram}G
#SBATCH --time={time}
#SBATCH --output={logdir}/fp_{plex}_%j.out
#SBATCH --error={logdir}/fp_{plex}_%j.err
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

def write_submit(plex, paths, opts, out_path):
    os.makedirs(paths['logdir'], exist_ok=True)
    os.makedirs(paths['workdir'], exist_ok=True)
    text = SUBMIT_TEMPLATE.format(plex=plex, partition=opts.partition, threads=opts.threads,
                                  ram=opts.ram, time=opts.time, logdir=paths['logdir'],
                                  java_home=opts.java_home, fragpipe_bin=opts.fragpipe_bin,
                                  workflow=paths['workflow'], manifest=paths['manifest'],
                                  workdir=paths['workdir'], tools_folder=opts.tools_folder)
    open(out_path, 'w').write(text)
    os.chmod(out_path, 0o755)


def main():
    ap = argparse.ArgumentParser(description='Prep one FragPipe plex end-to-end + submit script.')
    ap.add_argument('raw_dir', help='dir holding this plex\'s .raw (or pre-made .mzML)')
    ap.add_argument('--plex', required=True, help='plex token, e.g. acgb1 / pooled / aorta')
    ap.add_argument('--species', required=True, choices=sorted(SPECIES_TAG))
    ap.add_argument('--channels', type=int, default=None,
                    help='TMT channel count (e.g. 10/16/18). Optional: if omitted, derived from the '
                         'sample_map row count. Pass it only to assert an expected count.')
    ap.add_argument('--workflow', required=True, help='FragPipe .workflow template to patch')
    ap.add_argument('--sample-map', required=True, help='this plex\'s sample_map .xlsx')
    ap.add_argument('--fasta-dir', required=True, help='dir of *_fragpipe.fasta (per-plex)')
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

    plex = a.plex.strip().lower()   # token may contain underscores (cortex_keele); keep them
    if not os.path.isdir(a.raw_dir):
        sys.exit(f'raw_dir not found: {a.raw_dir}')

    wf_dir = os.path.join(a.out_dir, 'workflows')
    mf_dir = os.path.join(a.out_dir, 'manifests')
    annot_dir = os.path.join(a.out_dir, 'annotations')
    submit_dir = os.path.join(a.out_dir, 'submit')
    for d in (wf_dir, mf_dir, annot_dir, submit_dir):
        os.makedirs(d, exist_ok=True)

    print(f'[{plex}] species={a.species} channels={a.channels or "auto"}')

    # 1. annotation (also yields the authoritative channel count)
    annot_path = os.path.join(annot_dir, f'{plex}_annotation.txt')
    channels = write_annotation(a.sample_map, annot_path, a.channels)
    print(f'  {channels} channels')

    # 2. raw -> mzML
    if a.no_convert:
        mzmls = sorted(glob.glob(os.path.join(a.raw_dir, '*.mzML')))
        print(f'  --no-convert: using {len(mzmls)} existing .mzML')
    else:
        mzmls = convert_raws(a.raw_dir, a.trfp)
    if not mzmls:
        sys.exit(f'{plex}: no mzML to stage')

    # 3. stage
    stage_spectra(mzmls, annot_path, a.spectra_root, plex)
    staged = sorted(glob.glob(os.path.join(a.spectra_root, plex, '*.mzML')))

    # 4. workflow + manifest
    fasta = find_fasta(a.fasta_dir, plex)
    wf_path = os.path.join(wf_dir, f'{plex}.workflow')
    mf_path = os.path.join(mf_dir, f'{plex}.fp-manifest')
    write_workflow(a.workflow, fasta, channels, wf_path)
    write_manifest(staged, plex, mf_path)
    print(f'  workflow -> {wf_path}\n  manifest -> {mf_path} ({len(staged)} files)\n  fasta = {os.path.basename(fasta)}')

    # 5. submit script
    submit_path = os.path.join(submit_dir, f'submit_{plex}.sh')
    paths = dict(workflow=wf_path, manifest=mf_path,
                 workdir=os.path.join(a.out_dir, 'results', plex),
                 logdir=os.path.join(a.out_dir, 'logs'))
    write_submit(plex, paths, a, submit_path)
    print(f'  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
