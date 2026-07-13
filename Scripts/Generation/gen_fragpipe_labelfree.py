#!/usr/bin/env python3
"""
gen_fragpipe_labelfree.py

One dataset, one call: prepare everything FragPipe needs for a single label-free DIA run and
write a standalone SLURM submit script you sbatch individually. Adding a new dataset requires NO
edits to any script -- you just call this with that dataset's parameters on the command line.

This is the label-free DIA sibling of gen_fragpipe_plex.py (TMT). The differences that matter:
  - no TMT: no channels, no annotation.txt, no tmtintegrator.channel_num patching.
  - the manifest is built from a metadata.csv (ID,Condition,...): experiment=Condition,
    bioreplicate=unique-per-subject, datatype=DIA. That is exactly the layout FragPipe's
    MSstats export (diann.generate-msstats=true) wants for a two-group LFQ comparison --
    each file becomes its own quant column, grouped by Condition.
  - the search enzyme (Trypsin/LysC/GluC) is patched into MSFragger per run, since a dataset
    is typically split into one protease per subdir (PD_2026/{Trypsin,LysC,GluC}/).

For the given dataset this does, in order:
  1. metadata          reads metadata.csv, matches each raw/mzML to a subject ID, and derives
                       (experiment=Condition, bioreplicate) for every file.
  2. raw -> mzML       converts every .raw in RAW_DIR with ThermoRawFileParser (-f=2 plain mzML),
                       skipping any that already exist (idempotent). Skipped if --no-convert.
  3. spectra staging   builds <spectra-root>/<dataset>/ with symlinks to the mzML.
  4. workflow+manifest patches the given .workflow template (db-path, and --enzyme if given)
                       and writes the .fp-manifest pointing at the staged spectra.
  5. submit script     writes submit_<dataset>.sh that runs FragPipe headless for this set.

Example:
  python gen_fragpipe_labelfree.py /scratch/maropakis.a/MQ_raw/PD_2026/Trypsin \
    --dataset    pd_trypsin \
    --enzyme     trypsin \
    --metadata   /scratch/maropakis.a/MQ_raw/PD_2026/metadata.csv \
    --workflow   /home/maropakis.a/scripts/FragPipe/templates/DIA_SpecLib_Quant.workflow \
    --fasta-dir  /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
    --out-dir    /scratch/maropakis.a/Frag_outputs \
    --spectra-root /scratch/maropakis.a/spectra

Then: sbatch /scratch/maropakis.a/Frag_outputs/submit/submit_pd_trypsin.sh

Run after 2_buildFragFASTA.py (needs the FragPipe FASTAs) and once metadata.csv exists.
"""

import argparse
import glob
import os
import re
import subprocess
import sys

SPECIES_TAG = {
    'human': ('Homo sapiens', 9606),
    'mouse': ('Mus musculus', 10090),
}

# MSFragger enzyme presets, keyed by --enzyme. Each patches the primary-enzyme block of the
# workflow. Names match FragPipe's enzyme dropdown; strict* variants ignore the proline rule.
ENZYME_SPECS = {
    'trypsin': dict(name='stricttrypsin', cut='KR', nocut='', sense='C'),
    'lysc':    dict(name='lysc',          cut='K',  nocut='', sense='C'),
    'gluc':    dict(name='gluc',          cut='DE', nocut='', sense='C'),
}

# Helpers
def norm_col(col):
    """'Sample name'/'Condition' -> 'sample_name'/'condition'."""
    return re.sub(r'\s+', '_', str(col).strip().lower())

def load_metadata(path):
    """metadata.csv -> ordered list of (id, condition) rows. Requires ID + Condition columns."""
    import pandas as pd
    df = pd.read_csv(path)
    df.columns = [norm_col(c) for c in df.columns]
    if not {'id', 'condition'} <= set(df.columns):
        sys.exit(f'{path}: need ID + Condition columns, got {list(df.columns)}')
    df = df.dropna(subset=['id', 'condition'])
    df['id'] = df['id'].astype(str).str.strip()
    df['condition'] = df['condition'].astype(str).str.strip()
    return list(zip(df['id'], df['condition']))

def match_id(stem, ids):
    """Match an mzML filename stem to a metadata ID (whole-token match; longest wins).

    IDs are short tokens (e.g. subject numbers) that appear in the filename; \b keeps 2264 from
    matching a bare 264 and vice-versa. Exact stem == id is preferred over an embedded match.
    """
    if stem in ids:
        return stem
    hits = [i for i in ids if re.search(rf'(?<![0-9A-Za-z]){re.escape(i)}(?![0-9A-Za-z])', stem)]
    return max(hits, key=len) if hits else None

def build_manifest_rows(mzmls, meta):
    """Assign (path, experiment=condition, bioreplicate) to each mzML from metadata.

    bioreplicate is a globally-unique integer per file (1..N over files sorted by name), so
    FragPipe/MSstats treats every file as a distinct biological replicate within its Condition.
    Files with no metadata match are dropped with a loud warning (never silently included).
    """
    cond_by_id = dict(meta)
    ids = list(cond_by_id)
    rows, unmatched = [], []
    biorep = 0
    for f in sorted(mzmls, key=os.path.basename):
        stem = os.path.splitext(os.path.basename(f))[0]
        sid = match_id(stem, ids)
        if sid is None:
            unmatched.append(os.path.basename(f))
            continue
        biorep += 1
        rows.append((f, sid, cond_by_id[sid], biorep))
    if unmatched:
        print(f'  WARN: {len(unmatched)} mzML had no metadata match and were DROPPED: '
              + ', '.join(unmatched))
    return rows

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

def stage_spectra(mzmls, spectra_root, dataset):
    """Symlink mzML into <spectra-root>/<dataset>/. Return staged dir."""
    dst = os.path.join(spectra_root, dataset)
    os.makedirs(dst, exist_ok=True)
    for f in mzmls:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    print(f'  staged {len(mzmls)} mzML -> {dst}')
    return dst

def find_fasta(fasta_dir, dataset):
    """Find the FragPipe FASTA for this dataset (built per-dataset by 2_buildFragFASTA.py)."""
    # why: names look like S1_PD_fragpipe.fasta or S9_cortex_keele_fragpipe.fasta; the label
    # between S#_ and _fragpipe is the dataset token (underscores kept), matched case-insensitively.
    for p in sorted(glob.glob(os.path.join(fasta_dir, '*_fragpipe.fasta'))):
        label = os.path.basename(p)
        label = re.sub(r'^S\d+_', '', label, flags=re.I)
        label = re.sub(r'(?:_MTP)?_fragpipe\.fasta$', '', label, flags=re.I)
        if label.lower() == dataset.lower():
            return p
    sys.exit(f'no *_fragpipe.fasta in {fasta_dir} for dataset {dataset!r}')

def patch_line(text, key, value):
    """Replace `key=...` in a workflow, or append it if absent."""
    pat = re.compile(rf'^{re.escape(key)}=.*$', re.MULTILINE)
    line = f'{key}={value}'
    out = pat.sub(line, text)
    return out if out != text else text.rstrip('\n') + f'\n{line}\n'

def write_workflow(template, fasta, enzyme, out_path):
    wf = open(template).read()
    wf = patch_line(wf, 'database.db-path', os.path.abspath(fasta))
    if enzyme is not None:
        spec = ENZYME_SPECS[enzyme]
        wf = patch_line(wf, 'msfragger.search_enzyme_name_1', spec['name'])
        wf = patch_line(wf, 'msfragger.search_enzyme_cut_1', spec['cut'])
        wf = patch_line(wf, 'msfragger.search_enzyme_nocut_1', spec['nocut'])
        wf = patch_line(wf, 'msfragger.search_enzyme_sense_1', spec['sense'])
        wf = patch_line(wf, 'msfragger.misc.fragger.enzyme-dropdown-1', spec['name'])
    open(out_path, 'w').write(wf)

def write_manifest(rows, out_path):
    """FragPipe manifest: <mzML>\t<experiment=condition>\t<bioreplicate>\tDIA."""
    with open(out_path, 'w') as fh:
        for path, _sid, cond, biorep in rows:
            fh.write(f'{os.path.abspath(path)}\t{cond}\t{biorep}\tDIA\n')

def write_sample_table(rows, out_path):
    """Human-readable record of file -> id / condition / bioreplicate (not used by FragPipe)."""
    with open(out_path, 'w') as fh:
        fh.write('file\tid\tcondition\tbioreplicate\n')
        for path, sid, cond, biorep in rows:
            fh.write(f'{os.path.basename(path)}\t{sid}\t{cond}\t{biorep}\n')

SUBMIT_TEMPLATE = """\
#!/usr/bin/env bash
#SBATCH --job-name=fp_{dataset}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={ram}G
#SBATCH --time={time}
#SBATCH --output={logdir}/fp_{dataset}_%j.out
#SBATCH --error={logdir}/fp_{dataset}_%j.err
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

def write_submit(dataset, paths, opts, out_path):
    os.makedirs(paths['logdir'], exist_ok=True)
    os.makedirs(paths['workdir'], exist_ok=True)
    text = SUBMIT_TEMPLATE.format(dataset=dataset, partition=opts.partition, threads=opts.threads,
                                  ram=opts.ram, time=opts.time, logdir=paths['logdir'],
                                  java_home=opts.java_home, fragpipe_bin=opts.fragpipe_bin,
                                  workflow=paths['workflow'], manifest=paths['manifest'],
                                  workdir=paths['workdir'], tools_folder=opts.tools_folder)
    open(out_path, 'w').write(text)
    os.chmod(out_path, 0o755)


def main():
    ap = argparse.ArgumentParser(description='Prep one label-free DIA dataset end-to-end + submit script.')
    ap.add_argument('raw_dir', help='dir holding this dataset\'s .raw (or pre-made .mzML)')
    ap.add_argument('--dataset', required=True, help='dataset token, e.g. pd_trypsin / pd_lysc')
    ap.add_argument('--species', default='human', choices=sorted(SPECIES_TAG),
                    help='kept for record/parity with the TMT generator; does not alter the workflow')
    ap.add_argument('--enzyme', choices=sorted(ENZYME_SPECS), default=None,
                    help='protease to patch into MSFragger (trypsin/lysc/gluc). Omit to keep the '
                         'template\'s enzyme.')
    ap.add_argument('--metadata', required=True,
                    help='metadata.csv with ID + Condition columns; drives experiment/bioreplicate')
    ap.add_argument('--workflow', required=True, help='FragPipe DIA .workflow template to patch')
    ap.add_argument('--fasta-dir', required=True, help='dir of *_fragpipe.fasta (per-dataset)')
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

    dataset = a.dataset.strip().lower()   # token may contain underscores (pd_trypsin); keep them
    if not os.path.isdir(a.raw_dir):
        sys.exit(f'raw_dir not found: {a.raw_dir}')

    wf_dir = os.path.join(a.out_dir, 'workflows')
    mf_dir = os.path.join(a.out_dir, 'manifests')
    annot_dir = os.path.join(a.out_dir, 'annotations')
    submit_dir = os.path.join(a.out_dir, 'submit')
    for d in (wf_dir, mf_dir, annot_dir, submit_dir):
        os.makedirs(d, exist_ok=True)

    print(f'[{dataset}] species={a.species} enzyme={a.enzyme or "template-default"}')

    # 1. metadata (drives experiment=Condition, bioreplicate)
    meta = load_metadata(a.metadata)
    print(f'  {len(meta)} subjects in metadata')

    # 2. raw -> mzML
    if a.no_convert:
        mzmls = sorted(glob.glob(os.path.join(a.raw_dir, '*.mzML')))
        print(f'  --no-convert: using {len(mzmls)} existing .mzML')
    else:
        mzmls = convert_raws(a.raw_dir, a.trfp)
    if not mzmls:
        sys.exit(f'{dataset}: no mzML to stage')

    # 3. stage
    staged_dir = stage_spectra(mzmls, a.spectra_root, dataset)
    staged = sorted(glob.glob(os.path.join(staged_dir, '*.mzML')))

    # 4. manifest rows from metadata (matched against the staged files)
    rows = build_manifest_rows(staged, meta)
    if not rows:
        sys.exit(f'{dataset}: no mzML matched metadata; nothing to run')
    conds = sorted({c for _p, _s, c, _b in rows})
    print(f'  {len(rows)} files matched -> conditions: {", ".join(conds)}')

    # 5. workflow + manifest + sample table
    fasta = find_fasta(a.fasta_dir, dataset)
    wf_path = os.path.join(wf_dir, f'{dataset}.workflow')
    mf_path = os.path.join(mf_dir, f'{dataset}.fp-manifest')
    tbl_path = os.path.join(annot_dir, f'{dataset}_samples.txt')
    write_workflow(a.workflow, fasta, a.enzyme, wf_path)
    write_manifest(rows, mf_path)
    write_sample_table(rows, tbl_path)
    print(f'  workflow -> {wf_path}\n  manifest -> {mf_path} ({len(rows)} files)'
          f'\n  samples  -> {tbl_path}\n  fasta = {os.path.basename(fasta)}')

    # 6. submit script
    submit_path = os.path.join(submit_dir, f'submit_{dataset}.sh')
    paths = dict(workflow=wf_path, manifest=mf_path,
                 workdir=os.path.join(a.out_dir, 'results', dataset),
                 logdir=os.path.join(a.out_dir, 'logs'))
    write_submit(dataset, paths, a, submit_path)
    print(f'  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
