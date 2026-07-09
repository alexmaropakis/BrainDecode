#!/usr/bin/env python3
"""
2_buildFragFASTA.py

Build FragPipe-ready FASTAs from *_MTP.fasta using the per-plex CSVs from 1_PrepFASTA.py.

STRICT PER-PLEX. Each *_MTP.fasta (one plex/tissue) is built using ONLY its own plex CSV
('{token}.csv' in --csv-dir, e.g. acgb1.csv, pooled.csv, aorta.csv). An MTP sequence kept in one
plex is NOT carried into another plex's FASTA; the keep-set is loaded per plex, never merged.

Reference entries pass through unchanged. Kept MTP entries get a Philosopher-safe mock-UniProt
header whose accession is <acc>-<mtp_id>-<plex_token> (made unique per entry). Reversed rev_ decoys
are appended for *every* target, including MTPs.

Species (for the OS/OX header tag) is read from the per-plex CSV's 'species' column written by
1_PrepFASTA.py — never guessed from the filename, so a new dataset needs no edit here.

Example:
  python3 2_buildFragFASTA.py \
    --mtp-dir /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --csv-dir /scratch/maropakis.a/Dependencies/mtp_maps/ \
    --out-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe/

Run after 1_PrepFASTA.py
"""

import argparse
import csv
import os
import re

SPECIES_TAG = {
    'human': ('Homo sapiens', 9606),
    'mouse': ('Mus musculus', 10090),
}
MTP_ID_RE = re.compile(r'MTP\|(\d+)')   # real header format is '>MTP|7998_0_base...'
DECOY_PREFIX = 'rev_'

## Helper functions
def plex_token(filename):
    """e.g. 'S5_ACGB5_MTP.fasta' -> 'acgb5' ; 'S9_cortex_keele_MTP.fasta' -> 'cortex_keele'.

    The token is the text between the S# sample prefix and _MTP, lowercased, underscores kept.
    Must match 1_PrepFASTA.mtp_token so '{token}.csv' is the right per-plex CSV.
    """
    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    stem = re.sub(r'^S\d+_', '', stem)
    return stem.lower()


def load_keep(csv_path):
    """Return ({sequence: (accession, gene)}, species) for MTPs marked 'keep' in one plex CSV."""
    keep, species = {}, None
    with open(csv_path, newline='') as fh:
        for row in csv.DictReader(fh):
            # why: species is written per row by 1_PrepFASTA (uniform within a plex); take it from
            # any row so the OS/OX header tag is never guessed from the filename.
            species = species or row.get('species')
            if row['status'] == 'keep':
                keep[row['sequence']] = (row['accession'], row['gene'])
    return keep, species

def parse_fasta(path):
    """Yield (header, sequence) tuples; header keeps its leading '>'."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        yield header, ''.join(seq)

def mtp_header(accession, gene, mid, species):
    os_name, ox = SPECIES_TAG[species]
    return (f'>sp|{accession}|{gene}-mut {gene} mistranslated {mid} '
            f'OS={os_name} OX={ox} GN={gene} PE=1 SV=1')

def unique_accession(base, seen):
    """Return base, or base-d2, base-d3, ... if already taken."""
    accession, k = base, 2
    while accession in seen:
        accession = f'{base}-d{k}'
        k += 1
    seen.add(accession)
    return accession

def build(src, dst, keep, species, tok):
    targets, seen_acc = [], set()
    n_ref = n_mtp = 0

    for header, seq in parse_fasta(src):
        if header.startswith('>MTP|'):
            if seq not in keep:
                continue
            acc, gene = keep[seq]
            m = MTP_ID_RE.search(header)
            # why: header is '>MTP|7998_0_base...'; take the numeric id as MTP7998. A bare pipe in
            # the accession breaks Philosopher header parsing (silent entry drop), so sanitize.
            mid = f'MTP{m.group(1)}' if m else re.sub(r'[^A-Za-z0-9]', '_', header[1:].split()[0])
            accession = unique_accession(f'{acc}-{mid}-{tok}', seen_acc)
            targets.append((mtp_header(accession, gene, mid, species), seq))
            n_mtp += 1
        else:
            targets.append((header, seq))
            n_ref += 1

    with open(dst, 'w') as out:
        # Forward targets.
        for header, seq in targets:
            out.write(f'{header}\n{seq}\n')
        # Reversed decoys for every target (references AND MTPs).
        for header, seq in targets:
            out.write(f'>{DECOY_PREFIX}{header[1:]}\n{seq[::-1]}\n')

    return n_ref, n_mtp

## Run
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True)
    ap.add_argument('--csv-dir', required=True,
                    help='dir of per-plex CSVs from 1_PrepFASTA.py ({token}.csv)')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)

    n = 0
    for fn in sorted(os.listdir(a.mtp_dir)):
        if not fn.endswith('_MTP.fasta'):
            continue
        tok = plex_token(fn)
        csv_path = os.path.join(a.csv_dir, f'{tok}.csv')
        if not os.path.exists(csv_path):
            # why: each plex must have its own CSV; a missing one would silently drop all its MTPs.
            raise SystemExit(f'{fn}: per-plex CSV not found: {csv_path}')
        keep, sp = load_keep(csv_path)
        if sp not in SPECIES_TAG:
            raise SystemExit(f'{tok}.csv: species {sp!r} not in {sorted(SPECIES_TAG)}')
        dst = os.path.join(a.out_dir, fn.replace('_MTP.fasta', '_fragpipe.fasta'))
        n_ref, n_mtp = build(os.path.join(a.mtp_dir, fn), dst, keep, sp, tok)
        print(f'{fn} [{sp}] <- {tok}.csv: {n_ref} ref + {n_mtp} MTP kept')
        n += 1

    print(f'Wrote {n} FASTAs')


if __name__ == '__main__':
    main()
