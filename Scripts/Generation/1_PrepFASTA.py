#!/usr/bin/env python3
"""
1_PrepFASTA.py

Create per-plex .csv files for 2_buildFragFASTA.py and build correct headers for MTPs.

SAAP list is the appended *_MTP.fasta, parent proteins are derived from MaxQuant DP search
txt outputs:
    BP_seq - evidence.txt --> 'Proteins'
    --> accession - proteinGroups.txt 'Fasta headers' (GN=)
    --> gene, description

*_MTP.fasta only has the peptide sequence, but BP is normal MQ identification, so it appears
in the 'Sequence' row; each MTP is matched to its BP by a single-residue difference (same length,
one mismatch) within the indexed evidence seqs, then resolved as above.

STRICT PER-PLEX RESOLUTION. Each *_MTP.fasta is one plex/tissue and is resolved ONLY against that
plex's own DP combined/txt evidence. There is NO crossing: an MTP sequence seen in multiple plexes
is resolved independently in each plex against that plex's evidence, and written into that plex's
own CSV. The FragPipe FASTAs built downstream are therefore plex-specific.

Matching: each MTP file and each DP txt dir collapse to the same token (plex/tissue, lowercased,
non-alphanumerics stripped), e.g. S1_ACGB1_MTP.fasta -> 'acgb1' <-> Ping2018_ACG_B1_DP -> 'acgb1';
S1_Pooled_MTP.fasta -> 'pooled' <-> Bai_2020_Pooled_DP. Only DP dirs are used (Val excluded).
A file matching zero or >1 DP dir is a hard error rather than a silent mis-resolution.

Every MTP whose base peptide resolves to a real accession (within its own plex) is kept
(status='keep'); the rest are 'unresolved'.

Emits one CSV per plex, named '{token}.csv', with the schema 2_buildFragFASTA.py consumes:
    sequence , accession , gene , description , bp_seq , all_accessions , status , n_base_candidates
  - `sequence` is the MTP peptide (2_buildFragFASTA.py keys keep[seq] on it, per plex)

Species routing (ACG/FC/POOLED -> human; else mouse) only selects which root set provides the DP
dirs; it does NOT pool plexes. Each plex still resolves against only its own DP dir.

Example usage:
-------
  python3 1_PrepFASTA.py \
    --mtp-dir     /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-root  /scratch/maropakis.a/MQ_outputs/Ping_2018 \
    --human-root  /scratch/maropakis.a/MQ_outputs/Bai_2020 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Keele_2025 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Tsumagari_2023 \
    --out-dir     /scratch/maropakis.a/Dependencies/mtp_maps/

Token disambiguation: tissue names repeat across studies (Takasugi kidney vs Keele kidney; Keele
cortex vs Tsumagari cortex), so datasets listed in DATASET_SUFFIX get a dataset suffix in their
token (Keele kidney -> 'kidney_keele', Tsumagari cortex rep1 -> 'cortex_1_tsumagari'). Datasets
with globally-unique plex names (Ping ACG/FC, Bai pooled, Takasugi tissues) keep the bare token.
To add a disambiguated dataset, add one entry to DATASET_SUFFIX below.

Each root is walked recursively for every <...>_DP/combined/txt/ dir holding
evidence.txt + proteinGroups.txt. The DP dir bound to each MTP file is decided by token match.
"""

import argparse
import csv
import os
import re
from collections import defaultdict

## header parsers
HDR_RE = re.compile(r'(?:sp|tr)\|([^|]+)\|\S+\s+(.+?)\s+OS=')
GN_RE  = re.compile(r'GN=(\S+)')

## Helper functions
def mtp_token(filename):
    """S1_ACGB1_MTP.fasta -> 'acgb1' ; S9_cortex_keele_MTP.fasta -> 'cortex_keele'.

    The token is the text between the S# sample prefix and _MTP, lowercased, underscores kept so
    disambiguated tissue names (cortex_keele, cortex_1_tsumagari) match dir_token / on-disk names.
    Compact plex names (ACGB1, Pooled) carry no underscore and pass through unchanged.
    """
    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    stem = re.sub(r'^S\d+_', '', stem)
    return stem.lower()

# Datasets whose plex token gets a dataset suffix to avoid tissue-name collisions across studies
# (e.g. Keele cortex vs Tsumagari cortex vs Takasugi has no cortex). Bare datasets need no suffix
# because their plex/tissue names are globally unique. Add a new disambiguated dataset here only.
DATASET_SUFFIX = {'keele2025': 'keele', 'tsumagari_2023': 'tsumagari'}
DATASET_RE = re.compile(r'^([A-Za-z]+(?:_?\d{4}))_(.+)$')   # <Name><year> or <Name>_<year> prefix

def token_from_rest(dataset, rest):
    """Build the canonical plex token from a dataset name + the plex/tissue remainder.

    Bare datasets (Ping/Bai/Takasugi): 'ACG_B1'->'acgb1', 'Pooled'->'pooled', 'Aorta'->'aorta'.
    Suffixed datasets: Keele 'cortex'->'cortex_keele'; Tsumagari 'Cortex_1'->'cortex_1_tsumagari'.
    Matches the names already used in annotations/ and sample_map/ on disk.
    """
    key = dataset.lower()
    if key in DATASET_SUFFIX:
        suf = DATASET_SUFFIX[key]
        # why: Tsumagari keeps the replicate number with underscores (cortex_1_tsumagari);
        # Keele has no replicate (cortex_keele).
        if '_' in rest:
            return '_'.join(rest.lower().split('_') + [suf])
        return f'{rest.lower()}_{suf}'
    return re.sub(r'[^A-Za-z0-9]', '', rest).lower()

def dir_token(txt_dir):
    """.../Ping2018_ACG_B1_DP/combined/txt -> 'acgb1' ; Keele2025_kidney_DP -> 'kidney_keele'."""
    leaf = os.path.basename(os.path.dirname(os.path.dirname(txt_dir)))  # the *_DP folder
    leaf = re.sub(r'_DP$', '', leaf, flags=re.I)
    m = DATASET_RE.match(leaf)
    if not m:
        raise SystemExit(f'cannot parse dataset prefix from DP dir name: {leaf!r}')
    return token_from_rest(m.group(1), m.group(2))

def find_dp_txt_dirs(root):
    """Walk a dataset root; return DP combined/txt dirs (Val excluded) with evidence+proteinGroups."""
    txt_dirs = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if 'evidence.txt' in filenames and 'proteinGroups.txt' in filenames:
            # why: BP identifications come from the DP (dependent-peptide) search only.
            leaf = os.path.basename(os.path.dirname(os.path.dirname(dirpath)))
            if leaf.upper().endswith('_DP'):
                txt_dirs.append(dirpath)
    return sorted(txt_dirs)

def build_token_index(roots, species):
    """Expand dataset roots into {plex/tissue token: DP txt dir}; error on token collision."""
    token_to_dir = {}
    for root in roots:
        if not os.path.isdir(root):
            raise SystemExit(f'{species}: root not found: {root}')
        found = find_dp_txt_dirs(root)
        if not found:
            raise SystemExit(f'{species}: no *_DP/combined/txt under {root}')
        print(f'  {species}: {len(found)} DP txt dir(s) under {root}')
        for d in found:
            tok = dir_token(d)
            if tok in token_to_dir:
                raise SystemExit(f'{species}: duplicate DP token {tok!r}:\n'
                                 f'  {token_to_dir[tok]}\n  {d}')
            token_to_dir[tok] = d
    return token_to_dir

def parse_fasta(path):
    """Yield (header, sequence); header keeps its leading '>'."""
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


## proteinGroups.txt: accession -> (gene, description)
def parse_protein_groups(mq_dir):
    acc_gene, acc_desc = {}, {}
    pg = os.path.join(mq_dir, 'proteinGroups.txt')
    if not os.path.exists(pg):
        raise SystemExit(f'proteinGroups.txt not found: {pg}')
    with open(pg, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'Fasta headers' not in (reader.fieldnames or []):
            raise SystemExit(f"'Fasta headers' missing in {pg}. "
                             f'Available: {reader.fieldnames}')
        for row in reader:
            for entry in (row.get('Fasta headers') or '').split(';'):
                m = HDR_RE.search(entry)
                if not m:
                    continue
                acc_desc.setdefault(m.group(1), m.group(2).strip())
                g = GN_RE.search(entry)
                if g:
                    acc_gene.setdefault(m.group(1), g.group(1))
    return acc_gene, acc_desc


## evidence.txt: Sequence -> 'Proteins' string
def index_evidence(mq_dir):
    seq_to_proteins = {}
    ev = os.path.join(mq_dir, 'evidence.txt')
    if not os.path.exists(ev):
        raise SystemExit(f'evidence.txt not found: {ev}')
    with open(ev, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        cols = reader.fieldnames or []
        if 'Sequence' not in cols or 'Proteins' not in cols:
            raise SystemExit(f"'Sequence'/'Proteins' missing in {ev}. "
                             f'Available: {cols}')
        for row in reader:
            seq = (row.get('Sequence') or '').strip().upper()
            prot = (row.get('Proteins') or '').strip()
            if seq and prot and prot.lower() != 'nan':
                seq_to_proteins.setdefault(seq, prot)
    return seq_to_proteins

def proteins_to_accessions(protein_str):
    """Real UniProt accessions in a MaxQuant 'Proteins' string (drops MTP| entries)."""
    accs = []
    for tok in str(protein_str).split(';'):
        tok = tok.replace('CON__', '')
        if tok.startswith('MTP'):
            continue
        m = re.match(r'(?:sp|tr)\|([^|]+)\|', tok) or re.match(r'^([A-Z0-9][A-Z0-9-]+)$', tok)
        if m:
            accs.append(m.group(1).split('-')[0])
    return list(dict.fromkeys(accs))


## MTP -> base peptide by single-residue difference
def build_length_index(sequences):
    """Group evidence sequences by length for O(matches) base-peptide search."""
    by_len = defaultdict(list)
    for s in sequences:
        by_len[len(s)].append(s)
    return by_len

def find_base_peptides(mtp, by_len):
    """Base candidates = same-length evidence peptides differing at exactly one position."""
    out = []
    for cand in by_len.get(len(mtp), ()):
        if cand == mtp:
            continue
        diff = 0
        for a, b in zip(mtp, cand):
            if a != b:
                diff += 1
                if diff > 1:
                    break
        if diff == 1:
            out.append(cand)
    return out

def resolve_one_mtp(mtp, species, by_len, seq_to_proteins, acc_gene, acc_desc):
    """Resolve a single MTP seq against one plex's DP index; return a record dict."""
    bases = find_base_peptides(mtp, by_len)
    chosen_bp, accs = '', []
    for bp in bases:
        a = proteins_to_accessions(seq_to_proteins.get(bp, ''))
        if a:
            chosen_bp, accs = bp, a
            break
    if accs:
        acc = accs[0]
        return dict(sequence=mtp, species=species, accession=acc,
                    gene=acc_gene.get(acc, ''),
                    description=acc_desc.get(acc, ''),
                    bp_seq=chosen_bp,
                    all_accessions=';'.join(accs),
                    status='keep',
                    n_base_candidates=len(bases))
    return dict(sequence=mtp, species=species, accession='', gene='', description='',
                bp_seq='', all_accessions='', status='unresolved',
                n_base_candidates=len(bases))


## per-file MTP collection (each file is one plex; resolved against its own DP dir only)
def collect_mtp_files(mtp_dir):
    """Read all *_MTP.fasta; return [(filename, token, sorted MTP seqs), ...]."""
    files = []
    seen_tokens = {}
    for fn in sorted(os.listdir(mtp_dir)):
        if not fn.endswith('_MTP.fasta'):
            continue
        tok = mtp_token(fn)
        if tok in seen_tokens:
            # why: two MTP files collapsing to one plex token would silently overwrite one CSV.
            raise SystemExit(f'two MTP files share plex token {tok!r}: '
                             f'{seen_tokens[tok]} and {fn}')
        seen_tokens[tok] = fn
        seqs = set()
        for header, seq in parse_fasta(os.path.join(mtp_dir, fn)):
            if header.startswith('>MTP|'):
                s = seq.strip().upper()
                if s:
                    seqs.add(s)
        files.append((fn, tok, sorted(seqs)))
    return files

def write_plex(token, records, out_dir):
    """One CSV per plex, named '{token}.csv'."""
    path = os.path.join(out_dir, f'{token}.csv')
    cols = ['sequence', 'species', 'accession', 'gene', 'description',
            'bp_seq', 'all_accessions', 'status', 'n_base_candidates']
    counts = {}
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in records:
            counts[r['status']] = counts.get(r['status'], 0) + 1
            w.writerow(r)
    print(f'  {token}: {len(records)} MTP seqs -> {dict(sorted(counts.items()))}  ({path})')


## Run
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True,
                    help='dir of *_MTP.fasta (the appended FASTAs)')
    ap.add_argument('--human-root', action='append', default=[],
                    help='dataset root for human (e.g. .../Ping_2018, .../Bai_2020); '
                         'repeatable. Walked for all *_DP/combined/txt dirs.')
    ap.add_argument('--mouse-root', action='append', default=[],
                    help='dataset root for mouse (e.g. .../Takasugi_2024); repeatable. '
                         'Walked for all *_DP/combined/txt dirs.')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    files = collect_mtp_files(a.mtp_dir)
    print(f'Read {len(files)} *_MTP.fasta files ({len(files)} plexes)')

    # one token -> (DP dir, species) index spanning both root sets. species comes from which root
    # flag the dir was found under, so adding a dataset = passing its root, no token rule to edit.
    token_index = {}
    for species, roots in (('human', a.human_root), ('mouse', a.mouse_root)):
        if not roots:
            continue
        print(f'{species}:')
        for tok, d in build_token_index(roots, species).items():
            if tok in token_index:
                raise SystemExit(f'plex token {tok!r} found under two roots:\n'
                                 f'  {token_index[tok][0]}\n  {d}')
            token_index[tok] = (d, species)

    # resolve each plex STRICTLY against its own matched DP dir; one CSV per plex, no merging
    print('resolving per plex:')
    for fn, tok, seqs in files:
        if tok not in token_index:
            raise SystemExit(f'{fn}: plex token {tok!r} matched no DP dir. '
                             f'Available tokens: {sorted(token_index)}')
        mq_dir, species = token_index[tok]
        leaf = os.path.basename(os.path.dirname(os.path.dirname(mq_dir)))
        acc_gene, acc_desc = parse_protein_groups(mq_dir)
        seq_to_proteins = index_evidence(mq_dir)
        by_len = build_length_index(seq_to_proteins.keys())
        print(f'  {fn} [{species}] <- {tok} ({leaf}): '
              f'{len(acc_desc):,} proteins, {len(seq_to_proteins):,} peptides, {len(seqs)} MTP seqs')
        records = [resolve_one_mtp(m, species, by_len, seq_to_proteins, acc_gene, acc_desc) for m in seqs]
        write_plex(tok, records, a.out_dir)


if __name__ == '__main__':
    main()
