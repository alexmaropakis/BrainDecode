#!/usr/bin/env python3

"""
A script to map SAAP sequences to parent proteins via BLASTp, per species. 
Reads MTP_seq, BP_seq, and Dataset columns of a SAAP_quant_df.csv.
ACG*/FC* datasets are human, other tissues are mouse. 

Writes 4 CSVs (filtered drops trypsin / immunoglobulin / reference / no-parent):
  {prefix}_human_unfiltered.csv   {prefix}_human_filtered.csv
  {prefix}_mouse_unfiltered.csv   {prefix}_mouse_filtered.csv

The single-amino-acid swap is computed from MTP_seq vs BP_seq and stored in the
output (swap = MTP_aa->BP_aa, plus position / bp_aa / mtp_aa columns).
Each unique (MTP_seq, Dataset) is blasted once; output rows are per unique
(MTP_seq, BP_seq, Dataset) so distinct swaps on the same MTP are preserved.

Requires NCBI BLAST+ (makeblastdb, blastp) on PATH.

Example usage:
python3 BLAST_mtps.py \
  --csv        /scratch/maropakis.a/Dependencies/SAAP_quant_df.csv \
  --human-ref  /scratch/maropakis.a/Dependencies/FASTA/HUMAN.fasta \
  --mouse-ref  /scratch/maropakis.a/Dependencies/FASTA/MOUSE_UP000000589_10090.fasta \
  --out-dir    /scratch/maropakis.a/Dependencies/mtp_maps/ \
  --prefix     April26 \
  --threads 16

"""
import argparse, csv, os, re, subprocess

IG_GENE    = re.compile(r'^(IG[HKL]|JCHAIN|IGJ)', re.I)
TRYP_GENE  = re.compile(r'^(PRSS[123]|TRY\d*)$', re.I)
IG_TITLE   = re.compile(r'immunoglobulin', re.I)
TRYP_TITLE = re.compile(r'trypsin', re.I)
GN_RE      = re.compile(r'\bGN=(\S+)')

def species_for(dataset):
    # Human Datasets are tokens beginning with ACG or FC (e.g. ACGB1, FCB3);
    # everything else (Aorta, Brain, Heart, ...) is mouse.
    # Adjust this if your Dataset values use a different naming scheme.
    tokens = re.split(r'[._\-\s]', dataset.upper())
    is_human = any(t.startswith('ACG') or t.startswith('FC') for t in tokens)
    return 'human' if is_human else 'mouse'

def dataset_slug(dataset):
    return re.sub(r'[^A-Za-z0-9]', '', dataset).lower()

def compute_swap(mtp, bp):
    # Single-AA swap between MTP and BP. Direction reported: BP -> MTP.
    if not bp:
        return dict(pos='', bp_aa='', mtp_aa='', swap='no_bp')
    if len(mtp) != len(bp):
        return dict(pos='', bp_aa='', mtp_aa='', swap='len_mismatch')
    diffs = [(i, m, b) for i, (m, b) in enumerate(zip(mtp, bp)) if m != b]
    if not diffs:
        return dict(pos='', bp_aa='', mtp_aa='', swap='identical')
    if len(diffs) == 1:
        i, m, b = diffs[0]
        return dict(pos=i + 1, bp_aa=b, mtp_aa=m, swap=f'{b}to{m}')
    # >1 difference: keep all, e.g. 'multi:F3L;G7A' (BP_aa pos MTP_aa)
    return dict(pos=';'.join(str(i + 1) for i, _, _ in diffs),
                bp_aa=';'.join(b for _, _, b in diffs),
                mtp_aa=';'.join(m for _, m, _ in diffs),
                swap='multi:' + ';'.join(f'{b}{i+1}{m}' for i, m, b in diffs))

def collect_records(rows, seq_col, bp_col, dataset_col):
    """Output records per unique (MTP_seq, BP_seq, Dataset);
    one blast query (qid) per unique (MTP_seq, Dataset)."""
    records = []
    seen = set()              # (slug, mtp, bp) -> output dedup
    qids = {}                 # (slug, mtp)     -> qid (one blast per unique MTP/Dataset)
    counters = {}
    for row in rows:
        mtp = re.sub(r'\s+', '', row.get(seq_col, '') or '').upper()
        bp  = re.sub(r'\s+', '', row.get(bp_col, '') or '').upper()
        dataset = (row.get(dataset_col, '') or '').strip()
        if not mtp or not dataset:
            continue
        slug = dataset_slug(dataset)
        key = (slug, mtp, bp)
        if key in seen:
            continue
        seen.add(key)
        qkey = (slug, mtp)
        if qkey not in qids:
            counters[slug] = counters.get(slug, 0) + 1
            qids[qkey] = f'{slug}_{counters[slug]}'
        records.append(dict(qid=qids[qkey], dataset=dataset, slug=slug,
                            seq=mtp, bp=bp, species=species_for(dataset),
                            **compute_swap(mtp, bp)))
    return records

def ensure_db(ref):
    if not os.path.exists(ref + '.phr'):
        subprocess.run(['makeblastdb', '-in', ref, '-dbtype', 'prot',
                        '-parse_seqids'], check=True)

def run_blast(query, ref, out_tab, threads):
    fmt = ('6 qseqid sseqid pident length mismatch gapopen '
           'qstart qend qlen stitle bitscore')
    subprocess.run(['blastp', '-task', 'blastp-short', '-query', query,
                    '-db', ref, '-outfmt', fmt, '-out', out_tab,
                    '-evalue', '1000', '-comp_based_stats', '0',
                    '-max_target_seqs', '50', '-num_threads', str(threads)],
                   check=True)

def acc_gene(sseqid, stitle):
    acc, entry = sseqid, ''
    parts = sseqid.split('|')
    if len(parts) >= 3:                 # sp|P04637|P53_HUMAN
        acc, entry = parts[1], parts[2]
    m = GN_RE.search(stitle)
    gene = m.group(1) if m else (entry.split('_')[0] if entry else acc)
    return acc, gene

def classify(hits):
    full = [h for h in hits if h['gapopen'] == 0 and h['qstart'] == 1
            and h['qend'] == h['qlen'] and h['length'] == h['qlen']]
    if any(h['mismatch'] == 0 for h in full):
        return {'status': 'drop_reference'}
    cand = sorted((h for h in full if h['mismatch'] == 1),
                  key=lambda h: -h['bitscore'])
    if not cand:
        return {'status': 'drop_no_parent'}
    best = cand[0]
    a, g = acc_gene(best['sseqid'], best['stitle'])
    if TRYP_GENE.match(g) or TRYP_TITLE.search(best['stitle']):
        return {'status': 'drop_trypsin', 'accession': a, 'gene': g}
    if IG_GENE.match(g) or IG_TITLE.search(best['stitle']):
        return {'status': 'drop_ig', 'accession': a, 'gene': g}
    genes = {acc_gene(h['sseqid'], h['stitle'])[1] for h in cand}
    return {'status': 'keep', 'accession': a, 'gene': g,
            'ambiguous': len(genes) > 1}

def process_species(species, records, ref, out_dir, threads, prefix):
    if not records:
        print(f'{species}: no records, skipping'); return
    query = os.path.join(out_dir, f'{prefix}_mtp_query_{species}.fasta')
    tab   = os.path.join(out_dir, f'{prefix}_mtp_blast_{species}.tsv')
    written = set()                         # one FASTA entry per qid
    with open(query, 'w') as f:
        for r in records:
            if r['qid'] not in written:
                written.add(r['qid'])
                f.write(f">{r['qid']}\n{r['seq']}\n")
    ensure_db(ref)
    run_blast(query, ref, tab, threads)

    by_q = {}
    for line in open(tab):
        p = line.rstrip('\n').split('\t')
        by_q.setdefault(p[0], []).append(dict(
            sseqid=p[1], length=int(p[3]), mismatch=int(p[4]),
            gapopen=int(p[5]), qstart=int(p[6]), qend=int(p[7]),
            qlen=int(p[8]), stitle=p[9], bitscore=float(p[10])))

    unfilt = os.path.join(out_dir, f'{prefix}_{species}_unfiltered.csv')
    filt   = os.path.join(out_dir, f'{prefix}_{species}_filtered.csv')
    counts = {}
    cols = ['query_id', 'dataset', 'mtp_seq', 'bp_seq', 'swap', 'position',
            'bp_aa', 'mtp_aa', 'species', 'status', 'accession', 'gene',
            'ambiguous']
    with open(unfilt, 'w', newline='') as uf, open(filt, 'w', newline='') as ff:
        wu, wf = csv.writer(uf), csv.writer(ff)
        wu.writerow(cols); wf.writerow(cols)
        for r in records:
            res = classify(by_q.get(r['qid'], []))
            counts[res['status']] = counts.get(res['status'], 0) + 1
            row = [r['qid'], r['dataset'], r['seq'], r['bp'], r['swap'],
                   r['pos'], r['bp_aa'], r['mtp_aa'], species, res['status'],
                   res.get('accession', ''), res.get('gene', ''),
                   res.get('ambiguous', '')]
            wu.writerow(row)
            if res['status'] == 'keep':
                wf.writerow(row)
    uniq = len({r['seq'] for r in records})
    print(f'{species}: {len(records)} rows ({uniq} unique MTP seqs) ->',
          dict(sorted(counts.items())))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, help='SAAP_quant_df.csv')
    ap.add_argument('--human-ref', required=True)
    ap.add_argument('--mouse-ref', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--prefix', default='April26')
    ap.add_argument('--seq-col', default='MTP_seq')
    ap.add_argument('--bp-col', default='BP_seq')
    ap.add_argument('--dataset-col', default='Dataset')
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with open(a.csv, newline='') as fh:
        rows = list(csv.DictReader(fh))
    if rows:
        missing = [c for c in (a.seq_col, a.bp_col, a.dataset_col)
                   if c not in rows[0]]
        if missing:
            raise SystemExit(f'Columns not found: {missing}. '
                             f'Available: {list(rows[0])}')

    records = collect_records(rows, a.seq_col, a.bp_col, a.dataset_col)
    groups = {'human': [], 'mouse': []}
    for r in records:
        groups[r['species']].append(r)

    process_species('human', groups['human'], a.human_ref, a.out_dir,
                    a.threads, a.prefix)
    process_species('mouse', groups['mouse'], a.mouse_ref, a.out_dir,
                    a.threads, a.prefix)

if __name__ == '__main__':
    main()
