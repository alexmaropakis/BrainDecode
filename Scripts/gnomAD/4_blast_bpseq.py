#!/usr/bin/env python3
"""
4_blast_bpseq.py

Map base peptides (BP) to Ensembl-110 protein coordinates with BLASTp.

Following Tsour et al. 2024, BP peptides are searched against the Ensembl 110
proteome and only hits with 100% identity, no mismatches, no gaps, spanning the
full peptide are kept. For such a hit the protein position of each substitution
is exact (subject_start + query_offset), and the matched transcript (ENST) is
carried through for the codon lookup in stage 5.

INPUTS
    --csv      saap_input.csv from stage 3 (cols: source_dataset, bp_seq, mtp_seq, ...)
    --ref      Ensembl 110 proteome FASTA (Homo_sapiens.GRCh38.pep.all.fa)
    --out      output directory
    --threads  BLAST threads

OUTPUT
    saap_mapped.csv with columns:
        dataset, gene, transcript, protein, protein_pos, bp_aa, mtp_aa,
        bp_seq, mtp_seq

Ensembl pep.all header format (parsed for gene_symbol: and transcript: tokens):
  >ENSP00000493376.2 pep ... gene:ENSG... transcript:ENST00000632684.1 ...
   gene_symbol:TRBV20OR9-2 ...
A header lacking transcript: is treated as a fatal format mismatch
"""

import argparse
import csv
import os
import subprocess
import tempfile


def parse_ensembl_header(stitle):
    # Ensembl pep.all stitle -> (gene_symbol, transcript_id); transcript required
    gene, tx = "", ""
    for t in stitle.split():
        if t.startswith("gene_symbol:"):
            gene = t.split(":", 1)[1]
        elif t.startswith("transcript:"):
            tx = t.split(":", 1)[1]
    return gene, tx


def write_bp_fasta(rows, path):
    # rows, path -> write BP sequences as query FASTA (q{i} identifiers)
    with open(path, "w") as f:
        for i, r in enumerate(rows):
            f.write(f">q{i}\n{r['bp_seq']}\n")


def ensure_db(ref):
    # build BLAST protein db if the .phr index is missing
    if not os.path.exists(ref + ".phr"):
        subprocess.run(
            ["makeblastdb", "-in", ref, "-dbtype", "prot", "-parse_seqids"],
            check=True,
        )


def run_blast(query, ref, out, threads):
    # blastp-short with subject coordinates and full subject title
    fmt = ("6 qseqid sseqid pident length mismatch gapopen "
           "qstart qend sstart send qlen stitle bitscore")
    subprocess.run(
        ["blastp", "-task", "blastp-short", "-query", query, "-db", ref,
         "-outfmt", fmt, "-out", out, "-evalue", "1000",
         "-comp_based_stats", "0", "-max_target_seqs", "5",
         "-num_threads", str(threads)],
        check=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    rows = list(csv.DictReader(open(a.csv)))

    tmp = tempfile.mkdtemp()
    fasta = os.path.join(tmp, "bp.fasta")
    blast = os.path.join(tmp, "bp.tsv")

    write_bp_fasta(rows, fasta)
    ensure_db(a.ref)
    run_blast(fasta, a.ref, blast, a.threads)

    # collect hits per query
    hits = {}
    for line in open(blast):
        (q, s, pid, length, mm, go, qs, qe, ss, se, ql, st, bs) = \
            line.rstrip("\n").split("\t")
        hits.setdefault(q, []).append({
            "sseqid": s, "pident": float(pid), "length": int(length),
            "mismatch": int(mm), "gapopen": int(go),
            "qstart": int(qs), "qend": int(qe), "sstart": int(ss),
            "send": int(se), "qlen": int(ql), "stitle": st,
            "bitscore": float(bs),
        })

    out = os.path.join(a.out, "saap_mapped.csv")
    n_exact = n_no_hit = n_bad_header = 0

    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dataset", "gene", "transcript", "protein", "protein_pos",
                    "bp_aa", "mtp_aa", "bp_seq", "mtp_seq"])

        for i, r in enumerate(rows):
            qid = f"q{i}"
            bp, mtp = r["bp_seq"], r["mtp_seq"]
            if qid not in hits:
                n_no_hit += 1
                continue

            # 100% identity, 0 mismatch, 0 gap, full-length,
            # query anchored at residue 1. Among qualifying hits take best score.
            cand = [h for h in hits[qid]
                    if h["pident"] == 100.0 and h["mismatch"] == 0
                    and h["gapopen"] == 0 and h["length"] == h["qlen"]
                    and h["qstart"] == 1 and h["qend"] == h["qlen"]]
            if not cand:
                n_no_hit += 1
                continue
            h = max(cand, key=lambda x: x["bitscore"])

            gene, tx = parse_ensembl_header(h["stitle"])
            if not tx:
                n_bad_header += 1     # wrong proteome format (e.g. UniProt)
                continue
            prot = h["sseqid"].split("|")[-1] if "|" in h["sseqid"] else h["sseqid"]
            sstart = h["sstart"]

            # one row per differing residue; protein_pos exact
            for j, (b, m) in enumerate(zip(bp, mtp)):
                if b == m:
                    continue
                protein_pos = sstart + j
                w.writerow([r.get("source_dataset", ""), gene, tx, prot,
                            protein_pos, b, m, bp, mtp])
                n_exact += 1

    print(out)
    print(f"{n_exact} substitution rows from exact BP hits")
    print(f"{n_no_hit} queries with no qualifying (100%/0mm/0gap/full) hit")
    if n_bad_header:
        print(f"WARNING: {n_bad_header} hits lacked transcript: in header "
              f"(wrong proteome format?)")


if __name__ == "__main__":
    main()
