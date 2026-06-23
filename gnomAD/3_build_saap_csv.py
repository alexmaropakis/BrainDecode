#!/usr/bin/env python3
"""
3_build_saap_csv.py

Walk experiment directories, extract MTP/BP sequence pairs from the
MTP_quant_dict.p pickles produced by the AAS pipeline (validation2 +
quantification must have been run first), and consolidate them into a single
CSV used as the BLAST input for stage 4.

INPUTS
    --root    Top-level directory holding experiment folders with .p files.
    --prefix  Folder-name prefix filter (e.g. Ping_, ACG_, FC_).
    --out     Output CSV path.
    --target  Pickle filename to read (default: MTP_quant_dict.p).

OUTPUT (CSV)
    source_dataset  TMT plex identifier inferred from the folder (e.g. ACGb1)
    saap_id         key index in the quant dictionary
    bp_seq          base peptide sequence
    mtp_seq         mistranslated / variant peptide sequence
    pickle_file     origin file path (traceability)

EXAMPLE
    python3 3_build_saap_csv.py \
      --root   /scratch/maropakis.a/AAS_Pipeline/ \
      --prefix Ping_ \
      --out    /scratch/maropakis.a/Dependencies/gnomAD_pipeline/mtp_maps/saap_input.csv
"""

import argparse
import csv
import glob
import os
import pickle
import re

## Helper functions 
def load_pickle(path):
    # path -> loaded dict of SAAP entries
    with open(path, "rb") as f:
        return pickle.load(f)


def extract_plex_id(path):
    # full pickle path -> TMT plex identifier (e.g. ACGb1, FCB2), else folder name
    folder = os.path.basename(os.path.dirname(path))
    m = re.search(r"(ACG|FC)[A-Za-z]*b\d+", folder)
    return m.group(0) if m else folder


def find_pickles(root, prefix):
    # root, prefix -> list of .p paths under matching experiment folders
    paths = []
    for d in os.listdir(root):
        full = os.path.join(root, d)
        if not os.path.isdir(full):
            continue
        if not d.startswith(prefix):
            continue
        paths.extend(glob.glob(os.path.join(full, "**", "*.p"), recursive=True))
    return paths

## Main processing 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Where experimental dirs are located")
    ap.add_argument("--prefix", required=True, help="Prefix for experimental dirs")
    ap.add_argument("--out", required=True, help="where to put output csv")
    ap.add_argument("--target", default="MTP_quant_dict.p", help="Contains all SAAP-BP pairs for a dataset")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    pickle_files = find_pickles(a.root, a.prefix)

    n_written = 0
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source_dataset", "saap_id", "bp_seq", "mtp_seq", "pickle_file"])

        for pkl in pickle_files:
            # restrict to the main quant dictionary file
            # useful if multiple .p files in experimental dir
            if a.target and a.target not in pkl:
                continue
            try:
                data = load_pickle(pkl)
            except Exception as e:
                print(f"skip {pkl}: {e}")
                continue

            source_dataset = extract_plex_id(pkl)
            for k, v in data.items():
                bp = v.get("BP_seq", "")
                mtp = v.get("MTP_seq", "")
                if not bp or not mtp:        # skip incomplete entries
                    print("incomplete entires found...")
                    continue
                w.writerow([source_dataset, k, bp, mtp, pkl])
                n_written += 1

    print(f"written: {a.out} ({n_written} SAAP rows)")


if __name__ == "__main__":
    main()
