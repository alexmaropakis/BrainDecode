# gen_fragpipe_fasta.py
# Convert MTP FASTA files from AAS pipeline Validation1 output into FragPipe-compatible FASTA with reversed decoys appended
#
# Intended usage:
#   python3 gen_fragpipe_fasta.py \
#   --fasta-dir /scratch/maropakis.a/Dependencies/FASTA_appended/ \
#   --output-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe/
#
# Input:
#   Directory containing *_MTP.fasta files
#
# Output:
#   *_fragpipe.fasta files with cleaned headers + decoys appended

import argparse
import os
import sys

# Clean MTP headers into FragPipe-safe format
def clean_header(h):
    if not h.startswith(">"):
        return h
    h = h[1:]
    if h.startswith("MTP|"):
        core = h.split("|", 1)[1]
        core = core.split("_base")[0]
        core = core.replace("|", "_")
        return ">MTP_" + core
    return ">" + h.split()[0].replace(" ", "_")

# Parse FASTA into (header, sequence) tuples
def parse_fasta(f):
    entries = []
    h = None
    s = []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if h is not None:
                entries.append((h, "".join(s)))
            h = line
            s = []
        else:
            s.append(line)
    if h is not None:
        entries.append((h, "".join(s)))
    return entries

# Write targets + reversed decoys
def write_fragpipe(entries, out_path, decoy_prefix="rev_"):
    with open(out_path, "w") as out:
        cleaned = []
        for h, seq in entries:
            ch = clean_header(h)
            cleaned.append((ch, seq))
            out.write(ch + "\n" + seq + "\n")
        for h, seq in cleaned:
            out.write(">" + decoy_prefix + h[1:] + "\n" + seq[::-1] + "\n")

# Process single file
def process_file(in_path, out_path, force=False):
    if os.path.exists(out_path) and not force:
        return False
    with open(in_path) as f:
        entries = parse_fasta(f)
    write_fragpipe(entries, out_path)
    return True

# Main
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", required=True, help="Input directory with *_MTP.fasta files")
    parser.add_argument("--output-dir", required=True, help="Output directory for FragPipe FASTAs")
    parser.add_argument("--force", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()

    if not os.path.isdir(args.fasta_dir):
        sys.exit("ERROR: fasta-dir not found")

    os.makedirs(args.output_dir, exist_ok=True)

    files = [f for f in os.listdir(args.fasta_dir) if f.endswith("_MTP.fasta")]
    if not files:
        sys.exit("ERROR: no *_MTP.fasta files found")

    written = 0
    for f in sorted(files):
        src = os.path.join(args.fasta_dir, f)
        dst = os.path.join(args.output_dir, f.replace("_MTP.fasta", "_fragpipe.fasta"))
        if process_file(src, dst, args.force):
            written += 1

    print(f"Done. Written: {written} FASTA files")

if __name__ == "__main__":
    main()