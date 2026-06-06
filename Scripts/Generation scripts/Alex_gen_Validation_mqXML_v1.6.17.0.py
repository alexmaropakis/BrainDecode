#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re

"""
Generate MaxQuant XML file and SLURM script
Last updated by: Alex Maropakis, 05-22-2026

Example usage:
python Alex_gen_Validation_mqXML.py /home/maropakis.a/scripts/XML/templates/tmt10_Val_MS3_template.xml \
    /scratch/maropakis.a/MQ_raw/Ping_2018/ACG/b3/ \
    -o Ping2018_ACG_B3_Val.xml \
    -e Ping_2018 \
    -s S3 \
    --partition short \
    -t 21
"""

print("Defining arguments...")
parser = argparse.ArgumentParser(description='Generate MaxQuant XML and SLURM script')
parser.add_argument('input_xml', type=argparse.FileType('r', encoding='UTF-8'))
parser.add_argument('raw_file_folders', type=str, nargs='+')
parser.add_argument('-o', '--outfile', type=str, required=True)
parser.add_argument('-e', '--experiment', type=str, required=True, help='Experiment name, e.g. Takasugi_2024')
parser.add_argument('-s', '--sample', type=str, required=True, help='Sample / MQ experiment label, e.g. S3 or 3')
parser.add_argument('-mq', '--mq-version', type=str, default='2_8_0_0', help='MaxQuant version, e.g. 2_8_0_0')
parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads to use')
parser.add_argument('--species', choices=['human', 'mouse'], default='human')
parser.add_argument('--search-type', type=str, default='Val', choices=['DP', 'Val'])
parser.add_argument('--partition', type=str, default='slavov')
parser.add_argument('--mem', type=str, default='16Gb')

args = parser.parse_args()

mqpar_text = args.input_xml.read()
args.input_xml.close()


def natural_key(path):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)', os.path.basename(path))]


def raw_fraction(path, fallback):
    match = re.search(r'(?:^|[_-])FR(\d+)(?:[_\.-]|$)', os.path.basename(path), flags=re.IGNORECASE)
    return int(match.group(1)) if match else fallback


def mq_experiment_label(sample):
    return sample[1:] if sample.upper().startswith('S') and sample[1:].isdigit() else sample


# Replace FASTA path
appended_dir = '/scratch/maropakis.a/Dependencies/FASTA_appended/'

fasta_candidates = [x for x in os.listdir(appended_dir) if x.endswith('_MTP.fasta')]

outfile_key = os.path.basename(args.outfile).replace('.xml', '')

matched_fastas = []

# Ping datasets
if 'Ping2018_ACG_B' in outfile_key:
    batch = re.search(r'B(\d+)', outfile_key).group(1)
    target = f'S{batch}_ACGB{batch}_MTP.fasta'
    matched_fastas = [x for x in fasta_candidates if x == target]

elif 'Ping2018_FC_B' in outfile_key:
    batch = re.search(r'B(\d+)', outfile_key).group(1)
    target = f'S{batch}_FCB{batch}_MTP.fasta'
    matched_fastas = [x for x in fasta_candidates if x == target]

# Takasugi datasets
elif 'Takasugi_2024' in outfile_key:
    tissue = outfile_key.replace('Takasugi_2024_', '').replace('_Val', '')
    matched_fastas = [
        x for x in fasta_candidates
        if args.sample in x and tissue.lower() in x.lower()
    ]

if len(matched_fastas) == 0:
    raise FileNotFoundError(f'No appended FASTA found for outfile={args.outfile}')

if len(matched_fastas) > 1:
    raise ValueError(f'Multiple appended FASTAs matched: {matched_fastas}')

fasta_path = os.path.join(appended_dir, matched_fastas[0])

identifier_rule = '>([^\s]*)'

fasta_block = f"""<fastaFiles>
   <FastaFileInfo>
      <fastaFilePath>{fasta_path}</fastaFilePath>
      <identifierParseRule>{identifier_rule}</identifierParseRule>
      <descriptionParseRule>>(.*)</descriptionParseRule>
      <taxonomyParseRule></taxonomyParseRule>
      <variationParseRule></variationParseRule>
      <modificationParseRule></modificationParseRule>
      <taxonomyId></taxonomyId>
   </FastaFileInfo>
</fastaFiles>"""

print(f"Replacing FASTA path... {fasta_path}")
start_tag = "<fastaFiles>"
end_tag = "</fastaFiles>"

start_idx = mqpar_text.find(start_tag)
end_idx = mqpar_text.find(end_tag) + len(end_tag)

if start_idx == -1 or end_idx == -1:
    raise ValueError("Could not locate <fastaFiles> block in XML template")

mqpar_text = mqpar_text[:start_idx] + fasta_block + mqpar_text[end_idx:]

# Turn DotNet = True if not already 
print("Turning DotNet = True...")
mqpar_text = re.sub(
    r'<useDotNetCore>False</useDotNetCore>',
    '<useDotNetCore>True</useDotNetCore>',
    mqpar_text,
    flags=re.DOTALL
)

# Collect raw files in order of fraction number 
print("Collecting raw files...")
raw_files = []

for folder in args.raw_file_folders:
    for f in os.listdir(folder):
        path = os.path.join(folder, f)

        if os.path.isfile(path) and f.lower().endswith('.raw'):
            raw_files.append(path)

raw_files.sort(key=natural_key)

file_counter = len(raw_files)

if file_counter == 0:
    raise FileNotFoundError("No .raw files found in the provided raw_file_folders.")

file_paths_xml = '<filePaths>\n'

for path in raw_files:
    file_paths_xml += f'<string>{path}</string>\n'

file_paths_xml += '</filePaths>'

mqpar_text = re.sub(
    r'<filePaths>.*?</filePaths>',
    file_paths_xml,
    mqpar_text,
    flags=re.DOTALL
)

# Adjust experiments, fractions, ptms, paramGroupIndices 
print("Adjusting experiments, fractions, ptms, paramGroupIndices...")

label = mq_experiment_label(args.sample)

experiments_text = '<experiments>\n' + ''.join(
    [f'<string>{label}</string>\n' for _ in range(file_counter)]
) + '</experiments>'

fractions_text = '<fractions>\n' + ''.join(
    [f'<short>{raw_fraction(path, i + 1)}</short>\n' for i, path in enumerate(raw_files)]
) + '</fractions>'

ptms_text = '<ptms>\n' + ''.join(['<boolean>False</boolean>\n' for _ in range(file_counter)]) + '</ptms>'

group_inds_text = '<paramGroupIndices>\n' + ''.join(['<int>0</int>\n' for _ in range(file_counter)]) + '</paramGroupIndices>'

reference_channels_text = '<referenceChannel>\n' + ''.join(['<string></string>\n' for _ in range(file_counter)]) + '</referenceChannel>'

mqpar_text = re.sub(
    r'<experiments>.*?</referenceChannel>',
    experiments_text + fractions_text + ptms_text + group_inds_text + reference_channels_text,
    mqpar_text,
    flags=re.DOTALL
)

# Set up output folder 
search_name = os.path.splitext(os.path.basename(args.outfile))[0]

output_folder = f'/scratch/maropakis.a/MQ_outputs/{args.experiment}/{args.search_type}/{search_name}'

print(f"Setting up output folder... {output_folder}")

os.makedirs(output_folder, exist_ok=True)

mqpar_text = re.sub(
    r'<fixedCombinedFolder>.*?</fixedCombinedFolder>',
    f'<fixedCombinedFolder>{output_folder}</fixedCombinedFolder>',
    mqpar_text,
    flags=re.DOTALL
)

mqpar_text = re.sub(
    r'<numThreads>.*?</numThreads>',
    f'<numThreads>{args.threads}</numThreads>',
    mqpar_text,
    flags=re.DOTALL
)

# Write XML file 
xml_dir = os.path.join(os.path.expanduser('~'), 'scripts', 'XML', args.search_type)

os.makedirs(xml_dir, exist_ok=True)

xml_path = os.path.join(xml_dir, os.path.basename(args.outfile))

with open(xml_path, 'w') as f:
    f.write(mqpar_text)

# Generate SLURM script 
print(f"Generating SLURM script... {search_name}")

mq_version_dots = args.mq_version.replace('_', '.')

slurm_script = f"""#!/bin/bash
#SBATCH --job-name={search_name}
#SBATCH --output={search_name}.out
#SBATCH --error={search_name}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem}
#SBATCH --partition={args.partition}

export DOTNET_ROOT=$HOME/dotnet8
MAXQUANT_EXE=/home/maropakis.a/MQ/MaxQuant_{mq_version_dots}/bin/MaxQuantCmd.exe

srun dotnet $MAXQUANT_EXE {xml_path}
"""

slurm_path = os.path.join(
    os.path.expanduser('~'),
    'scripts',
    'Batch',
    'searches',
    search_name + '.sh')

os.makedirs(os.path.dirname(slurm_path), exist_ok=True)

with open(slurm_path, 'w') as f:
    f.write(slurm_script)

# Summary 
print(f"XML and SLURM script successfully created! {search_name}")
print(f"  FASTA:   {fasta_path}")
print(f"  Sample: {args.sample}")
print(f"  MQ experiment label: {label}")
print(f"  XML:    {xml_path}")
print(f"  SLURM:  {slurm_path}")
print(f"  Output: {output_folder}")
