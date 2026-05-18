#!/usr/bin/env python3

"""
Generate experiment-specific SAAP pipeline scripts from templates.

Scans every .py file in the templates directory, applies path
substitutions, and writes ready-to-run scripts plus a SLURM
submission script.

Template filtering by --mode:
  '_TMT_' in filename       -> included only in 'tmt' mode
  '_LabelFree_' in filename -> included only in 'labelfree' mode
  all other .py files       -> included in both modes

Species filtering by --species:
  human -> uses normal Detection and Validation1 templates
  mouse -> uses Mouse Detection and Mouse Validation1 templates
  all later steps, such as Validation2 and Quant, use the normal templates

MQ_dir path structure:
  MQ_outputs/<mq_folder>/DP/<search_name>_DP/combined/txt/
  MQ_outputs/<mq_folder>/Val/<search_name>_Val/combined/txt/
  The _DP or _Val suffix is auto-appended based on which placeholder
  the template contains.

FASTA naming (--fasta):
  Provide the base name (e.g. S1_ACGB1). The script appends:
    _noMTP.fasta  for the reference FASTA
    _MTP.fasta    for the appended output

TMT example:
  python gen_pipeline.py \
      --experiment Ping_2018_ACGb1 --mode tmt --species human \
      --mq_folder Ping_2018 --search_name Ping2018_ACG_B1 \
      --sample_map sample_map_acgb1 \
      --fasta S1_ACGB1

Mouse TMT example:
  python gen_pipeline.py \
      --experiment Takasugi_2024_Heart --mode tmt --species mouse \
      --mq_folder Takasugi_2024 --search_name Takasugi_2024_Heart \
      --sample_map sample_map_heart \
      --fasta S3_Takasugi2024

Last updated by: Alex Maropakis, 05-18-2026
"""

import argparse
import os
import re
import sys
from glob import glob


## Set directories

HOME = '/home/maropakis.a'
SCRATCH = '/scratch/maropakis.a'

DEFAULT_TEMPLATE_DIR = os.path.join(HOME, 'scripts', 'Pipeline', 'templates')
DEFAULT_OUTPUT_DIR = os.path.join(HOME, 'scripts', 'Pipeline', 'experiments')
DEFAULT_SLURM_DIR = os.path.join(HOME, 'scripts', 'Batch', 'python')


## Find templates

STEP_ORDER = ['Detection', 'Validation1', 'Validation2', 'Quant']


def _step_sort_key(filename):
    """Return (priority, filename) so scripts run in pipeline order."""
    name_upper = filename.upper()
    for i, keyword in enumerate(STEP_ORDER):
        if keyword.upper() in name_upper:
            return (i, filename)
    return (len(STEP_ORDER), filename)


def _is_mouse_template(filename):
    """Return True for species-specific mouse template names."""
    return '_Mouse_' in filename or '_mouse_' in filename


def _is_species_specific_step(filename):
    """Only Detection and Validation1 have species-specific templates."""
    name_upper = filename.upper()
    return 'DETECTION' in name_upper or 'VALIDATION1' in name_upper


def _species_template_ok(filename, species):
    """
    For Detection and Validation1:
      mouse -> keep Mouse templates only
      human -> keep non-Mouse templates only

    For all other steps:
      keep normal templates only
    """
    is_mouse = _is_mouse_template(filename)

    if _is_species_specific_step(filename):
        if species == 'mouse':
            return is_mouse
        return not is_mouse

    return not is_mouse


def discover_templates(template_dir, mode, species):
    """
    Find all .py files in template_dir.

    Filter by:
      mode    -> tmt / labelfree
      species -> human / mouse for Detection and Validation1 only

    Returns list sorted in pipeline execution order.
    """
    all_py = sorted(glob(os.path.join(template_dir, '*.py')))
    selected = []

    for path in all_py:
        name = os.path.basename(path)
        is_tmt = '_TMT_' in name or '_tmt_' in name
        is_lf = '_LabelFree_' in name or '_labelfree_' in name or '_Labelfree_' in name

        if is_tmt and mode != 'tmt':
            continue
        if is_lf and mode != 'labelfree':
            continue
        if not _species_template_ok(name, species):
            continue

        selected.append(name)

    selected.sort(key=_step_sort_key)
    return selected


def make_output_name(template_name, experiment):
    """Quant_v3.py -> Quant_v3_PD_2026.py"""
    stem, ext = os.path.splitext(template_name)
    return f'{stem}_{experiment}{ext}'


## Define substitution patterns for paths

def make_substitutions(args):
    """
    Build (compiled_regex, replacement) pairs for all placeholder paths.
    Patterns that don't match a given file are silently skipped.
    """
    subs = []

    # MQ_dir:
    # 'MQ_outputs/    /DP/combined/txt/'
    # -> 'MQ_outputs/<mq_folder>/DP/<search_name>_DP/combined/txt/'
    #
    # 'MQ_outputs/    /Val/combined/txt/'
    # -> 'MQ_outputs/<mq_folder>/Val/<search_name>_Val/combined/txt/'
    subs.append((
        re.compile(r"(MQ_outputs/)[ \t]+(/(DP|Val)/)combined/txt/"),
        lambda m: f"{m.group(1)}{args.mq_folder}{m.group(2)}{args.search_name}_{m.group(3)}/combined/txt/"
    ))

    # TMT detection templates append 'combined/txt/...' later.
    # 'MQ_outputs/    /' -> 'MQ_outputs/<mq_folder>/DP/<search_name>_DP/'
    subs.append((
        re.compile(r"(MQ_outputs/)[ \t]+(/)"),
        lambda m: f"{m.group(1)}{args.mq_folder}/DP/{args.search_name}_DP{m.group(2)}"
    ))

    # aas_dir: 'AAS_Pipeline/  '
    subs.append((
        re.compile(r"(AAS_Pipeline/)\s+'"),
        rf"\g<1>{args.experiment}/'"
    ))

    # sample_map xlsx: 'Dependencies/sample_map/           .xlsx'
    if args.sample_map:
        subs.append((
            re.compile(r"(Dependencies/sample_map/)\s+(\.xlsx)"),
            rf"\g<1>{args.sample_map}\g<2>"
        ))

    # tissues txt: 'Dependencies/sample_map/    .txt'
    if args.tissue_list:
        subs.append((
            re.compile(r"(Dependencies/sample_map/)\s+(\.txt)"),
            rf"\g<1>{args.tissue_list}\g<2>"
        ))

    # FASTA paths derive from --fasta.
    if args.fasta:
        subs.append((
            re.compile(r"(database_dir\s*\+\s*')\s+(\.fasta')"),
            rf"\g<1>{args.fasta}_noMTP\g<2>"
        ))

        subs.append((
            re.compile(r"(\{s\}_)\s+(MTP\.fasta)"),
            rf"\g<1>{args.fasta}_\g<2>"
        ))

    return subs


def apply_substitutions(text, subs):
    for pattern, replacement in subs:
        text = pattern.sub(replacement, text)
    return text


## Generate slurm script

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}_%j.out
#SBATCH --error={job_name}_%j.err
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --partition={partition}

echo "=============================="
echo "{job_name}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "=============================="

cd {script_dir}
srun python {script_name}

echo "=============================="
echo "Finished: $(date)"
echo "=============================="
"""


def generate_slurm(args, script_name, script_dir):
    """Generate a single SLURM .sh file for one pipeline script."""
    job_name = os.path.splitext(script_name)[0]
    return SLURM_TEMPLATE.format(
        job_name=job_name,
        mem=args.mem,
        cpus=args.cpus,
        partition=args.partition,
        script_dir=script_dir,
        script_name=script_name,
    )


## Generate pipeline files

def main():
    parser = argparse.ArgumentParser(
        description='Generate experiment-specific SAAP pipeline scripts from templates.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    parser.add_argument('--experiment', required=True,
        help='Experiment name (e.g. PD_2026, Bai_2020)')
    parser.add_argument('--mode', required=True, choices=['tmt', 'labelfree'],
        help="Pipeline mode: 'tmt' or 'labelfree'")
    parser.add_argument('--species', default='human', choices=['human', 'mouse'],
        help="Species for Detection and Validation1 templates only (default: human)")
    parser.add_argument('--mq_folder', required=True,
        help='Top-level MQ output folder (e.g. Ping_2018)')
    parser.add_argument('--search_name', required=True,
        help='MQ search subfolder base name (e.g. Ping2018_ACG_B1). '
             '_DP and _Val suffixes are appended automatically.')

    # Mode-specific
    parser.add_argument('--sample_map', default=None,
        help='[TMT] Sample map filename without .xlsx (under Dependencies/sample_map/)')
    parser.add_argument('--tissue_list', default=None,
        help='[Label-free] Tissue list filename without .txt (under Dependencies/sample_map/)')
    parser.add_argument('--fasta', default=None,
        help='FASTA base name without suffixes (e.g. S1_ACGB1).')

    # Paths
    parser.add_argument('--template_dir', default=DEFAULT_TEMPLATE_DIR,
        help=f'Template directory (default: {DEFAULT_TEMPLATE_DIR})')
    parser.add_argument('--output_dir', default=None,
        help=f'Output directory for .py scripts (default: {DEFAULT_OUTPUT_DIR}/<experiment>/)')
    parser.add_argument('--slurm_dir', default=DEFAULT_SLURM_DIR,
        help=f'Output directory for SLURM .sh script (default: {DEFAULT_SLURM_DIR})')

    # SLURM options
    parser.add_argument('--mem', default='100G',
        help='SLURM memory (default: 100G)')
    parser.add_argument('--cpus', default='8',
        help='SLURM CPUs per task (default: 8)')
    parser.add_argument('--partition', default='short',
        help='SLURM partition (default: short)')

    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(DEFAULT_OUTPUT_DIR, args.experiment)

    if args.mode == 'tmt' and not args.sample_map:
        parser.error('--sample_map is required for TMT mode')
    if args.mode == 'labelfree' and not args.tissue_list:
        parser.error('--tissue_list is required for label-free mode')

    # -- Discover templates --
    templates = discover_templates(args.template_dir, args.mode, args.species)
    if not templates:
        print(f'ERROR: No .py templates found in {args.template_dir} for mode "{args.mode}" and species "{args.species}"')
        sys.exit(1)

    print(f'\nFound {len(templates)} template(s) for mode "{args.mode}" and species "{args.species}":')
    for t in templates:
        print(f'  {t}')
    print()

    # -- Build substitutions --
    subs = make_substitutions(args)

    # -- Process each template --
    os.makedirs(args.output_dir, exist_ok=True)
    output_scripts = []

    for template_name in templates:
        template_path = os.path.join(args.template_dir, template_name)
        with open(template_path, 'r') as f:
            text = f.read()

        text = apply_substitutions(text, subs)

        out_name = make_output_name(template_name, args.experiment)
        out_path = os.path.join(args.output_dir, out_name)
        output_scripts.append(out_name)

        with open(out_path, 'w') as f:
            f.write(text)
        print(f'  {template_name}  ->  {out_name}')

    # -- Generate SLURM scripts --
    os.makedirs(args.slurm_dir, exist_ok=True)
    slurm_paths = []

    for script_name in output_scripts:
        slurm_name = os.path.splitext(script_name)[0] + '.sh'
        slurm_path = os.path.join(args.slurm_dir, slurm_name)
        slurm_text = generate_slurm(args, script_name, os.path.abspath(args.output_dir))

        with open(slurm_path, 'w') as f:
            f.write(slurm_text)
        slurm_paths.append(slurm_path)
        print(f'  SLURM  ->  {slurm_path}')

    ## Summary 
    print(f'''
{"=" * 60}
  Experiment:    {args.experiment}
  Mode:          {args.mode}
  Species:       {args.species}
  MQ folder:     {args.mq_folder}
  Search name:   {args.search_name}
  Scripts:       {args.output_dir}
  SLURM:         {args.slurm_dir}
  SLURM files:   {len(slurm_paths)}
{"=" * 60}
''')


if __name__ == '__main__':
    main()