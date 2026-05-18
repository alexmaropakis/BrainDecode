# BrainDecode

BrainDecode is a workflow for detecting, validating, and analyzing substituted amino acid peptides in aging and neurodegeneration proteomics datasets.

The repository contains the code, templates, sample maps, reusable dependency tables, and analysis notebooks. Large raw/search inputs, MaxQuant output, and plot exports are kept outside the Git repository in `Project_BrainDecode`.

## Repository Layout

```text
BrainDecode/
├── Dependencies/
│   ├── Analysis_Outputs/
│   │   ├── Bai_2020/
│   │   ├── Ping_2018/
│   │   └── Takasugi_2024/
│   ├── Bai_2020/
│   ├── PD_2026/
│   ├── Ping_2018/
│   │   ├── acg/
│   │   └── fc/
│   ├── Sample_maps/
│   └── Takasugi_2024/
├── MQ_templates/
├── Scripts/
│   ├── Analysis scripts/
│   ├── Generation scripts/
│   └── Pipeline scripts/
├── LICENSE
└── README.md
```

## External Project Data

Large data files are expected in a sibling folder:

```text
Project_BrainDecode/
├── Analysis_Inputs/
│   ├── Bai_2020/
│   ├── Ping_2018/
│   └── Takasugi_2024/
├── mq_output/
└── Plots/
```

The notebooks currently use these absolute project roots:

```python
CODE_DIR = "/Users/alexmaropakis/Projects/BrainDecode/"
PROJECT_DIR = "/Users/alexmaropakis/Projects/Project_BrainDecode/"
```

`Analysis_Inputs`, `mq_output`, and `Plots` should stay in `Project_BrainDecode` rather than this Git repository.

## Key Folders

`Scripts/Generation scripts/` contains scripts that generate MaxQuant XML files, pipeline scripts, and translation resources.

`Scripts/Pipeline scripts/` contains the pipeline template scripts used for SAAP detection, validation, and quantification.

`Scripts/Analysis scripts/` contains the Jupyter notebooks for dependency generation and downstream analyses.

`MQ_templates/` contains MaxQuant XML templates.

`Dependencies/Sample_maps/` contains sample map spreadsheets used by the notebooks.

`Dependencies/Ping_2018/`, `Dependencies/Takasugi_2024/`, and `Dependencies/Bai_2020/` contain reusable analysis dependency files such as fragment dictionaries, PTM heatmap data, dataset metrics, and validation summaries.

`Dependencies/Analysis_Outputs/` contains non-plot analysis outputs such as `.xlsx`, `.tsv`, and `.p` files. Plot files should not be stored there.

## Plot Policy

All `.pdf` and `.png` outputs should go to:

```text
Project_BrainDecode/Plots/
```

`Dependencies/Analysis_Outputs/` should contain no `.pdf` or `.png` files.

## MaxQuant Output Policy

MaxQuant output should use lowercase `mq_output`:

```text
Project_BrainDecode/mq_output/
```

The notebooks expect dataset subfolders under `mq_output`, including `Ping_2018`, `Takasugi_2024`, and `Bai_2020`.

## Typical Workflow

1. Place or generate raw analysis inputs in `Project_BrainDecode/Analysis_Inputs/`.
2. Place MaxQuant output in `Project_BrainDecode/mq_output/`.
3. Use scripts in `Scripts/Generation scripts/` to generate XML or pipeline files as needed.
4. Run notebooks in `Scripts/Analysis scripts/`.
5. Save tables and reusable outputs to `Dependencies/Analysis_Outputs/`.
6. Save plots to `Project_BrainDecode/Plots/`.

## Git Notes

Avoid committing large raw data, MaxQuant output, or generated plots. Keep those in `Project_BrainDecode` or another external storage location.

Recommended exclusions:

```gitignore
.DS_Store
__pycache__/
*.pyc
Analysis_Inputs/
mq_output/
Plots/
Project_BrainDecode/
*.raw
```
