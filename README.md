# BrainDecode: Computational workflow for detecting, validating, and analyzing amino acid substitutions from alternate RNA decoding in aging and neurodegeneration.

This repository provides tools to identify, validate, and quantify amino acid substitutions in LC-MS proteomics data that arise from alternative RNA decoding. The described pipelines evaluate multiplexed LC-MS proteomics data described in [Bai et al. (2020)](10.1016/j.neuron.2019.12.015), [Ping et al. (2018)](https://doi.org/10.1038/sdata.2018.36), and [Takasugi et al. (2024)](10.1038/s41467-024-52845-x).

This project builds on the work of the **Slavov Laboratory**:

* [Decode Website](https://decode.slavovlab.net) &nbsp; | &nbsp; [Preprint article](https://doi.org/10.1101/2024.08.26.609665)
* [Decode_Pipeline](https://github.com/SlavovLab/decode/tree/main/decode_pipeline)

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
├── frag_output/
├── Supplementary_Tables/
└── Plots/
```

## Key Folders

`Scripts/Generation scripts/` contains scripts that generate MaxQuant XML files, pipeline scripts, and translation resources.

`Scripts/Pipeline scripts/` contains the pipeline template scripts used for SAAP detection, validation, and quantification.

`Scripts/Analysis scripts/` contains the Jupyter notebooks for dependency generation and downstream analyses.

`MQ_templates/` contains MaxQuant XML templates.

`Dependencies/Sample_maps/` contains sample map spreadsheets used by the notebooks.

`Dependencies/Ping_2018/`, `Dependencies/Takasugi_2024/`, and `Dependencies/Bai_2020/` contain reusable analysis dependency files such as fragment dictionaries, PTM heatmap data, dataset metrics, and validation summaries.

`Dependencies/Analysis_Outputs/` contains non-plot analysis outputs such as `.xlsx`, `.tsv`, and `.p` files. Plot files should not be stored there.

## Analysis Workflows
1. Place or generate raw analysis inputs in `Project_BrainDecode/Analysis_Inputs/`.
2. Place MaxQuant output in `Project_BrainDecode/mq_output/`.
3. Use scripts in `Scripts/Generation scripts/` to generate XML or pipeline files as needed.
4. Run notebooks in `Scripts/Analysis scripts/`.
5. Save tables and reusable outputs to `Dependencies/Analysis_Outputs/`.
6. Save plots to `Project_BrainDecode/Plots/`.
7. Save supplementary tables to `Project_BrainDecode/Supplementary_Tables/`. 

## Data Generation Workflows
### Step 1: Custom protein databases 
Use RNA-seq data matched to LC-MS proteomics data to create sample-specific protein databases.

The code for this step is in [custom_protein_database_pipeline](https://github.com/SlavovLab/decode/tree/main/custom_protein_database_pipeline) and the [README.md](custom_protein_database_pipeline/README.md) in that directory contains detailed instructions for running the code.

If no matched RNA-seq data is available, this step can be skipped, but caution should be taken in interpreting quantified amino acid substitutions as there is lower confidence that they are not encoded in the genome.

### Step 2: Identifying modified peptides with MaxQuant or DIA-NN

[Guide for running MaxQuant in Linux.](https://github.com/alexmaropakis/MaxQuant_Linux_Tutorial)
The dependent peptide search algorithm in MaxQuant is used to identify peptides with modifications in LC-MS proteomics data. 

The LC-MS prteomics data is ideally searched against the sample-specific database generated in Step 1. If not available, species-specific UniProt FASTA can be used. 

A sample MaxQuant parameter file is provided in MaxQuant_templates, along with a script to create a new parameter file with user-defined parameters (raw files, fasta, etc.)

The output from this dependent peptide search is required to proceed with the next steps of the pipeline.

### Step 3: Identifying candidate alternate translation events
Search for modified peptides in dependent peptide search results that may represent amino acid substitutions. Add candidate peptides to custom protein sequence databases for validation search.


### Step 4. Validation search with MaxQuant, FragPipe, or another proteomics data search engine
Run a standard database search using the protein databases appended with candidate substituted peptides (step 3).

A sample MaxQuant parameter file is provided in MaxQuant_templates. The output from this validation search is required to proceed with the next steps of the pipeline.


### Step 5. Quantify alternative RNA decoding events


### Step 6. Downstream data analysis
