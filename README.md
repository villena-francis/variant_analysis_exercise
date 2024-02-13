# Analysis of genomic variants of a patient with advanced squamous cell lung cancer (Part I)

## Overview
In this case, we will process multiple FASTQ files from Whole Exome Sequencing (WES) of tumour and normal samples to generate Variant Call Format (VCF) files.

The starting files for this training were not provided due to their large size, initial state of working directory should look something like this:

```
.
├── env
│   └── variant_analysis_exercise.yml
├── genome
│   └── hg38_chr12.fa
├── intervals
│   └── S07604514_Padded_v6_chr12.bed
├── raw_data
│   ├── normal_1.fq.gz
│   ├── normal_2.fq.gz
│   ├── tumor_1.fq.gz
│   └── tumor_2.fq.gz
├── README.md
└── scripts
    ├── cleanup.sh
    └── pipeline.sh
```

## Conda environment setup
The specifications of the conda environment with the necessary programs for this exercise are in `env/variant_analysis_exercise.yml`. To install it, use the following command:

```
conda env create -f env/variant_analysis_exercise.yml
```

## Usage

To run the pipeline, navigate to the main directory and execute `pipeline.sh` with bash:

```bash
bash scripts/pipeline.sh
```
To rerun the pipeline with different samples, execute `cleanup.sh` to clear the working directory and manually organize the new files into their respective folders.

```bash
bash scripts/cleanup.sh
```