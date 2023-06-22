# cfDNA_midpoint_coverage 

This repository contains a snakemake workflow to calculate different scores through the nucleosome positioning of cell-free DNA signals. It is based on the repository https://github.com/kircherlab/cfDNA and was adjusted and updated as part of a master thesis. This repository originally calculated the Windowed Protection Score published by Snyder *et al.* (Cell, 2016). As part of the thesis additionally the option was added to calculate the Midpoint Coverage Score by Doebley *et al.* (Nature Communications, 2022).

## Table of Contents

- [cfDNA\_midpoint\_coverage](#cfdna_midpoint_coverage)
  - [Table of Contents](#table-of-contents)
  - [Usage](#usage)
    - [Step 1: Obtain a copy of this workflow](#step-1-obtain-a-copy-of-this-workflow)
    - [Step 2: Configure workflow](#step-2-configure-workflow)
    - [Step 3: Install Snakemake](#step-3-install-snakemake)
    - [Step 4: Execute workflow](#step-4-execute-workflow)
  - [Workflows](#workflows)
    - [Transcription Factor Binding Site Filtering](#transcription-factor-binding-site-filtering)
    - [Nucleosome Positioning Scores](#nucleosome-positioning-scores)
      - [Input](#input)
      - [Output](#output)



## Usage

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) the repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, `samples.tsv` to specify your sample setup and `regions.tsv` to specify target regions.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

The workflows are executed from the repository root folder. The different analyses have to be executed separately. To specify the respective workflow use the `-s` switch followed by the path of the `Snakefile` (e.g.: `./snakefile_WPS.smk`)

Activate the conda environment:

```bash
conda activate snakemake
```

Test your configuration by performing a dry-run via

```bash
snakemake -s path/to/Snakefile --use-conda -n
```

Execute the workflow locally via

```bash
snakemake -s path/to/Snakefile --use-conda --cores $N
```

using `$N` cores or run it in a cluster environment via

```bash
snakemake -s path/to/Snakefile --use-conda --cluster qsub --jobs 100
```

or

```bash
snakemake -s path/to/Snakefile --use-conda --drmaa --jobs 100
```

If you not only want to fix the software stack but also the underlying OS, use

```bash
snakemake --use-conda --use-singularity
```

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for further details.

## Workflows

### Transcription Factor Binding Site Filtering

The workflow to calculate nucleosome positioning scores (explained below) needs as input samples (e.g. cfDNA sequencing data in .bam format) and regions of interest. These regions can for example be TFBSs. This repository contains a workflow to preprocess and filter specific TFBSs downloaded from the Gene Transcription Regulatory Database (GTRD). These preprocessed TFBSs can then be used as regions of interest in the nucleosome positioning score workflow.

**Note:** More information about the download and preprocessing of the GTRD files [here](resources/TFBS_download/download.md/)

### Nucleosome Positioning Scores

The workflow analyzes the nucleosome positioning based on fragmentation patterns of cfDNA at regions of interest. It calculates the Windowed Protection Score based on the publication by Snyder *et al.* (Cell, 2016) and the Midpoint Coverage Score based on the publication by Doebley *et al.* (Nature Communications, 2022). It takes as input one or more sample files in .bam format and on or more region file in .bed format. During the analysis, the fragmentation pattern of every sample is determined at every region. Every region is extendd to an additional 1,000 bp up- and downstream. These enlarged windows are intersected with an exclusion list and all sites not overlapping that exclusion list are kept.

#### Input

- configured by the user ([samples.tsv](config/samples.tsv)):
    - analysis ID
    - samples
    - path to sample .bam files
    - reference samples for plotting
    - genome build per sample
    - weights (-g indicates that the sample .bam files were corrected for GC bias with https://github.com/kircherlab/cfDNA-UniFlow)
- configured by the user ([regions.tsv](config/regions.tsv)):
    - bed file containing regions of interest (e.g. TFBS), all having the same length

**Note:** More information about config files [here](config/README.md)

#### Output

- table containing bp specific Windowed Protection Score and Midpoint Coverage Score for regions listed in bed
- line plot showing normalized Windowed Protection Score and Midpoint Coverage Score of multiple samples