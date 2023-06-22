# Download and preprocessing of transcription factor binding sites

## Download the TFBS from GTRD

Bash commands to download the meta cluster data set in either version 19.10 or 20.06. 

**GTRD version 19.10:**
wget gtrd.biouml.org/downloads/19.10/chip-seq/Homo%20sapiens_meta_clusters.interval.gz -O resources/TFBS_download/Homo_sapiens_meta_clusters_v19.interval.gz

**GTRD version 20.06:**
wget gtrd.biouml.org/downloads/20.06/intervals/chip-seq/Homo_sapiens_meta_clusters.zip -O resources/TFBS_download/Homo_sapiens_meta_clusters_v20.interval.gz


## Preprocessing of the downloaded files

For the preprocessing and filtering of the transcription factor binding sites, specify the version in the config/config_filtering.yml file. The workflow is executed with snakemake. Be sure that snakemake is installed and the correct conda environment is activated.

Execute a dry-run to test configurations with:
```
snakemake -s snakefile_TFBS_filtering.smk -n
```
If the dry run works, execute the workflow locally by specifiying the number of cores the workflow can use.
```
snakemake -s snakefile_TFBS_filtering.smk -j <cores>
```

The workflow only works for GTRD version 19.10 or 20.06.