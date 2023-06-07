### Preprocessing workflow ###

import glob, os.path
import pandas as pd

configfile: "config/config_filtering.yml"

rule all:
	input:
		expand("results/TFBS_bedFile/Homo_sapiens_meta_clusters_{version}.bed.gz", version=config["version"]),
		regions=expand("config/regions_GTRD_{number}_{version}.tsv",number=config["number_of_sites"],version=config["version"])

rule prepare_19:
	input:
		download_file=config["file"]
	output:
		bed_file="results/TFBS_bedFile/Homo_sapiens_meta_clusters_v19.bed.gz"
	shell:
		"""
		zcat {input.download_file} | awk 'BEGIN{{OFS="\t";FS="\t"}} {{print $1,$2,$3,$6,$13,"."}}' | tail -n +2 | sort -k1,1 -k2,2n | gzip -c > {output.bed_file}
		"""

rule prepare_20:
	input:
		download_file=config["file"]
	output:
		bed_file="results/TFBS_bedFile/Homo_sapiens_meta_clusters_v20.bed.gz"
	shell:
		"""
		zcat {input.download_file} | awk 'BEGIN{{OFS="\t";FS="\t"}} {{print $1,$2,$3,$14,$11,"."}}' | tail -n +2 | sort -k1,1 -k2,2n | gzip -c > {output.bed_file}
		"""

rule separateTFs:
	input:
		bed_file=expand("results/TFBS_bedFile/Homo_sapiens_meta_clusters_{version}.bed.gz",version=config["version"])
	output:
		regions=expand("config/regions_GTRD_{number}_{version}.tsv",number=config["number_of_sites"],version=config["version"])
	params:
		out_dir=expand("results/tf_{number}_{version}/",number=config["number_of_sites"],version=config["version"]),
		verbose="True",
		sites=config["valid_sites"],
		number=config["number_of_sites"],
	shell:
		"""
		python3 workflow/scripts/TFBS_filtering/separateInputForTFs_griffin.py -i {input.bed_file} -o {params.out_dir} \
		-v {params.verbose} -s {params.sites} -n {params.number} -r {output.regions}
		"""
#scripts/tf_linecount.sh {output} {params.out_dir}

# rule filterMetaCluster:
# 	input:
# 		gtrd19="resources/Homo_sapiens_meta_clusters_v19.interval.gz"
# 	output:
# 		bed19="results/input/Homo_sapiens_meta_clusters_v19.bed.gz"
# 	shell:
# 		"""
# 		zcat {input.gtrd19} | awk 'BEGIN{{OFS="\t";FS="\t"}} {{print $1,$2,$3,$6,$12,$13,"."}}' | tail -n +2 | sort -k1,1 -k2,2n | gzip -c > {output.bed19}
# 		"""

# rule subtract_blacklist:
# 	input:
# 		bed="results/input/{sample}.bed.gz",
# 		blacklist="results/exclusion_data/blacklist_exclusionFile.bed.gz"
# 	output:
# 		bed_sub="results/subtract/{sample}.subtract_blacklist.bed.gz"
# 	shell:
# 		"""
# 		bedtools subtract -A -a {input.bed} -b {input.blacklist} | sort -k1,1 -k2,2n | gzip -c  > {output.bed_sub}
# 		"""

# rule separateTFs_blacklist:
# 	input:
# 		bed_sub="results/subtract/Homo_sapiens_meta_clusters_v19.subtract_blacklist.bed.gz"
# 	output:
# 		"results/Homo_sapiens_meta_clusters_v19_blacklist.counts.txt"
# 	params:
# 		out_dir="results/tf_v19_bl",
# 		verbose="True"
# 	shell:
# 		"""
# 		python3 scripts/separateInputForTFs_griffin.py -i {input.bed_sub} -o {params.out_dir} -v {params.verbose}
# 		scripts/tf_linecount.sh {output} {params.out_dir}
# 		"""
