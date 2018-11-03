./get_file_list.py \
--server newblue4.acrc.bris.ac.uk \
--user gh13047 \
--dirs /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/cleaned_for_elastic \
/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/ukbb_broad/cleaned_for_elastic \
/projects/MRC-IEU/research/data/ukbiobank/summary/gwas/dev/release_candidate/data/cleaned_for_elastic \
--outdir ../studies \
--neo4j-bolt bolt://ieu-db-interface.epi.bris.ac.uk:27687 \
--neo4j-user gib \
--neo4j-password aW4tNBhWKxhd


ID = [ name for name in os.listdir('studies') if os.path.isdir(os.path.join('studies', name)) ]

import json
import os
import subprocess

for id in ID:
	print(id)
	j = json.load(open(os.path.join("studies", id, "metadata.json"), "rt"))
	cmd = 'ln -s ' + os.path.join(j['elastic_file_path'], j['elastic_file']) + ' ' + os.path.join('studies', id, 'harmonised.gz')
	subprocess.call(cmd, shell=True)




for i in studies/*/metadata.json
do
echo $i
sed -i 's@: 9@: ""@g' $i
done


#!/bin/bash

./prepare_ref.sh

./clump.py \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ~/snakemake/GUGC_MetaAnalysis_Results_UA.csv.tab.gz \
--snp-col 1 \
--pval-col 7 \
--delimiter $'\t' \
--header 1 \
--gzipped 1 \
--clean \
--out gug

MarkerName,n_total,A1,A2,beta,se,p_gc

./clump.py \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas-info temp.json \
--outdir temp \
--clean

./clump.py \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ../sandpit/jointGwasMc_HDL.txt.gz \
--gzipped \
--snp-col 3 \
--pval-col 9 \
--delimiter $'\t' \
--header


./master_list.py \
--dirs temp \
--output temp/instrument-master.txt \
--bfile ../ref/data_maf0.01_rs_snps


./extract_master.r \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas-info temp.json \
--snplist temp/instrument-master.txt \
--outdir temp \
--clean

./extract_master.r \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ../sandpit/jointGwasMc_HDL.txt.gz \
--gzipped \
--snplist ../sandpit/instrument-master.txt \
--snp-col 3 \
--ncontrol-col 8 \
--oa-col 5 \
--ea-col 4 \
--pval-col 9 \
--beta-col 6 \
--se-col 7 \
--delimiter $'\t' \
--header



snakemake -pr -j 3 --cluster-config bc4-cluster.json --cluster "sbatch --partition {cluster.partition} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --mem {cluster.mem} --output {cluster.output}"


1. Create json files for all datasets
	- get filename
	- match to elastic
	- if it is not present call it an orphan
	- if it is present
		- get mr-base id
		- path
		- filename
		- gzipped or not
		- header or not
		- delimiter
		- columns

2. 



for f in *; do
if [ ! -r $f ]
then
echo $f
fi
done



