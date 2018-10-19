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

./clump.py \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ~/snakemake/2014-02-13_EAGLE_Eczema_GWAMA_CEU_fixed_without_23andMe.rsid.tab.gz \
--snp-col 1 \
--pval-col 7 \
--delimiter $'\t' \
--header 0 \
--gzipped 1 \
--out eczema \
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
--dirs ../sandpit \
--output ../sandpit/instrument-master.txt \
--bfile ../ref/data_maf0.01_rs_snps


./extract_master.r \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ../sandpit/GUGC_MetaAnalysis_Results_UA.csv \
--snplist ../sandpit/instrument-master.txt \
--snp-col 1 \
--ncontrol-col 2 \
--oa-col 4 \
--ea-col 3 \
--pval-col 7 \
--beta-col 5 \
--se-col 6 \
--delimiter , \
--header

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



snakemake -j 200 --cluster-config bc4-cluster.json --cluster "sbatch --partition {cluster.partition} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --mem {cluster.mem} --output {cluster.output}"




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



./get_file_list.py \
--server newblue4.acrc.bris.ac.uk \
--user gh13047 \
--dirs /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/cleaned_for_elastic \
/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/ukbb_broad/cleaned_597_files \
/projects/MRC-IEU/research/data/ukbiobank/summary/gwas/dev/release_candidate/data/ukb-pipeline/cleaned \
--outdir ../data \
--neo4j-bolt bolt://ieu-db-interface.epi.bris.ac.uk:27687 \
--neo4j-user gib \
--neo4j-password aW4tNBhWKxhd


