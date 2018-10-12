#!/bin/bash

./prepare_ref.sh

./clump.py \
--bfile ../ref/data_maf0.01_rs_snps \
--gwas ../sandpit/GUGC_MetaAnalysis_Results_UA.csv \
--snp-col 1 \
--pval-col 7 \
--delimiter , \
--header


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
