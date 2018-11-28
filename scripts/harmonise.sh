#!/bin/bash

set -e

id=$1
ref=$2
# ref="../reference/1000g/1kg_v3_nomult.bcf"


gunzip -c ../studies/$id/elastic.gz | cut -f 1 > ../studies/$id/snplist.txt
wc -l ../studies/$id/snplist.txt

time bcftools view -i"ID=@../studies/$id/snplist.txt" $ref | bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' | sed '1 i\
CHROM\tPOS\tID\tREF\tALT\tAF
' | gzip -c > ../studies/$id/ref_extract.txt.gz

rm -f ../studies/$id/snplist.txt

Rscript harmonise_against_ref.r \
--ref-file ../studies/$id/ref_extract.txt.gz \
--ref-build b37 \
--gwas-file ../studies/$id/elastic.gz \
--gwas-header F \
--gwas-snp 1 \
--gwas-ref 3 \
--gwas-alt 2 \
--gwas-af 4 \
--gwas-beta 5 \
--gwas-se 6 \
--gwas-pval 7 \
--gwas-n0 8 \
--gwas-n1 0 \
--out ../studies/2/harmonised.bcf

rm -f ../studies/$id/ref_extract.txt.gz

