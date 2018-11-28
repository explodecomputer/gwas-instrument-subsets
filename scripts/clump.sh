#!/bin/bash

set -e

id=$1
ldref=$2

cat <(echo -e "SNP P") <(bcftools view -i 'PVAL<5e-8' ../studies/$id/harmonised.bcf | bcftools query -f'%ID %PVAL\n') > ../studies/$id/derived/instruments/tophits.txt
wc -l ../studies/$id/derived/instruments/tophits.txt

plink \
--bfile $ldref \
--clump ../studies/$id/derived/instruments/tophits.txt \
--clump-kb 10000 \
--clump-r2 0.001 \
--clump-p1 5e-8 \
--clump-p2 5e-8 \
--out ../studies/$id/derived/instruments/tophits.txt 

awk '{ print $3 }' ../studies/$id/derived/instruments/tophits.txt.clumped | sed '/^[[:space:]]*$/d' > ../studies/$id/derived/instruments/clump.txt

rm -f ../studies/$id/derived/instruments/tophits.txt*

