#!/bin/bash

set -e

id=$1
gwasdir=$2
ldref=$3


cat <(echo -e "SNP P") <(bcftools view -i 'PVAL<5e-8' $gwasdir/$id/harmonised.bcf | bcftools query -f'%ID %PVAL\n') > $gwasdir/$id/derived/instruments/tophits.txt
wc -l $gwasdir/$id/derived/instruments/tophits.txt

plink \
--bfile $ldref \
--clump $gwasdir/$id/derived/instruments/tophits.txt \
--clump-kb 10000 \
--clump-r2 0.001 \
--clump-p1 5e-8 \
--clump-p2 5e-8 \
--out $gwasdir/$id/derived/instruments/tophits.txt

awk '{ print $3 }' $gwasdir/$id/derived/instruments/tophits.txt.clumped | sed '/^[[:space:]]*$/d' > $gwasdir/$id/derived/instruments/snplist.txt

bcftools view -i "ID=@$gwasdir/$id/derived/instruments/snplist.txt" $gwasdir/$id/harmonised.bcf | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\tb37\t%ID\n' > $gwasdir/$id/derived/instruments/clump.txt
rm $gwasdir/$id/derived/instruments/snplist.txt

rm -f $gwasdir/$id/derived/instruments/tophits.txt*

