#!/bin/bash

awk '{ if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") && ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G"))
	print $2 }' ../ref/data_maf0.01_rs.bim > snps

plink --bfile ../ref/data_maf0.01_rs --extract snps --make-bed --out ../ref/data_maf0.01_rs_snps

rm snps
