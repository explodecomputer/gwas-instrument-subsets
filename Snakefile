import os.path

REF = 'ref/data_maf0.01_rs_snps'
RDSF_CONFIG = 'rdsf_config.json'

ID = [ name for name in os.listdir('studies') if os.path.isdir(os.path.join('studies', name)) ]

ID = ['2', '6', '7']

# Create a rule defining all the final files

rule all:
	input: 
		expand('studies/{id}/master_list.csv.gz', id=ID)


# Step 1: clump each GWAS

rule clump:
	input:
		'studies/{id}/metadata.json'
	output:
		'studies/{id}/clump.txt'
	shell:
		"./scripts/clump.py --bfile {REF} --gwas-info {input} --rdsf-config {RDSF_CONFIG}"


# Step 2: Create a master list of all unique instrumenting SNPs

rule master_list:
	output:
		'studies/instruments.txt'
	input:
		expand('studies/{id}/clump.txt', id=ID)
	shell:
		'./scripts/master_list.py --bfile {REF} --dirs studies --output {output}'

# Step 3: Extract the master list from each GWAS

rule extract_master:
	output:
		'studies/{id}/master_list.csv.gz'
	input:
		a = rules.master_list.output, b = 'studies/{id}/metadata.json'
	shell:
		"./scripts/extract_master.r --bfile {REF} --gwas-info {input.b} --snplist {input.a} --rdsf-config {RDSF_CONFIG}"

