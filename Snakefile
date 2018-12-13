import os.path
import re

# Define some variables
# REF = '../reference/1000g_filtered/data_maf0.01_rs_snps'
REF = '../vcf-reference-datasets/1000g_filtered/data_maf0.01_rs_snps'
VCFREF = '../../vcf-reference-datasets/1000g/1kg_v3_nomult.bcf'
GWASDIR = '../gwas-files'


# configfile: 'config.json'

# Find all the initial study files
# ID = [1237]
ID = [ name for name in os.listdir('../gwas-files') if os.path.isdir(os.path.join('../gwas-files', name)) ]
ID1 = list(filter(lambda x: re.search('^UKB-a', x), ID))
ID2 = list(filter(lambda x: re.search('^[0-9]', x), ID))
ID = ID1 + ID2


ID = ['2', '6', '7']


# Setup SFTP
#from snakemake.remote.SFTP import RemoteProvider
#SFTP = RemoteProvider(username=config['user'], password=config['password'])


# Create a rule defining all the final files

rule all:
	input: 
		expand('{GWASDIR}/{id}/derived/instruments/ml.csv.gz', GWASDIR=GWASDIR,id=ID)

# Step 1: clump each GWAS

rule clump:
	input:
		'{GWASDIR}/{id}/elastic.gz'
	output:
		'{GWASDIR}/{id}/derived/instruments/clump.txt'
	shell:
		"mkdir -p studies/{wildcards.id}/derived/instruments/; ./clump.sh {wildcards.id} ../../gwas-files {REF}"

# Step 2: Create a master list of all unique instrumenting SNPs

rule master_list:
	output:
		'{GWASDIR}/instruments.txt'
	input:
		expand('{GWASDIR}/{id}/derived/instruments/clump.txt', GWASDIR=GWASDIR,id=ID)
	shell:
		'./scripts/master_list.py --dirs {GWASDIR} --output {output}'

# Step 3: Extract the master list from each GWAS

rule extract_master:
	output:
		'studies/{id}/derived/instruments/ml.csv.gzi'
	input:
		a = rules.master_list.output, b = 'studies/{id}/elastic.gz', c = rules.clump.output
	shell:
		"Rscript scripts/extract_masterlist.r --snplist {GWASDIR}/instruments.txt --bcf-dir {GWASDIR} --out {output} --bfile {REF} --vcf-ref {VCFREF} --gwas-id {wildcards.id} --instrument-list --get-proxies yes"
