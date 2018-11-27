import os.path
import re

# Define some variables
REF = 'ref/data_maf0.01_rs_snps'
# configfile: 'config.json'

# Find all the initial study files
# ID = ['2', '6', '7']
# ID = [1237]
ID = [ name for name in os.listdir('studies') if os.path.isdir(os.path.join('studies', name)) ]
ID1 = list(filter(lambda x: re.search('^UKB-a', x), ID))
ID2 = list(filter(lambda x: re.search('^[0-9]', x), ID))
ID = ID1 + ID2


# Setup SFTP
#from snakemake.remote.SFTP import RemoteProvider
#SFTP = RemoteProvider(username=config['user'], password=config['password'])


# Create a rule defining all the final files

rule all:
	input: 
		expand('studies/{id}/derived/instruments/master_list.csv.gz', id=ID)

# Step 1: clump each GWAS

rule clump:
	input:
		'studies/{id}/elastic.gz'
	output:
		'studies/{id}/derived/instruments/clump.txt'
	shell:
		"mkdir -p studies/{wildcards.id}/derived/instruments/; ./scripts/clump.py --bfile {REF} --gwas {input} --out {output}"


# Step 2: Create a master list of all unique instrumenting SNPs

rule master_list:
	output:
		'studies/instruments.txt'
	input:
		expand('studies/{id}/derived/instruments/clump.txt', id=ID)
	shell:
		'./scripts/master_list.py --bfile {REF} --dirs studies --output {output}'

# Step 3: Extract the master list from each GWAS

rule extract_master:
	output:
		'studies/{id}/derived/instruments/master_list.csv.gz'
	input:
		a = rules.master_list.output, b = 'studies/{id}/elastic.gz'
	shell:
		"./scripts/extract_master.r --bfile {REF} --snplist {input.a} --gwas {input.b} --out {output}"

