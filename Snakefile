import os.path

# Define some variables
REF = 'ref/data_maf0.01_rs_snps'
configfile: 'config.json'

# Find all the initial study files
# ID = [ name for name in os.listdir('studies') if os.path.isdir(os.path.join('studies', name)) ]
ID = ['2', '6', '7']


# Setup SFTP
from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username=config['user'], password=config['password'])


# Create a rule defining all the final files

rule all:
	input: 
		expand('studies/{id}/master_list.csv.gz', id=ID)

# Step 1: clump each GWAS

rule clump:
	input:
		a='studies/{id}/metadata.json', b=SFTP.remote('newblue4.acrc.bris.ac.uk/panfs/panasas01/sscm/gh13047/repo/gwas-instrument-subsets/studies/{id}/harmonised.gz')
	output:
		'studies/{id}/clump.txt'
	shell:
		"./scripts/clump.py --bfile {REF} --gwas-info {input.a} --gwas {input.b}"


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
		a = rules.master_list.output, b = 'studies/{id}/metadata.json', c = SFTP.remote('newblue4.acrc.bris.ac.uk/panfs/panasas01/sscm/gh13047/repo/gwas-instrument-subsets/studies/{id}/harmonised.gz')
	shell:
		"./scripts/extract_master.r --bfile {REF} --gwas-info {input.b} --snplist {input.a} --gwas {input.c}"

