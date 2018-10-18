

REF = 'ref/data_maf0.01_rs_snps'
SERVER = 'newblue4.acrc.bris.ac.uk'

import paramiko
import os.path
from snakemake.remote.SFTP import RemoteProvider
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect(SERVER, username=config['user'])
sftp = client.open_sftp()

def find_files_sftp(DIR, gwas_dict):
	for filename in sftp.listdir(DIR):
		if filename.endswith('.gz'):
			gwas_name = filename.replace('.gz', '')
			gwas_dict[gwas_name] = os.path.join(DIR, str(filename))
	return gwas_dict

gwas_dict = {}
gwas_dict = find_files_sftp('/panfs/panasas01/sscm/gh13047/snakemake', gwas_dict)

# gwas_dict = find_files_sftp('/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/cleaned_for_elastic', gwas_dict)
# gwas_dict = find_files_sftp('/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/ukbb_broad/cleaned_597_files', gwas_dict)
# gwas_dict = find_files_sftp('/projects/MRC-IEU/research/data/ukbiobank/summary/gwas/dev/release_candidate/data/ukb-pipeline/cleaned', gwas_dict)

# Gather cleaned gwas files if the dataset is local:

# import glob
# from pathlib import Path   

# def find_files_local(DIR, gwas_dict):
# 	gwas_path = Path(DIR).glob('*.gz')
# 	for f in gwas_path:
# 		gwas_name = f.name.replace('.gz', '')
# 		gwas_dict[gwas_name] = str(f)
# 	return gwas_dict

# gwas_dict = {}
# gwas_dict = find_files('/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/cleaned_for_elastic', gwas_dict)
# gwas_dict = find_files('/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/ukbb_broad/cleaned_597_files', gwas_dict)
# gwas_dict = find_files('/projects/MRC-IEU/research/data/ukbiobank/summary/gwas/dev/release_candidate/data/ukb-pipeline/cleaned', gwas_dict)


SFTP = RemoteProvider(username=config['user'], password=config['password'])

# Create directories



# Create a rule defining all the final files

rule all:
	input: 
		expand('data/{gwas_name}.harmonised.csv.gz', gwas_name=gwas_dict.keys())



# Step 1: clump each GWAS

rule clump:
	input:
		lambda wildcards: SFTP.remote(SERVER + gwas_dict[wildcards.gwas_name])
	output:
		'data/{gwas_name}.snplist'
	shell:
		"./scripts/clump.py --bfile {REF} --gwas {input} --out {output} --snp-col 1 --pval-col 7 --delimiter $'\\t' --gzipped 1 --header 0"


# Step 2: Create a master list of all unique instrumenting SNPs

rule master_list:
	output:
		'data/instrument-master.txt'
	input:
		expand('data/{gwas_name}.snplist', gwas_name=gwas_dict.keys())
	shell:
		'./scripts/master_list.py --bfile {REF} --dirs data --output {output}'

# Step 3: Extract the master list from each GWAS

rule extract_master:
	output:
		'data/{gwas_name}.harmonised.csv.gz'
	input:
		a = rules.master_list.output, b = lambda wildcards: SFTP.remote(SERVER + gwas_dict[wildcards.gwas_name])
	shell:
		"./scripts/extract_master.r --bfile {REF} --gwas {input.b} --out {output} --snplist {input.a} --snp-col 1 --ea-col 2 --oa-col 3 --eaf-col 4 --beta-col 5 --se-col 6 --pval-col 7 --ncontrol-col 8 --delimiter $'\\t' --gzipped 1 --header 0"



