#!/usr/bin/env python3

import argparse
import paramiko
import os.path
import re
import json
from py2neo import Graph

parser = argparse.ArgumentParser(description = 'List directory')
parser.add_argument('--server', required=True)
parser.add_argument('--user', required=True)
# parser.add_argument('--password', required=True)
parser.add_argument('--dirs', nargs='+', required=True)
parser.add_argument('--regex', default='.+\.gz')
parser.add_argument('--outdir', required=True)
parser.add_argument('--neo4j-bolt', required=True)
parser.add_argument('--neo4j-user', required=True)
parser.add_argument('--neo4j-password', required=True)



args = parser.parse_args()

if not os.path.exists(vars(args)['outdir']):
	os.makedirs(vars(args)['outdir'])


client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect(vars(args)['server'], username=vars(args)['user'])
sftp = client.open_sftp()

def find_files_sftp(DIR, regex, gwas_dict):
	files = list(filter(regex.search, sftp.listdir(DIR)))
	print("Found " + str(len(files)) + " files")
	for gwas_name in files:
		gwas_dict[gwas_name] = DIR
	return gwas_dict

gwas_dict = {}
for DIR in vars(args)['dirs']:
	print("Searching directory: " + DIR)
	gwas_dict = find_files_sftp(DIR, re.compile(vars(args)['regex']), gwas_dict)

print("Total files: " + str(len(gwas_dict)))

# These are all the files
with open("000_rdsf_full_list.txt", 'wt') as out:
	for f in gwas_dict:
		out.write(os.path.join(gwas_dict[f], f) + "\n")


graph = Graph(vars(args)['neo4j_bolt'], user = vars(args)['neo4j_user'], password = vars(args)['neo4j_password'])

# ids = graph.run("match (n:Study) return n.id as id, n.filename as filename, n.path as path" ).data()

ids = graph.run("match (a:Group)-[r:ACCESS_TO]->(n:Study) where NOT a.name contains 'GSK' and NOT a.name = 'GTEx' return a.name as group, n.id as id, n.filename as filename").data()

import itertools
for key, group in itertools.groupby(ids, key=lambda x:x['group']):
	print(key, len(list(group)))


# Find the mr-base id for each filename
fine = []
need_gz = []
missing = []


for d in ids:
	if d['filename'] in gwas_dict:
		fine.append(d['filename'])
	elif d['filename'] + '.gz' in gwas_dict:
		need_gz.append(d['filename'])
	else:
		missing.append(d['filename'])


len(missing)
len(need_gz)
len(fine)

with open(os.path.join(vars(args)['outdir'],"000_missing_file.txt"), "wt") as o:
	for i in missing:
		o.write(i + "\n")

present_id = []
missing_id = []
mrbase_files = [x['filename'] for x in ids]
for d in gwas_dict:
	if d in mrbase_files:
		present_id.append(os.path.join(gwas_dict[d], d))
	else:
		missing_id.append(os.path.join(gwas_dict[d], d))

with open(os.path.join(vars(args)['outdir'],"000_missing_id.txt"), 'wt') as o:
	for i in missing_id:
		o.write(i + '\n')


mrbase = []
for d in ids:
	if d['filename'] in gwas_dict:
		gzipflag = d['filename'].endswith('.gz')
		headerflag = 'cleaned_597' in gwas_dict[d['filename']]
		mrbase.append({
			'id': d['id'],
			'elastic_file': d['filename'],
			'elastic_file_path': gwas_dict[d['filename']],
			'gzipped': gzipflag,
			'header': headerflag,
			'delimiter': '\t',
			'snp_col': 1,
			'ea_col': 2,
			'oa_col': 3,
			'eaf_col': 4,
			'beta_col': 5,
			'se_col': 6,
			'pval_col': 7,
			'ncontrol_col': 8,
			'ncase_col': 9
		})
		outdir = os.path.join(vars(args)['outdir'], d['id'])
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		with open(os.path.join(outdir, 'metadata.json'), 'wt') as o:
			json.dump(
				mrbase[-1],
				o, 
				indent=2, 
				separators=(',', ': ')
			)


