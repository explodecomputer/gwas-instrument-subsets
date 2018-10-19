#!/usr/bin/env python3

import argparse
import subprocess
import gzip
import os
import logging
import json

parser = argparse.ArgumentParser(description = 'Extract and clump top hits')
parser.add_argument('--bfile', required=True)
parser.add_argument('--gwas-info', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--pval-threshold', type=float, default=5e-8)
parser.add_argument('--clump-r2', type=float, default=0.001)
parser.add_argument('--clump-kb', type=float, default=1000)
parser.add_argument('--clean', action='store_true', default=False)

args = parser.parse_args()
# args = parser.parse_args(['--bfile', '../ref/data_maf0.01_rs', '--gwas', '../GUGC_MetaAnalysis_Results_UA.csv', '--snp-col', '1', '--pval-col', '7', '--header', '--delimiter', ','])

if not os.path.exists(vars(args)['outdir']):
	os.makedirs(vars(args)['outdir'])

gwas_info = json.load(open(vars(args)['gwas_info'], 'rt'))
rootname = os.path.join(vars(args)['outdir'], gwas_info['id'])

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler(rootname+'.clump-log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info(json.dumps(vars(args), indent=1))


snplist=rootname

if gwas_info['gzipped'] == 1:
	f = gzip.open(vars(args)['gwas'], 'rt')
else:
	f = open(os.path.join(gwas_info['elastic_file_path'], gwas_info['elastic_file']), 'rt')

if gwas_info['header'] == 1:
	f.readline()

o = open(snplist + '.tophits', 'wt')
o.write('SNP P\n')

n=0
for line in f:
	x = line.strip().split(gwas_info['delimiter'])
	try:
		if float(x[gwas_info['pval_col'] - 1]) < vars(args)['pval_threshold']:
			o.write(x[gwas_info['snp_col'] - 1] + ' ' + x[gwas_info['pval_col'] - 1] + '\n')
			n+=1
	except:
		logger.info("Error parsing line")

o.close()
f.close()

logger.info("found " + str(n) + " hits")

if n > 0:
	x = ('plink --bfile ' + vars(args)['bfile'] +
	' --clump ' + snplist + '.tophits' +
	' --clump-kb ' + str(vars(args)['clump_kb']) +
	' --clump-r2 ' + str(vars(args)['clump_r2']) +
	' --clump-p1 ' + str(vars(args)['pval_threshold']) +
	' --clump-p2 ' + str(vars(args)['pval_threshold']) +
	' --out ' + snplist)
	subprocess.call(x, shell=True)
	f = open(snplist + '.clumped', 'rt')
	o = open(snplist, 'wt')
	f.readline()
	n=0
	for line in f:
		try:
			o.write(line.strip().split()[2] + '\n')
			n+=1
		except:
			pass
	o.close()
	f.close()
	logger.info("found " + str(n) + " clumps")
else:
	try:
		os.remove(snplist)
	except OSError:
		pass
	open(snplist, 'a').close()

if vars(args)['clean'] is True:
	for f in [snplist + x for x in ['.tophits', '.log', '.nosex', '.clumped']]:
		try:
			os.remove(f)
		except OSError:
			pass

