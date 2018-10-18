#!/usr/bin/env python3

import argparse
import subprocess
import gzip
import os
import logging
import json

parser = argparse.ArgumentParser(description = 'Extract and clump top hits')
parser.add_argument('--bfile', required=True)
parser.add_argument('--gwas', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--pval-threshold', type=float, default=5e-8)
parser.add_argument('--snp-col', type=int, required=True)
parser.add_argument('--pval-col', type=int, required=True)
parser.add_argument('--delimiter', default=' ')
parser.add_argument('--gzipped', type=int, default=0)
parser.add_argument('--header', type=int, default=0)
parser.add_argument('--clump-r2', type=float, default=0.001)
parser.add_argument('--clump-kb', type=float, default=1000)

args = parser.parse_args()
# args = parser.parse_args(['--bfile', '../ref/data_maf0.01_rs', '--gwas', '../GUGC_MetaAnalysis_Results_UA.csv', '--snp-col', '1', '--pval-col', '7', '--header', '--delimiter', ','])


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler(vars(args)['out']+'.clump-log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info(json.dumps(vars(args), indent=1))


snplist=vars(args)['out']

if vars(args)['gzipped'] == 1:
	f = gzip.open(vars(args)['gwas'], 'rt')
else:
	f = open(vars(args)['gwas'], 'rt')

if vars(args)['header'] == 1:
	f.readline()

o = open(snplist, 'wt')
o.write('SNP P\n')

n=0
for line in f:
	x = line.strip().split(vars(args)['delimiter'])
	try:
		if float(x[vars(args)['pval_col'] - 1]) < vars(args)['pval_threshold']:
			o.write(x[vars(args)['snp_col'] - 1] + ' ' + x[vars(args)['pval_col'] - 1] + '\n')
			n+=1
	except:
		logger.info("Error parsing line")

o.close()
f.close()

logger.info("found " + str(n) + " hits")

if n > 0:
	x = ('plink --bfile ' + vars(args)['bfile'] +
	' --clump ' + snplist +
	' --clump-kb ' + str(vars(args)['clump_kb']) +
	' --clump-r2 ' + str(vars(args)['clump_r2']) +
	' --clump-p1 ' + str(vars(args)['pval_threshold']) +
	' --clump-p2 ' + str(vars(args)['pval_threshold']) +
	' --out ' + snplist)
	subprocess.call(x, shell=True)
	f = open(snplist + '.clumped', 'rt')
	o = open(snplist + '.clumped.snplist', 'wt')
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
		os.remove(snplist + '.clumped.snplist')
	except OSError:
		pass
	open(snplist + '.snplist', 'a').close()


