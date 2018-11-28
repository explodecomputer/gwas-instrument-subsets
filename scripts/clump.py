#!/usr/bin/env python3

import argparse
import subprocess
import gzip
import os
import logging
import json
import glob
import paramiko

parser = argparse.ArgumentParser(description = 'Extract and clump top hits')
parser.add_argument('--bfile', required=True)
parser.add_argument('--gwas', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--pval-threshold', type=float, default=5e-8)
parser.add_argument('--clump-r2', type=float, default=0.001)
parser.add_argument('--clump-kb', type=float, default=1000)
parser.add_argument('--no-clean', action='store_true', default=False)

args = parser.parse_args()

gwasfile = vars(args)['gwas']
outfile = vars(args)['out']
outdir = os.path.dirname(outfile)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler(os.path.join(outdir, 'clump.log'))
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info(json.dumps(vars(args), indent=1))

# Expect gzipped file
f = gzip.open(gwasfile, 'rt')
# Read header
f.readline().strip().split("\t")

o = open(outfile + '.tophits', 'wt')
o.write('SNP P\n')

n=0
h=0
p=0
for line in f:
	x = line.strip().split("\t")
	n+=1
	try:
		if float(x[6]) < vars(args)['pval_threshold']:
			o.write(x[0] + ' ' + x[6] + '\n')
			h+=1
	except:
		p+=1

o.close()
f.close()


logger.info("found " + str(n) + " lines")
logger.info("found " + str(p) + " problem lines")
logger.info("found " + str(h) + " hits")

if h > 0:
	x = ('plink --bfile ' + vars(args)['bfile'] +
	' --clump ' + outfile + '.tophits' +
	' --clump-kb ' + str(vars(args)['clump_kb']) +
	' --clump-r2 ' + str(vars(args)['clump_r2']) +
	' --clump-p1 ' + str(vars(args)['pval_threshold']) +
	' --clump-p2 ' + str(vars(args)['pval_threshold']) +
	' --out ' + outfile)
	subprocess.call(x, shell=True)
	try:
		f = open(outfile + '.clumped', 'rt')
		o = open(outfile, 'wt')
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
	except:
		open(outfile, 'a').close()
else:
	try:
		os.remove(outfile)
	except OSError:
		pass
	open(outfile, 'a').close()

if vars(args)['no_clean'] is False:
	for f in glob.glob(outfile + ".*"):
		os.remove(f)


