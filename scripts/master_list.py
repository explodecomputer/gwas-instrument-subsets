#!/usr/bin/env python3

import argparse
import logging
import itertools
import glob
import subprocess
import json
import os

parser = argparse.ArgumentParser(description = 'Extract and clump top hits')
parser.add_argument('--dirs', nargs='+', required=True)
parser.add_argument('--output', required=True)


###args=parser.parse_args(['--dirs', '../../gwas-files', '--output', '../../gwas-files/instrument-master.txt'])

args = parser.parse_args()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler(vars(args)['output']+'.clump-log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info(json.dumps(vars(args), indent=1))

filelist = []
for d in vars(args)['dirs']:
	filelist += glob.glob(d + '/*/derived/instruments/clump.txt')

logger.info("found " + str(len(filelist)) + " snp lists")

snps = {}
for f in filelist:
	n = set(line.strip() for line in open(f, 'rt'))
	logging.info(f)
	logging.info("number of loaded snps: " + str(len(n)))
	count = len(snps)
	snps = set(itertools.chain(snps, n))
	new_count = len(snps)
	logging.info("number of new snps: " + str(new_count - count))

logging.info("total instrument count: " + str(new_count))

o = open(vars(args)['output'], 'wt')
[o.write(x + '\n') for x in snps]
o.close()


