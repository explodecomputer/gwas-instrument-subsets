#!/usr/bin/env python3

import json
import paramiko
import argparse
import subprocess
import pysftp
import os
from pathlib import Path

parser = argparse.ArgumentParser(description = 'Download file')
parser.add_argument('--rdsf-config', required=False, default='')
parser.add_argument('--dir', required=False, default='')

args = parser.parse_args()

def download_gwas(rdsf_config, remotepath, localpath):
	config = json.load(open(rdsf_config, 'rt'))
	client = paramiko.SSHClient()
	client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	client.connect(config['server'], username=config['user'], password=config['password'])
	sftp = client.open_sftp()
	sftp.get(remotepath, localpath)
	sftp.close()
	client.close()

def scp(rdsf_config, remotepath, localpath):
	config = json.load(open(rdsf_config, 'rt'))
	cmd = 'sshpass -p "' + config['password'] + '" scp ' + config['user'] + '@' + config['server'] + ':' + remotepath + ' ' + localpath
	a=subprocess.call(cmd, shell=True)
	return(a)

def pysftp_method(rdsf_config, remotepath, localpath):
	config = json.load(open(rdsf_config, 'rt'))
	localdir=os.path.dirname(localpath)
	with pysftp.Connection(config['server'], username=config['user'], password=config['password']) as sftp:
		with pysftp.cd(localdir):
			sftp.get(remotepath)

DIR = vars(args)['dir']

ID = [ os.path.join(DIR, name) for name in os.listdir(DIR) if os.path.isdir(os.path.join(DIR, name)) ]

# ID = ['../studies/2', '../studies/7']

prob=open(os.path.join(DIR, '000_scp_problem.txt'), 'wt')

for id in ID:
	print(id)
	gwas_info=json.load(open(os.path.join(id, 'metadata.json'), 'rt'))
	fn1 = os.path.join(id, gwas_info['elastic_file'])
	fn2 = os.path.join(id, 'elastic.gz')
	if Path(fn1).is_file():
		print("File exists. Renaming")
		os.rename(fn1, fn2)
	elif Path(fn2).is_file():
		print("Finished")
	else:
		flag=scp(vars(args)['rdsf_config'], os.path.join(gwas_info['elastic_file_path'], gwas_info['elastic_file']), os.path.join(id, gwas_info['elastic_file']))
		if(flag != 0):
			prob.write(str(id)+"\n")

# download_gwas(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])
# scp(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])
# pysftp_method(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])


