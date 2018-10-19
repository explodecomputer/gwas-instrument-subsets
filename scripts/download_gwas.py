#!/usr/bin/env python3

import json
import paramiko
import argparse

parser = argparse.ArgumentParser(description = 'Download file')
parser.add_argument('--rdsf-config', required=False, default='')
parser.add_argument('--remotepath', required=False, default='')
parser.add_argument('--localpath', required=False, default='')

args = parser.parse_args()

def download_gwas(rdsf_config, remotepath, localpath):
	config = json.load(open(rdsf_config, 'rt'))
	client = paramiko.SSHClient()
	client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	client.connect(config['server'], username=config['user'])
	sftp = client.open_sftp()
	sftp.get(remotepath, localpath)
	sftp.close()
	client.close()

download_gwas(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])


