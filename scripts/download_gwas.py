#!/usr/bin/env python3

import json
import paramiko
import argparse
import subprocess
import pysftp
import os

parser = argparse.ArgumentParser(description = 'Download file')
parser.add_argument('--rdsf-config', required=False, default='')
parser.add_argument('--remotepath', required=False, default='')
parser.add_argument('--localpath', required=False, default='')

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
	subprocess.call(cmd, shell=True)


def pysftp_method(rdsf_config, remotepath, localpath):
	config = json.load(open(rdsf_config, 'rt'))
	localdir=os.path.dirname(localpath)
	with pysftp.Connection(config['server'], username=config['user'], password=config['password']) as sftp:
		with pysftp.cd(localdir):
			sftp.get(remotepath)

# download_gwas(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])
# scp(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])
pysftp_method(vars(args)['rdsf_config'], vars(args)['remotepath'], vars(args)['localpath'])


