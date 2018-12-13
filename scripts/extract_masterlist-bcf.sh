#!/bin/bash

#SBATCH --job-name=ml
#SBATCH --array=0-1000
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-00:30:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=8G

echo "Running on ${HOSTNAME}"
module add languages/r/3.4.4

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}
i=$((i + 1000))

cd ${HOME}/mr-eve/gwas-instrument-subsets/scripts
ids=($(cat ids.txt))

echo ${#ids[@]}
id=`echo ${ids[$i]}`
echo $id
if [ -z "$id" ]
then
	echo "outside range"
	exit
fi

Rscript extract_masterlist.r \
--snplist ../../gwas-files/instruments.txt \
--bcf-dir ../../gwas-files \
--out ../../gwas-files/$id/derived/instruments/ml.csv.gz \
--bfile ../../vcf-reference-datasets/ukb/ukb_ref \
--vcf-ref ../../vcf-reference-datasets/1000g/1kg_v3_nomult.bcf \
--gwas-id $id \
--instrument-list \
--get-proxies yes



