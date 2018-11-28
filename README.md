## 1. Get the datasets

```

python get_file_list.py

./get_file_list.py \
--server newblue4.acrc.bris.ac.uk \
--user gh13047 \
--dirs /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/cleaned_for_elastic \
/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/ukbb_broad/cleaned_for_elastic \
/projects/MRC-IEU/research/data/ukbiobank/summary/gwas/dev/release_candidate/data/cleaned_for_elastic \
--outdir ../studies \
--neo4j-bolt bolt://ieu-db-interface.epi.bris.ac.uk:27687 \
--neo4j-user gib \
--neo4j-password aW4tNBhWKxhd


Example run

```
snakemake -r \
-j 4 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --partition {cluster.partition} \
  --nodes {cluster.nodes} \
  --cpus-per-task {cluster.cpus-per-task} \
  --time {cluster.time} \
  --mem {cluster.mem} \
  --output {cluster.output}"
```
