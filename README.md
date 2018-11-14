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
