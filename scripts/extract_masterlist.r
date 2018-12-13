if(!require(gwasvcftools))
{
	if(!required(devtools)) install.packages("devtools")
	devtools::install_github("MRCIEU/gwasvcftools")
}
library(gwasvcftools)
library(argparse)

# create parser object
parser <- ArgumentParser()
parser$add_argument('--snplist', required=TRUE)
parser$add_argument('--bcf-dir', required=TRUE)
parser$add_argument('--gwas-id', required=TRUE)
parser$add_argument('--out', required=TRUE)
parser$add_argument('--bfile', required=TRUE)
parser$add_argument('--get-proxies', default='yes')
parser$add_argument('--vcf-ref', required=FALSE)
parser$add_argument('--tag-r2', type="double", default=0.6)
parser$add_argument('--tag-kb', type="double", default=5000)
parser$add_argument('--tag-nsnp', type="double", default=5000)
parser$add_argument('--palindrome-freq', type="double", default=0.4)
parser$add_argument('--no-clean', action="store_true", default=FALSE)
parser$add_argument('--rdsf-config', required=FALSE, default='')
parser$add_argument('--instrument-list', default=FALSE, action="store_true")


args <- parser$parse_args()
# args <- parser$parse_args(c("--bfile", "../../vcf-reference-datasets/ukb/ukb_ref", "--gwas-id", "2", "--snplist", "temp.snplist", "--no-clean", "--out", "out", "--bcf-dir", "../../gwas-files", "--vcf-ref", "../../vcf-reference-datasets/1000g/1kg_v3_nomult.bcf", "--get-proxies"))
print(args)
tempname <- tempfile(pattern="extract", tmpdir=dirname(args[['out']]))
bcf <- file.path(args[['bcf_dir']], args[['gwas_id']], "harmonised.bcf")
snplist <- fread(args[['snplist']], header=FALSE, sep="\t")

o1 <- extract(bcf, snplist, tempname, args[['get_proxies']], args[["bfile"]], args[["vcf_ref"]])
o1$mrbaseid <- args[['gwas_id']]

if(!is.null(args[['instrument_list']]))
{
	clumpfile <- file.path(args[['bcf_dir']], args[['gwas_id']], "derived/instruments/clump.txt")
	clump <- fread(clumpfile, header=FALSE, sep="\t")
	clump$id <- paste(clump$V1, clump$V2)
	o1$instrument <- FALSE
	id <- paste(o1$CHROM, o1$POS)
	o1$instrument[id %in% clump$id] <- TRUE
}
write_out(o1, basename=args[["out"]], header=FALSE)

