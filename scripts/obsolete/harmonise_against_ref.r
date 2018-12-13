#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vcfR))
library(methods)
library(utils)

# create parser object
parser <- ArgumentParser()

parser$add_argument('--ref_file', required=TRUE)
parser$add_argument('--ref_build', required=TRUE)
parser$add_argument('--ref_info', required=TRUE)
parser$add_argument('--mrbase_id', required=TRUE)
parser$add_argument('--gwas_file', required=TRUE)
parser$add_argument('--gzipped', required=TRUE, type="integer", default=1)
parser$add_argument('--delimiter', default="\t", required=TRUE)
parser$add_argument('--skip', required=TRUE, type="integer", default=0)
parser$add_argument('--dbsnp_field', type="integer", required=TRUE)
parser$add_argument('--ea_field', type="integer", required=FALSE)
parser$add_argument('--nea_field', type="integer", required=TRUE)
parser$add_argument('--ea_af_field', type="integer", required=FALSE)
parser$add_argument('--effect_field', type="integer", required=FALSE)
parser$add_argument('--se_field', type="integer", required=FALSE)
parser$add_argument('--pval_field', type="integer", required=FALSE)
parser$add_argument('--n_field', type="integer", required=FALSE)
parser$add_argument('--out_type', required=TRUE, default="bcf")
parser$add_argument('--out', required=TRUE)

args <- parser$parse_args()
args[["gwas_header"]] <- as.logical(args[["gwas_header"]])
jlog <<- list(
	args=args,
	counts=list()
)

print(args)

read_dat <- function(filename, skip, delimiter, gzipped, snp, nea, ea, ea_af, effect, se, pval, n, jlog)
{
	if(gzipped)
	{
		dat <- data.table::fread(paste0("gunzip -c ", filename), header=FALSE, skip=skip, sep=delimiter)
	} else {
		dat <- data.table::fread(filename, header=FALSE, skip=skip, sep=delimiter)
	}
	nc <- ncol(dat)
	jlog[['counts']][['total_variants']] <<- nrow(dat)

	if(snp == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		snp <- ncol(dat)
	}
	if(nea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		nea <- ncol(dat)
	}
	if(ea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea <- ncol(dat)
	}
	if(ea_af == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea_af <- ncol(dat)
	}
	if(effect == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		effect <- ncol(dat)
	}
	if(se == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		se <- ncol(dat)
	}
	if(pval == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		pval <- ncol(dat)
	}
	if(n == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		n <- ncol(dat)
	}
	o <- format_data(
		dat, 
		type="outcome", 
		phenotype_col="outcome",
		snp_col=names(dat)[snp],
		beta_col=names(dat)[effect],
		se_col=names(dat)[se],
		effect_allele_col=names(dat)[ea],
		other_allele_col=names(dat)[nea],
		eaf_col=names(dat)[ea_af],
		pval_col=names(dat)[pval],
		samplesize_col=names(dat)[n],
	)
	o$beta.outcome[!is.finite(o$beta.outcome)] <- NA
	o$se.outcome[!is.finite(o$se.outcome)] <- NA
	o$pval.outcome[!is.finite(o$pval.outcome)] <- NA
	o$eaf.outcome[!is.finite(o$pval.outcome)] <- NA
	jlog[['counts']][['variants_not_read']] <<- nrow(dat) - nrow(o)
	ind <- is.finite(o$beta.outcome) & is.finite(o$pval.outcome)
	o <- o[ind,]
	jlog[['counts']][['variants_with_missing_stats']] <<- sum(!ind)
	jlog[['counts']][['variants_with_missing_pvals']] <<- sum(is.na(o$pval.outcome))
	o <- subset(o, !is.na(pval.outcome))

	return(o)
}



get_ref <- function(reference_file, snplist, outfile, remove_dup_rsids=TRUE)
{
	out1 <- paste0(outfile, ".snplist")
	out2 <- paste0(outfile, ".ref")
	write.table(snplist, file=out1, row=F, col=F, qu=F)

	# Extract relevent info
	cmd <- paste0("time bcftools view -i'ID=@", out1, "' ", reference_file, " | bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' > ", out2)
	print(cmd)
	system(cmd)
	unlink(out1)
	a <- data.table::fread(out2, header=FALSE, sep="\t")
	unlink(out2)
	names(a) <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")

	if(remove_dup_rsids)
	{
		a <- subset(a, !duplicated(ID))
	}
	a$beta <- 1
	a$se <- 0.1
	a$pval <- 0.1
	return(a)
}

# Read in gwas data
gwas <- read_dat(
	args[["gwas_file"]],
	skip=args[["skip"]],
	snp=args[["dbsnp_field"]],
	gzipped=args[["gzipped"]],
	delimiter=args[["delimiter"]],
	ea=args[["ea_field"]],
	nea=args[["nea_field"]],
	ea_af=args[["ea_af_field"]],
	effect=args[["effect_field"]],
	se=args[["se_field"]],
	pval=args[["pval_field"]],
	n=args[["n_field"]],
	jlog
)
gwas$mr_keep <- TRUE

# Read in ref
ref <- get_ref(args[["ref_file"]], gwas$SNP, args[["out"]]) 
a <- TwoSampleMR::format_data(
	ref,
	type="exposure",
	snp_col="ID",
	effect_allele_col="ALT",
	other_allele_col="REF",
	eaf_col="AF"
)

# Check strand
action <- TwoSampleMR::is_forward_strand(gwas$SNP, gwas$effect_allele.outcome, gwas$other_allele.outcome, ref$ID, ref$REF, ref$ALT, threshold=0.9)

# Harmonise
dat <- TwoSampleMR::harmonise_data(a, gwas, action, jlog=jlog)
jlog[['counts']][['harmonised_variants']] <- sum(dat$mr_keep)
jlog[['counts']][['variants_not_harmonised']] <- sum(!dat$mr_keep)

gwas_h <- dat %>% subset(mr_keep) %$%
	dplyr::data_frame(
		ID=SNP,
		ALT=effect_allele.exposure,
		REF=other_allele.exposure,
		BETA=beta.outcome,
		SE=se.outcome,
		PVALUE=pval.outcome,
		AF=eaf.outcome,
		N=samplesize.outcome) %>%
	dplyr::inner_join(subset(ref, select=c(ID,REF,ALT,CHROM,POS)), by=c("ID", "REF", "ALT"))
jlog[['counts']][['total_remaining_variants']] <- nrow(gwas_h)


# Create vcf format
vcf <- TwoSampleMR::make_vcf(
	ID = gwas_h$ID,
	ALT = gwas_h$ALT,
	REF = gwas_h$REF,
	B = gwas_h$BETA,
	SE = gwas_h$SE,
	PVAL = gwas_h$PVALUE,
	N = gwas_h$N,
	CHROM = gwas_h$CHROM,
	POS = gwas_h$POS,
	AF = gwas_h$AF,
	QUAL = rep(NA, nrow(gwas_h)),
	FILTER = rep("PASS", nrow(gwas_h)),
	build = args[["ref_build"]],
	meta_data = jlog
)

# Write vcf
TwoSampleMR::write_vcf(vcf, paste0(args[["out"]], ".", args[["out_type"]]))
