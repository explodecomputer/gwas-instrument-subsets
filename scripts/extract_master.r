#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(jsonlite))

write_out <- function(x, basename, header=FALSE)
{
	g <- gzfile(paste0(basename), "w")
	write.table(x, g, row.names=FALSE, col.names=FALSE, na="", sep=",")
	close(g)
	if(header) write.table(x[0,], file=paste0(basename, "_header.csv"), row.names=FALSE, col.names=TRUE, sep=",")
}

is_palindrome <- function(a1, a2)
{
	(a1 == "A" & a2 == "T") |
	(a1 == "T" & a2 == "A") |
	(a1 == "C" & a2 == "G") |
	(a1 == "G" & a2 == "C")
}

# create parser object
parser <- ArgumentParser()
parser$add_argument('--bfile', required=TRUE)
parser$add_argument('--gwas', required=FALSE)
parser$add_argument('--gwas-id', required=TRUE)
parser$add_argument('--instrument-list', required=TRUE)
parser$add_argument('--out', required=TRUE)
parser$add_argument('--snplist', required=TRUE)
parser$add_argument('--tag-r2', type="double", default=0.6)
parser$add_argument('--tag-kb', type="double", default=5000)
parser$add_argument('--tag-nsnp', type="double", default=5000)
parser$add_argument('--palindrome-freq', type="double", default=0.4)
parser$add_argument('--no-clean', action="store_true", default=FALSE)
parser$add_argument('--rdsf-config', required=FALSE, default='')


args <- parser$parse_args()
# args <- parser$parse_args(c("--bfile", "../ref/data_maf0.01_rs_snps", "--gwas-info", "../studies/2/metadata.json", "--snplist", "../studies/master_list.txt", "--no-clean"))

print(args)

rootname <- dirname(args[['out']])
tempname <- file.path(rootname, "temp")
outname <- args[['out']]
# rootname <- file.path(args[['outdir']], gwas_info[['id']])
# rootname <- gsub(".csv.gz$", "", args[["out"]])

# read gwas

input_file <- args[['gwas']]

input <- paste0("gunzip -c ", input_file)
gwas <- fread(input, sep="\t", header=FALSE)

# Rename gwas columns
names(gwas) <- c("snp_col", "ea_col", "oa_col", "eaf_col", "beta_col", "se_col", "pval_col", "samplesize_col")

snplist <- scan(args[["snplist"]], what="character")
gwas_1 <- gwas[gwas[["snp_col"]] %in% snplist, ]

# harmonise

ref <- fread(paste0(args[["bfile"]], ".master.frq"))
ref$beta <- 1
ref$se <- 0.1
ref$pval <- 0.1
a <- format_data(
	ref,
	type="exposure",
	snp_col="SNP",
	effect_allele_col="A1",
	other_allele_col="A2",
	eaf_col="MAF"
)

b <- format_data(gwas_1, type="outcome", 
	snp_col="snp_col",
	beta_col="beta_col",
	se_col="se_col",
	effect_allele_col="ea_col",
	other_allele_col="oa_col",
	eaf_col="eaf_col",
	samplesize_col="samplesize_col",
	pval_col="pval_col"
)

action <- is_forward_strand(b$SNP, b$effect_allele.outcome, b$other_allele.outcome, a$SNP, a$effect_allele.exposure, a$other_allele.exposure)

ab <- harmonise_data(a, b, action=action)

check_cols <- c("samplesize.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "pval.outcome", "beta.outcome", "se.outcome", "SNP")

for(i in check_cols)
{
	if(!i %in% names(ab))
	{
		message("filling in: ", i)
		ab[[i]] <- NA
	}
}

gwas_1 <- ab %>% filter(mr_keep) %$%
	data_frame(
		SNP=SNP,
		ea=effect_allele.outcome,
		oa=other_allele.outcome,
		beta=beta.outcome,
		se=se.outcome,
		eaf=eaf.outcome,
		samplesize=samplesize.outcome,
		pval=pval.outcome
	)


## missing snps

missing_snps <- snplist[!snplist %in% gwas_1[["SNP"]]]
if(length(missing_snps) > 0)
{
	# write gwas snplist + missing_list
	# write missing_list
	# run plink
	write.table(c(missing_snps, gwas[["snp_col"]]), file=paste0(tempname, ".searchspace"), row=FALSE, col=FALSE, qu=FALSE)
	write.table(missing_snps, file=paste0(tempname, ".targets"), row=FALSE, col=FALSE, qu=FALSE)
	cmd <- paste0(
		"plink --bfile ", args[["bfile"]], 
		" --extract ", tempname, ".searchspace",
		" --r2 in-phase with-freqs gz",
		" --ld-snp-list ", tempname, ".targets",
		" --ld-window-kb ", args[["tag_kb"]],
		" --ld-window-r2 ", args[["tag_r2"]],
		" --ld-window ", args[["tag_nsnp"]],
		" --out ", tempname, ".targets"
	)
	system(cmd)

	ld <- fread(paste0("gunzip -c ", tempname, ".targets.ld.gz"), header=TRUE) %>%
		filter(SNP_A != SNP_B) %>%
		mutate(PHASE=gsub("/", "", PHASE))
	temp <- do.call(rbind, strsplit(ld$PHASE, "")) %>% as_data_frame
	names(temp) <- c("x1", "y1", "x2", "y2")
	ld <- cbind(ld, temp)
	ld$palindrome <- is_palindrome(ld$y1, ld$y2)
	ld <- arrange(ld, palindrome, desc(R2)) %>%
		filter(!duplicated(SNP_A))
	gwas_2 <- gwas[gwas[["snp_col"]] %in% ld$SNP_B, ]

	a <- ld %$% data_frame(
		SNP=SNP_B,
		effect_allele=y1,
		other_allele=y2,
		eaf=MAF_B,
		beta=1,
		se=0.1,
		pval=0.1
	) %>% format_data(type="exposure")
	b <- format_data(gwas_2, type="outcome", 
		snp_col="snp_col",
		beta_col="beta_col",
		se_col="se_col",
		effect_allele_col="ea_col",
		other_allele_col="oa_col",
		eaf_col="eaf_col",
		samplesize_col="samplesize_col",
		pval_col="pval_col"
	)
	ab <- harmonise_data(a,b)

	gwas_2 <- ab %>% filter(mr_keep) %$%
		data_frame(
			proxy=SNP,
			proxy_ea=effect_allele.outcome,
			proxy_oa=other_allele.outcome,
			proxy_eaf=eaf.outcome,
			beta=beta.outcome,
			se=se.outcome,
			samplesize=samplesize.outcome,
			pval=pval.outcome
		)
	temp <- ld %$% 
	data_frame(
		SNP=SNP_A,
		proxy=SNP_B,
		proxy_r2=R2,
		eaf=MAF_A,
		ea=x1,
		oa=x2
	)
	gwas_2 <- inner_join(temp, gwas_2)
	gwas_1 <- bind_rows(gwas_1, gwas_2)
} else {
	gwas_1$proxy <- NA
	gwas_1$proxy_ea <- NA
	gwas_1$proxy_oa <- NA
	gwas_1$proxy_eaf <- NA
}

gwas_1$trait <- args[["gwas_id"]]

inst <- scan(args[["instrument_list"]], what=character())
gwas_1$instrument <- as.numeric(gwas_1$SNP %in% inst)

gwas_1 <- select(gwas_1,
	SNP, trait, instrument, ea, oa, beta, se, eaf, samplesize, pval, proxy, proxy_r2, proxy_ea, proxy_oa, proxy_eaf
)
names(gwas_1)[1] <- "snp"

message("Found SNPs: ", sum(is.na(gwas_1$proxy)))
message("Proxy SNPs: ", sum(!is.na(gwas_1$proxy)))
message("Lost SNPs: ", sum(!snplist %in% gwas_1$snp))

write_out(gwas_1, outname)

if(args[["no_clean"]] == FALSE)
{
	message("removing ancilliary files for ", tempname)
	unlink(paste0(tempname, c(".searchspace", ".targets", ".targets.ld.gz", ".targets.log", ".targets.nosex")))
}

