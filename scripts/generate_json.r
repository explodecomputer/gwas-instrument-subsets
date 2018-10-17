# for every study in available_outcomes
## get the study id
## check if the file exists
## if it does, create a yaml file that can be fed into a wrapper for clump etc

# library(TwoSampleMR)
# ao <- available_outcomes()
# ao <- subset(ao, access != "developer")

# fil <- expand.grid(filename=ao$filename, path=possible_paths, suffix=c("", ".gz"))
# fil$fullpath <- file.path(fil$path, paste0(fil$filename, fil$suffix))
# fil$exist <- file.exists(fil$fullpath)

# table(fil$exist)
# library(dplyr)
# subset(fil, !exist) %>% head



# missing <- subset(ao, !filename %in% subset(fil, exist)$filename)$filename

suppressPackageStartupMessages(library(argparse))
library(dplyr)
library(jsonlite)

write.json <- function(input, output, outdir)
{
	bn <- basename(input)
	fn <- file.path(outdir, output)
	out <- list()
	out$name <- bn
	out$file <- file.path(path, filename)
	out$delimiter <- "$'\t'"
	out$header <- FALSE
	out$gzipped <- TRUE
	out$snp_col <- 1
	out$ea_col <- 2
	out$oa_col <- 3
	out$eaf_col <- 4
	out$beta_col <- 5
	out$se_col <- 6
	out$pval_col <- 7
	out$ncontrol_col <- 8
	out %>% toJSON(pretty=TRUE)	%>%
		writeLines(fn)
}

parser <- ArgumentParser()
parser$add_argument('--outdir', required=TRUE)

objs <- lapply(possible_paths, function(x)
{
	fn <- list.files(path=x, pattern="*.gz")
	data_frame(path=x, filename=fn)
}) %>% bind_rows


