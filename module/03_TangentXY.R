# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Part of the TangentXY tool for normalizing and converting read 
#              depths from WES/SNP array data to log2 relative copy values.
#              Code to perform Tangent on autosomes and chrX using transformed
#              reference plane and on chrY using sex-matched reference plane.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

library('optparse')

option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-s", "--sif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-n", "--nt.matrix"), type="character", help="normal samples transformed signal matrix file", metavar="character"),
  make_option(c("-p", "--nsvd"), type="character", help="SVD processed normal samples signal file", metavar="character"),
  make_option(c("-t", "--tt.matrix"), type="character", help="tumor samples transformed signal matrix file", metavar="character"),
  make_option(c("-l", "--latent.n"), type="integer", help="number of latent factors to reconstruct normal subspace", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir <- opt$directory
setwd(dir)

suppressPackageStartupMessages(library('tidyverse'))
library('tidyverse')

suppressPackageStartupMessages(library('here'))
library('here')

if (!file.exists(here('output/TangentXY'))) {
  dir.create(here('output/TangentXY'), recursive=TRUE)
}


sif.file <- opt$sif
n.df.file <- opt$nt.matrix
n.autox.svd <- opt$nsvd
t.df.file <- opt$tt.matrix
num.lf <- opt$latent.n

sif <- read.delim(sif.file)
t.df <- read_delim(t.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')
n.autox.svd <- readRDS(n.autox.svd)

## Tangent on autosomes and chrX
cat('\nRunning Tangent on autosomes and chrX ...\n')

## Reconstruct a normal subspace with a specific number of latent factors
if (num.lf==0 | num.lf > length(n.autox.svd$d)) {
  stop('-l option need to be greater than 0, and less than or equal to the number of normal samples.')
} else if (num.lf==1) {
  N.autox <- as.matrix(n.autox.svd$u[, num.lf]) %*% n.autox.svd$d[num.lf] %*% t(n.autox.svd$v[, num.lf])
} else {
  N.autox <- n.autox.svd$u[, 1:num.lf] %*% diag(n.autox.svd$d[1:num.lf]) %*% t(n.autox.svd$v[,1:num.lf])  
}

T.autox <- t.df[!grepl('^Y', rownames(t.df)),] %>%
  as.matrix()

## Get the origin in the normal subspace
N.autox.means <- apply(N.autox, 1, mean)
N.autox0 <- N.autox - N.autox.means
T.autox0 <- T.autox - N.autox.means

## Run Tangent
Npi.autox <- pracma::pinv(N.autox0)
weights.autox <- Npi.autox %*% T.autox0
proj.autox <- N.autox0 %*% weights.autox
T.autox.norm <- T.autox0 - proj.autox

## Re-scaling after Tangent by median
T.autox.norm.medians <- T.autox.norm[!grepl('X', rownames(T.autox.norm)),] %>%
  apply(., 2, median)

T.autox.norm.rescaled <- t(t(T.autox.norm)- T.autox.norm.medians)

## Adjust male chrX so that it is relative to CN=2
male.tumor.samples <- sif %>%
  filter(gender=='male' & type=='tumor') %>%
  pull(sample.id)

T.autox.norm.rescaled[grepl('X', rownames(T.autox.norm.rescaled)), male.tumor.samples] <- T.autox.norm.rescaled[grepl('X', rownames(T.autox.norm.rescaled)), male.tumor.samples] - 1

cat('Done.\n')



## Tangent on male chrY
cat('\nRunning Tangent on male chrY ...\n')

n.df <- read_delim(n.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

male.normal.samples <- sif %>%
  filter(sample.id %in% colnames(n.df) & gender=='male') %>%
  pull(sample.id)

N.m <- n.df[, male.normal.samples] %>%
  as.matrix()

T.m <- t.df[, male.tumor.samples] %>%
  as.matrix()

## Get the origin in the normal subspace
N.m.means <- apply(N.m, 1, mean)
N.m0 <- N.m - N.m.means
T.m0 <- T.m - N.m.means

## Run Tangent
Npi.m <- pracma::pinv(N.m0)
weights.m <- Npi.m %*% T.m0
proj.m <- N.m0 %*% weights.m
T.m.norm <- T.m0 - proj.m

## Re-scaling after Tangent by median
T.m.norm.medians <- T.m.norm[!grepl('^X|^Y', rownames(T.m.norm)),] %>%
  apply(., 2, median)

T.m.norm.rescaled <- t(t(T.m.norm)- T.m.norm.medians)

T.m.y <- T.m.norm.rescaled[grepl('^Y', rownames(T.m.norm.rescaled)), ]

## Adjust chrY so that it is relative to CN=2
T.m.y.adj <- T.m.y - 1



## Combine "autosomes & chrX" and "chrY"
## Replace NAs in female chrY for downstream analysis (e.g. Circular Binary Segmentation)
T.norm <- as.data.frame(T.autox.norm.rescaled) %>%
  bind_rows(as.data.frame(T.m.y.adj))

rownames(T.norm) <- rownames(t.df)
T.norm <- T.norm %>%
  as.data.frame() %>%
  rownames_to_column('locus')

write.table(T.norm, file=here('output/TangentXY', paste0('TangentXYnormalized_tumor_log2RCN_', num.lf, 'latentFactors.txt')), sep='\t', row.names=FALSE, quote=FALSE)

cat('Done.\n')
