# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Part of the TangentXY tool for normalizing and converting read 
#              depths from WES/SNP array data to log2 relative copy values.
#              Code to perform singular value decomposition (SVD) on autosomal
#              and chrX signal matrix in normal samples to reconstruct 
#              dimension-reduced signal matrices.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

library('optparse')

option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-s", "--sif"), type="character", help="sample information file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir <- opt$directory
setwd(dir)

suppressPackageStartupMessages(library('tidyverse'))
library('tidyverse')

suppressPackageStartupMessages(library('here'))
library('here')

if (!file.exists(here('output/SVD'))) {
  dir.create(here('output/SVD'), recursive=TRUE)
}


sif.file <- opt$sif
n.lt.df.file <- here('output/LinearTransformation', 'normal_log2RCN_linearTrans.txt')

sif <- read.delim(sif.file)
n.lt.df <- read_delim(n.lt.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

n.autox <- n.lt.df[!grepl('^Y', rownames(n.lt.df)),] %>%
  as.matrix()

## SVD
n.autox.svd <- svd(n.autox)
colnames(n.autox.svd$u) <- colnames(n.autox)
rownames(n.autox.svd$u) <- rownames(n.autox)
rownames(n.autox.svd$v) <- colnames(n.autox)
saveRDS(n.autox.svd, file=here('output/SVD', 'n.autox.svd.rds'), compress=FALSE)

d.df <- n.autox.svd$d %>%
  as.data.frame() %>%
  setNames('d') %>%
  mutate(n=1:n())

g <- ggplot(d.df, aes(x=n, y=d)) +
  geom_line() +
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  labs(title='Importance of latent factors', y='r (Importance)', col='# of normals') +
  theme_bw(base_size=30) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())
ggsave(g, file=here('output/SVD', 'LatentFactors_Importance.pdf'), width=16, height=8)

if (nrow(d.df) > 100) {
  g <- ggplot(d.df %>% head(100), aes(x=n, y=d)) +
    geom_line() +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(title='Importance of top 100 latent factors', y='r (Importance)', col='# of normals') +
    ggtitle('Importance of top 100 latent factors') +
    theme_bw(base_size=30) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x=element_blank())
  ggsave(g, file=here('output/SVD', 'LatentFactors_Importance_top100.pdf'), width=16, height=8)
}
