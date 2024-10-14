# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Code to visualize genomic copy number signal intensity along
#              genomic loci.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

library('optparse')

option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-r", "--pif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-t", "--t.matrix"), type="character", help="tumor samples signal matrix file", metavar="character"),
  make_option(c("-i", "--soi"), type="character", help="sample of interest", metavar="character"),
  make_option(c("-l", "--latent.n"), type="character", help="number of latent factors to reconstruct normal subspace", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir <- opt$directory
setwd(dir)

suppressPackageStartupMessages(library('tidyverse'))
library('tidyverse')

suppressPackageStartupMessages(library('here'))
library('here')

if (!file.exists(here('output/SignalAlongLoci'))) {
  dir.create(here('output/SignalAlongLoci'), recursive=TRUE)
}

options(dplyr.summarise.inform = FALSE)


pif.file <- opt$pif
t.df.file <- opt$t.matrix
sample.of.interest <- opt$soi
num.lf <- opt$latent.n %>%
  strsplit(split=',') %>%
  unlist() %>%
  as.numeric() %>%
  sort()

pif <- read_delim(pif.file, progress=FALSE, show_col_types=FALSE)
t.df <- read_delim(t.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

t0 <- read_delim(opt$t.matrix, progress=FALSE, show_col_types=FALSE) %>%
    column_to_rownames('locus') %>%
    select(all_of(sample.of.interest)) %>%
    setNames('signal') %>%
    rownames_to_column('locus') %>%
    mutate(index=1:n()) %>%
    mutate(ln='Pre-normalization')

for (i in 1:length(num.lf)) {
  num.lf.i <- num.lf[i]
  t <- read_delim(here('output/TangentXY', paste0('TangentXYnormalized_tumor_log2RCN_', num.lf.i, 'latentFactors.txt')), progress=FALSE, show_col_types=FALSE) %>%
    column_to_rownames('locus') %>%
    select(all_of(sample.of.interest)) %>%
    setNames('signal') %>%
    rownames_to_column('locus') %>%
    mutate(index=1:n()) %>%
    mutate(ln=as.character(num.lf.i))
  if (i==1) {
    t.df <- bind_rows(t0, t)
  } else {
    t.df <- bind_rows(t.df, t)
  }
}


data.soi <- t.df %>%
  filter(!is.na(signal)) %>%
  left_join(pif, by='locus') %>%
  mutate(chr.class=case_when(chr %in% c(1,3,5,7,9,11,13,15,17,19,21,'X') ~ 'odd', TRUE ~ 'even')) %>%
  group_by(chr, arm, ln) %>%
  mutate(median=median(signal)) %>%
  ungroup() %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique() %>% gtools::mixedsort())) %>%
  arrange(chr, start) %>%
  mutate(ln=factor(.$ln, levels=c('Pre-normalization', as.character(num.lf))))

g <- ggplot(data.soi, aes(x=index, y=signal)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(aes(col=chr.class), show.legend=FALSE) +
  geom_line(aes(y=median), col='red') +
  scale_color_manual(values=c('odd'='black', 'even'='green')) +
  facet_wrap(~ln, nrow=1) +
  labs(y='Signal') +
  theme_bw(base_size=20)
ggsave(g, file=here('output/SignalAlongLoci', paste0(sample.of.interest, '_signal.pdf')), width=8 * (length(num.lf) + 1), height=2 * (length(num.lf) + 1))
