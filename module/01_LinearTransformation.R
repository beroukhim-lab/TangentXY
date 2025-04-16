# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Part of the TangentXY tool for normalizing and converting read 
#              depths from WES/SNP array data to log2 relative copy values.
#              Code to linearly transform male chrX signal in normal samples to 
#              a distribution similar to female chrX.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

library('optparse')

option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-s", "--sif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-n", "--n.matrix"), type="character", help="normal samples signal matrix", metavar="character"),
  make_option(c("-t", "--t.matrix"), type="character", help="tumor samples signal matrix", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir <- opt$directory
setwd(dir)

suppressPackageStartupMessages(library('tidyverse'))
library('tidyverse')

suppressPackageStartupMessages(library('here'))
library('here')

if (!file.exists(here('output/LinearTransformation'))) {
  dir.create(here('output/LinearTransformation'), recursive=TRUE)
}


sif.file <- opt$sif
n.df.file <- opt$n.matrix
t.df.file <- opt$t.matrix

sif <- read.delim(sif.file)
n.df <- read_delim(n.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')
t.df <- read_delim(t.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

## Linear transformation on male chrX in normal samples
female.normal.samples <- sif %>%
  filter(gender=='female' & type=='normal') %>%
  pull(sample.id)
male.normal.samples <- sif %>%
  filter(gender=='male' & type=='normal') %>%
  pull(sample.id)

female.x.mean <- n.df[grepl('^X', rownames(n.df)),] %>%
  select(all_of(female.normal.samples)) %>%
  as.matrix() %>%
  mean()
male.x.mean <- n.df[grepl('^X', rownames(n.df)),] %>%
  select(all_of(male.normal.samples)) %>%
  as.matrix() %>%
  mean()

n.df.x.transformed <- n.df[grepl('^X', rownames(n.df)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
  mutate(signal=case_when(sample.id %in% male.normal.samples ~ signal - male.x.mean + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='sample.id', values_from='signal') %>%
  unite(col=locus, c('chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

## Linear transformation on male chrY in normal samples
male.y.mean <- n.df[grepl('^Y', rownames(n.df)),] %>%
  select(all_of(male.normal.samples)) %>%
  as.matrix() %>%
  mean()

male.y.mean.mode <- n.df[grepl('^Y', rownames(n.df)),] %>%
  select(all_of(male.normal.samples)) %>%
  as.matrix() %>%
  apply(., 2, mean) %>%
  density() %>%
  {.$x[which.max(.$y)]}

n.df.y.transformed <- n.df[grepl('^Y', rownames(n.df)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
  group_by(sample.id) %>%
  mutate(sample_mean=mean(signal), sample_sd=sd(signal)) %>%
  ungroup() %>%
  mutate(signal=case_when(sample.id %in% male.normal.samples ~ signal - sample_mean + male.y.mean.mode,
                                      TRUE ~ signal)) %>%
  select(chr, pos, sample.id, signal) %>%
  pivot_wider(names_from='sample.id', values_from='signal') %>%
  unite(col=locus, c('chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

n.df.transformed <- n.df[!grepl('^X|^Y', rownames(n.df)), ] %>%
  bind_rows(n.df.x.transformed) %>%
  bind_rows(n.df.y.transformed) %>%
  rownames_to_column('locus')
write.table(n.df.transformed, file=here('output/LinearTransformation', 'normal_log2RCN_linearTrans.txt'), sep='\t', row.names=FALSE, quote=FALSE)

## Linear transformation on male chrX in tumor samples
male.tumor.samples <- sif %>%
  filter(gender=='male' & type=='tumor') %>%
  pull(sample.id)

t.df.x.transformed <- t.df[grepl('^X|^Y', rownames(t.df)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
  mutate(signal=case_when(sample.id %in% male.tumor.samples & chr=='X' ~ signal - male.x.mean + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='sample.id', values_from='signal') %>%
  unite(col=locus, c('chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

t.df.transformed <- t.df[!grepl('^X|^Y', rownames(t.df)), ] %>%
  bind_rows(t.df.x.transformed) %>%
  rownames_to_column('locus')
write.table(t.df.transformed, file=here('output/LinearTransformation', 'tumor_log2RCN_linearTrans.txt'), sep='\t', row.names=FALSE, quote=FALSE)



## Check the signal distribution of chrX before and after linear transformation
signaln.n.x.before <- n.df[grepl('^X', rownames(n.df)),] %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='Before transformation') %>%
  mutate(type='Normal')

signal.n.x.after <- n.df.transformed %>%
  filter(grepl('^X', locus)) %>%
  column_to_rownames('locus') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='After transformation') %>%
  mutate(type='Normal')

signal.t.x.before <- t.df[grepl('^X', rownames(t.df)),] %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='Before transformation') %>%
  mutate(type='Tumor')

signal.t.x.after <- t.df.transformed %>%
  filter(grepl('^X', locus)) %>%
  column_to_rownames('locus') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='After transformation') %>%
  mutate(type='Tumor')

signal.x <- signaln.n.x.before %>%
  bind_rows(signal.n.x.after) %>%
  bind_rows(signal.t.x.before) %>%
  bind_rows(signal.t.x.after) %>%
  mutate(transformation=factor(.$transformation, levels=c('Before transformation', 'After transformation'))) %>%
  mutate(gender=str_to_title(gender))

g <- ggplot(signal.x, aes(x=signal, group=sample.id)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_grid(type~transformation) +
  labs(title='Normal samples chrX signal distribution', fill='Gender') +
  theme_bw(base_size=20)
ggsave(g, file=here('output/LinearTransformation', 'ChrX_SignalDistribution.pdf'), width=12, height=8)

## Check the signal distribution of chrY before and after linear transformation
signaln.n.y.before <- n.df[grepl('^Y', rownames(n.df)),] %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='Before transformation') %>%
  mutate(type='Normal')

signal.n.y.after <- n.df.transformed %>%
  filter(grepl('^Y', locus)) %>%
  column_to_rownames('locus') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='After transformation') %>%
  mutate(type='Normal')

signal.t.y.before <- t.df[grepl('^Y', rownames(t.df)),] %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='Before transformation') %>%
  mutate(type='Tumor')

signal.t.y.after <- t.df.transformed %>%
  filter(grepl('^Y', locus)) %>%
  column_to_rownames('locus') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='After transformation') %>%
  mutate(type='Tumor')

signal.y <- signaln.n.y.before %>%
  bind_rows(signal.n.y.after) %>%
  bind_rows(signal.t.y.before) %>%
  bind_rows(signal.t.y.after) %>%
  mutate(transformation=factor(.$transformation, levels=c('Before transformation', 'After transformation'))) %>%
  mutate(gender=str_to_title(gender))

g <- ggplot(signal.y %>% filter(gender=='Male'), aes(x=signal, group=sample.id)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_grid(type~transformation) +
  scale_fill_manual(values=c('Male'='#00BFC4')) +
  labs(title='Normal samples chrY signal distribution', fill='Gender') +
  theme_bw(base_size=20)
ggsave(g, file=here('output/LinearTransformation', 'ChrY_SignalDistribution.pdf'), width=12, height=8)
