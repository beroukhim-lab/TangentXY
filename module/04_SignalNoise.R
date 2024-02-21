# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Code to visualize and evaluate noise reduction performance of
#              TangentXY.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

if (!requireNamespace('optparse', quietly=TRUE)) {
  print('Package "optparse" not installed. Installing...')
  install.packages('optparse')
}
library('optparse')


option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-s", "--sif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-r", "--pif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-t", "--t.matrix"), type="character", help="tumor samples signal matrix file", metavar="character"),
  make_option(c("-l", "--latent.n"), type="character", help="number of latent factors to reconstruct normal subspace", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# print(opt)


if (!requireNamespace('tidyverse', quietly=TRUE)) {
  print('Package "tidyverse" not installed. Installing...')
  install.packages('tidyverse')
}
suppressPackageStartupMessages(library('tidyverse'))
library('tidyverse')

dir <- opt$directory
setwd(dir)

if (!requireNamespace('here', quietly=TRUE)) {
  print('Package "here" not installed. Installing...')
  install.packages('here')
}
suppressPackageStartupMessages(library('here'))
library('here')

if (!file.exists(here('output/SignalNoise'))) {
  dir.create(here('output/SignalNoise'), recursive=TRUE)
}


sif.file <- opt$sif
pif.file <- opt$pif
t.df.file <- opt$t.matrix
num.lf <- opt$latent.n %>%
  strsplit(split=',') %>%
  unlist() %>%
  as.numeric() %>%
  sort()

sif <- read_delim(sif.file, progress=FALSE, show_col_types=FALSE)
pif <- read_delim(pif.file, progress=FALSE, show_col_types=FALSE)
t.df <- read_delim(t.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

male.tumors <- sif %>%
  filter(sample.id %in% colnames(t.df) & gender=='male') %>%
  pull(sample.id)

## Signal and noise
## Make function for calculating signal and noise
calc.signal.noise <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(pif)

  signal.auto <- data %>%
    filter(!chr %in% c('X', 'Y')) %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  noise.auto <- data %>% 
    filter(!chr %in% c('X', 'Y')) %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.x <- data %>% 
    filter(chr=='X') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.y <- data %>%
    filter(chr=='Y') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.auto=signal.auto,
                          noise.auto=noise.auto,
                          sn.auto=signal.auto/noise.auto,
                          noise.chrx=noise.x,
                          noise.chry=noise.y)

  return(result.df)
}

options(dplyr.summarise.inform = FALSE)

## Pre-normalization signal
cat('Calculating signal and noise in pre-normalization data...\n')

t0.signal.noise.list <- parallel::mclapply(t.df, calc.signal.noise, mc.cores=parallel::detectCores()-2)
t0.signal.noise.df <- t0.signal.noise.list %>%
  bind_rows() %>%
  mutate(sample.id=names(t0.signal.noise.list)) %>%
  mutate(lf='Pre-normalization')

## After TangentXY
for (i in 1:length(num.lf)) {
  num.lf.i <- num.lf[i]
  if (num.lf.i==1) {
    cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factor...\n'))
  } else {
    cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factors...\n'))
  }

  t.df.i.file <- here('output/TangentXY', paste0('TangentXYnormalized_tumor_log2RCN_', num.lf.i, 'latentFactors.txt'))
  t.df.i <- read_delim(t.df.i.file, , progress=FALSE, show_col_types=FALSE) %>%
    column_to_rownames('locus')

  t.i.signal.noise.list <- parallel::mclapply(t.df.i, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  t.i.signal.noise.df <- t.i.signal.noise.list %>%
    bind_rows() %>%
    mutate(sample.id=names(t.i.signal.noise.list)) %>%
    mutate(lf=as.character(num.lf.i))

  if (i==1) {
    signal.noise.df <- t0.signal.noise.df %>% bind_rows(t.i.signal.noise.df)
  } else {
    signal.noise.df <- signal.noise.df %>% bind_rows(t.i.signal.noise.df)
  }
}
write.table(signal.noise.df, file=here('output/SignalNoise', 'SignalNoise.txt'), sep='\t', row.names=FALSE, quote=FALSE)


signal.noise.df.l <- signal.noise.df %>%
  pivot_longer(cols=matches('^signal\\.|^noise\\.|^sn\\.'), names_pattern='(^signal\\.|^noise\\.|^sn\\.)(.*)$', names_to=c('metric', 'chr')) %>%
  mutate(lf=case_when(chr=='chry' & lf!='Pre-normalization' ~ 'TangentXY', TRUE ~ lf)) %>%
  distinct(sample.id, lf, metric, chr, .keep_all=TRUE) %>%
  mutate(lf=factor(.$lf, levels=c('Pre-normalization', num.lf, 'TangentXY'))) %>%
  mutate(metric=case_when(metric=='signal.' ~ 'Signal',
                          metric=='noise.' ~ 'Noise',
                          metric=='sn.' ~ 'SN')) %>%
  mutate(metric=factor(.$metric, levels=c('Signal', 'Noise', 'SN'))) %>%
  mutate(chr=case_when(chr=='auto' ~ 'Auto',
                        chr=='chrx' ~ 'ChrX',
                        chr=='chry' ~ 'ChrY')) %>%
  mutate(chr=factor(.$chr, levels=c('Auto', 'ChrX', 'ChrY'))) %>%
  left_join(sif, by='sample.id') %>%
  filter(!(chr=='ChrY' & gender=='female')) %>%
  filter(!is.na(value)) %>%
  mutate(gender=str_to_title(gender))

g <- ggplot(signal.noise.df.l, aes(x=lf, y=value)) +
  geom_violin() +
  geom_point(aes(col=gender), position=ggbeeswarm::position_beeswarm()) +
  ylim(0, NA) +
  ggh4x::facet_grid2(metric ~ chr, scales='free', space='free_x', independent='y') +
  labs(y='Value', col='Gender') +
  theme_bw(base_size=30) +
  theme(axis.text.x=element_text(angle=45, vjust =1, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/SignalNoise', 'SignalNoise.pdf'), width=24, height=20)
