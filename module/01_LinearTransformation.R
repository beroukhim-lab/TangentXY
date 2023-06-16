if (!requireNamespace('optparse', quietly=TRUE)) {
  print('Package "optparse" not installed. Installing...')
  install.packages('optparse')
}
library('optparse')


option_list = list(
  make_option(c("-d", "--directory"), type="character", help="working directory", metavar='character'),
  make_option(c("-s", "--sif"), type="character", help="sample information file", metavar="character"),
  make_option(c("-n", "--n.matrix"), type="character", help="normal samples signal matrix", metavar="character")
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

if (!file.exists(here('output/LinearTransformation'))) {
  dir.create(here('output/LinearTransformation'), recursive=TRUE)
}

sif.file <- opt$sif
n.df.file <- opt$n.matrix

sif <- read_delim(sif.file, progress=FALSE, show_col_types=FALSE)
n.df <- read_delim(n.df.file, progress=FALSE, show_col_types=FALSE) %>%
  column_to_rownames('locus')

## Linear transformation only on male chrX in normal samples
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
female.x.sd <- n.df[grepl('^X', rownames(n.df)),] %>%
  select(all_of(female.normal.samples)) %>%
  as.matrix() %>%
  sd()
male.x.sd <- n.df[grepl('X', rownames(n.df)),] %>%
  select(all_of(male.normal.samples)) %>%
  as.matrix() %>%
  sd()

n.df.x.transformed <- n.df[grepl('^X|^Y', rownames(n.df)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
  mutate(signal=case_when(sample.id %in% male.normal.samples & chr=='X' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='sample.id', values_from='signal') %>%
  unite(col=locus, c('chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

n.df.transformed <- n.df[!grepl('^X|^Y', rownames(n.df)), ] %>%
  bind_rows(n.df.x.transformed) %>%
  rownames_to_column('locus')
write.table(n.df.transformed, file=here('output/LinearTransformation', 'normal_log2RCN_linearTrans.txt'), sep='\t', row.names=FALSE, quote=FALSE)


## Check the signal distribution of chrX before and after linear transformation
signal.x.before <- n.df[grepl('^X', rownames(n.df)),] %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='Before transformation')

signal.x.after <- n.df.transformed %>%
  filter(grepl('^X', locus)) %>%
  column_to_rownames('locus') %>%
  pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
  left_join(sif, by='sample.id') %>%
  mutate(transformation='After transformation')

signal.x <- bind_rows(signal.x.before, signal.x.after) %>%
  mutate(transformation=factor(.$transformation, levels=c('Before transformation', 'After transformation'))) %>%
  mutate(gender=str_to_title(gender))

g <- ggplot(signal.x, aes(x=signal, group=sample.id)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_wrap(~transformation, nrow=1) +
  labs(title='Normal samples chrX signal distribution', fill='Gender') +
  theme_bw(base_size=20)
ggsave(g, file=here('output/LinearTransformation', 'NormalSamples_ChrX_ZscoreDistribution.pdf'), width=12, height=8)
