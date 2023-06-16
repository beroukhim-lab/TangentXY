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

  signal.x <- data %>%
    filter(chr=='X') %>%
    group_by(arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  noise.x <- data %>% 
    filter(chr=='X') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  signal.y <- data %>%
    filter(chr=='Y') %>%
    group_by(arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  noise.y <- data %>%
    filter(chr=='Y') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.auto=signal.auto,
                          noise.auto=noise.auto,
                          sn.auto=signal.auto/noise.auto,
                          signal.chrx=signal.x,
                          noise.chrx=noise.x,
                          sn.chrx=signal.x/noise.x,
                          signal.chry=signal.y,
                          noise.chry=noise.y,
                          sn.chry=signal.y/noise.y)

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
# signal.noise.df <- read.delim(file=here('output/SignalNoise', 'SignalNoise.txt'))


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
ggsave(g, file=here('output/SignalNoise', 'SignalNoise.png'), dpi=100, width=24, height=20)






# ## Signal (scatter plot)
# signal <- Y.signal.noise.df %>%
#   left_join(X.signal.noise.df %>% select(SampleID, signal), by='SampleID') %>%
#   left_join(sif, by='SampleID')

# g <- ggplot(signal %>% filter(gender!='na'), aes(x=signal.y, y=signal.x)) +
#   geom_point(aes(col=project, shape=gender), alpha=0.5) +
#   geom_abline(slope=1, intercept=0, linetype='dashed') +
#   facet_wrap(~dim, nrow=1) +
#   labs(title='Signal', x='Pre-Tangennt', y='Post-Tangent') +
#   coord_fixed(ratio = 1) +
#   guides(shape=guide_legend(order=1), color=guide_legend(order=2)) +
#   theme_bw(base_size=30) +
#   theme(legend.position='bottom')
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter.png'), dpi=100, width=24, height=8)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter.pdf'), width=24, height=8)
# # ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter.svg'), device=svg, width=24, height=8)

# NB.similarity <- readRDS(file=here('../20220725_Tangent_sex_TCGA_WES/output', 'NB.similarity.RData')) %>%
#   mutate(SampleID=paste(ID, tumor.sample, sep='.'))
# NT.similarity <- readRDS(file=here('../20220725_Tangent_sex_TCGA_WES/output', 'NT.similarity.RData')) %>%
#   mutate(SampleID=paste(ID, tumor.sample, sep='.'))
# similarity <- NB.similarity %>%
#   bind_rows(NT.similarity) %>%
#   filter(!is.na(similarity)) %>%
#   arrange(similarity) %>%
#   filter(!duplicated(SampleID))

# g <- ggplot(signal %>%
#       filter(dim==50) %>%
#       filter(gender!='na') %>%
#       left_join(similarity, by='SampleID') %>%
#       filter(!is.na(similarity)),
#     aes(x=signal.y, y=signal.x)) +
#   geom_point(aes(col=similarity, shape=gender)) +
#   geom_abline(slope=1, intercept=0, linetype='dashed') +
#   scale_shape_manual(values=c('female'=1, 'male'=2)) +
#   scale_color_gradientn(colors=rainbow(8), limits=c(0, 8)) +
#   facet_wrap(~project, nrow=5) +
#   labs(title='Signal', x='Pre-Tangennt', y='Post-Tangent') +
#   coord_fixed(ratio = 1) +
#   # guides(shape=guide_legend(order=1), color=guide_legend(order=2)) +
#   theme_bw(base_size=20) +
#   theme(legend.position='bottom')
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter_byTumorType.png'), dpi=100, width=18, height=12)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter_byTumorType.pdf'), width=18, height=12)
# # ggsave(g, file=here('output/12_SignalNoise', 'Signal_scatter_byTumorType.svg'), device=svg, width=18, height=12)


# signal.noise <- X.signal.noise.df %>%
#   bind_rows(Y.signal.noise.df %>% mutate(dim=as.character(dim))) %>%
#   mutate(sn=signal/noise) %>%
#   mutate(dim=factor(.$dim, levels=.$dim %>% unique())) %>%
#   left_join(sif, by='SampleID')

# ## Signal (violin plot)
# g <- ggplot(signal.noise %>% filter(gender!='na'), aes(x=dim, y=signal)) +
#   geom_violin(fill='red') +
#   geom_boxplot(outlier.shape=NA, width=0.05) +
#   # ggbeeswarm::geom_beeswarm(aes(col=project, shape=gender), show.legend=F) +
#   # geom_line(aes(group=SampleID)) +
#   ylim(0, NA) +
#   labs(x='Number of dimensions in reference plane', y='Signal') +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_violin.png'), dpi=100, width=10, height=6)
# ggsave(g, file=here('output/12_SignalNoise', 'Signal_violin.pdf'), width=10, height=6)
# # ggsave(g, file=here('output/12_SignalNoise', 'Signal_violin.svg'), device=svg, width=10, height=6)

# ## Noise (violin plot)
# g <- ggplot(signal.noise %>% filter(gender!='na'), aes(x=dim, y=noise)) +
#   geom_violin(fill='red') +
#   geom_boxplot(outlier.shape=NA, width=0.05) +
#   # ggbeeswarm::geom_beeswarm(aes(col=project, shape=gender), show.legend=F) +
#   ylim(0, NA) +
#   labs(x='Number of dimensions in reference plane', y='Noise') +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_violin.png'), dpi=100, width=10, height=6)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_violin.pdf'), width=10, height=6)
# # ggsave(g, file=here('output/12_SignalNoise', 'Noise_violin.svg'), device=svg, width=10, height=6)

# ## Signal-to-noise (violin plot)
# g <- ggplot(signal.noise %>% filter(gender!='na'), aes(x=dim, y=sn)) +
#   geom_violin(fill='red') +
#   geom_boxplot(outlier.shape=NA, width=0.05) +
#   # ggbeeswarm::geom_beeswarm(aes(col=project, shape=gender), show.legend=F) +
#   labs(x='Number of dimensions in reference plane', y='S/N') +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'SN_violin.png'), dpi=100, width=10, height=6)
# ggsave(g, file=here('output/12_SignalNoise', 'SN_violin.pdf'), width=10, height=6)
# # ggsave(g, file=here('output/12_SignalNoise', 'SN_violin.svg'), device=svg, width=10, height=6)

# dimensions <- signal.noise$dim %>% unique()
# g.signal.list <- vector(mode='list', length=length(dimensions)-1)
# g.noise.list <- vector(mode='list', length=length(dimensions)-1)
# g.sn.list <- vector(mode='list', length=length(dimensions)-1)
# for (i in 1:(length(dimensions)-1)) {
#   dim1 <- dimensions[i]
#   dim2 <- dimensions[i+1]
#   print(paste(dim1, 'and', dim2))

#   signal.i <- signal.noise %>%
#     filter(gender!='na') %>%
#     filter(dim %in% c(dim1, dim2)) %>%
#     select(signal, SampleID, dim, project, type, gender) %>%
#     pivot_wider(names_from=dim, values_from=signal) %>%
#     rename(dim1=5, dim2=6)

#   noise.i <- signal.noise %>%
#     filter(gender!='na') %>%
#     filter(dim %in% c(dim1, dim2)) %>%
#     select(noise, SampleID, dim, project, type, gender) %>%
#     pivot_wider(names_from=dim, values_from=noise) %>%
#     rename(dim1=5, dim2=6)

#   sn.i <- signal.noise %>%
#     filter(gender!='na') %>%
#     filter(dim %in% c(dim1, dim2)) %>%
#     select(sn, SampleID, dim, project, type, gender) %>%
#     pivot_wider(names_from=dim, values_from=sn) %>%
#     rename(dim1=5, dim2=6)

#   if (dim1!='Pre-Tangent') {
#     dim1 <- paste0('dim=', dim1)
#   }
#   dim2 <- paste0('dim=', dim2)

#   g.signal <- ggplot(signal.i, aes(x=dim1, y=dim2)) +
#     geom_point(aes(col=project, shape=gender), show.legend=FALSE) +
#     geom_abline(slope=1, intercept=0) +
#     labs(x=dim1, y=dim2) +
#     theme_bw(base_size=20)
#   # ggsave(g.signal, file=here('output/12_SignalNoise', 'Signal_compareNeighbors.svg'), device=svg, width=10, height=8)

#   g.noise <- ggplot(noise.i, aes(x=dim1, y=dim2)) +
#     geom_point(aes(col=project, shape=gender), show.legend=FALSE) +
#     geom_abline(slope=1, intercept=0) +
#     labs(x=dim1, y=dim2) +
#     theme_bw(base_size=20)

#   g.sn <- ggplot(sn.i, aes(x=dim1, y=dim2)) +
#     geom_point(aes(col=project, shape=gender), show.legend=FALSE) +
#     geom_abline(slope=1, intercept=0) +
#     labs(x=dim1, y=dim2) +
#     theme_bw(base_size=20)

#   g.signal.list[[i]] <- g.signal
#   g.noise.list[[i]] <- g.noise
#   g.sn.list[[i]] <- g.sn
# }

# g.signal.panel <- patchwork::wrap_plots(g.signal.list, nrow=4, ncol=4) + patchwork::plot_annotation('Signal', theme = theme(plot.title=element_text(size=30)))
# g.noise.panel <- patchwork::wrap_plots(g.noise.list, nrow=4, ncol=4) + patchwork::plot_annotation('Noise', theme = theme(plot.title=element_text(size=30)))
# g.sn.panel <- patchwork::wrap_plots(g.sn.list, nrow=4, ncol=4) + patchwork::plot_annotation('Signal-to-Noise', theme = theme(plot.title=element_text(size=30)))
# ggsave(g.signal.panel, file=here('output/12_SignalNoise', 'Signal_compareNeighbors.png'), dpi=100, width=24, height=24)
# ggsave(g.signal.panel, file=here('output/12_SignalNoise', 'Signal_compareNeighbors.pdf'), width=24, height=24)
# ggsave(g.noise.panel, file=here('output/12_SignalNoise', 'Noise_compareNeighbors.png'), dpi=100, width=24, height=24)
# ggsave(g.noise.panel, file=here('output/12_SignalNoise', 'Noise_compareNeighbors.pdf'), width=24, height=24)
# ggsave(g.sn.panel, file=here('output/12_SignalNoise', 'SN_compareNeighbors.png'), dpi=100, width=24, height=24)
# ggsave(g.sn.panel, file=here('output/12_SignalNoise', 'SN_compareNeighbors.pdf'), width=24, height=24)






# ## Noise on each chromosome. Run this code after merging Y.autox and Y.y in 14_combineAndDivideTangentOutput.R
# ## Make function for calculating noise
# calc.noise.each.chr <- function(list) {
#   data <- list %>%
#     as.data.frame() %>%
#     setNames('signal') %>%
#     bind_cols(probes)

#   noise <- data %>%
#     mutate(signal.lag=lag(signal)) %>%
#     group_by(chr, arm) %>%
#     mutate(n=1:n()) %>%
#     filter(n!=1) %>%
#     ungroup() %>%
#     mutate(abs.diff=abs(signal.lag - signal)) %>%
#     group_by(chr) %>%
#     summarize(noise=median(abs.diff)) %>%
#     mutate(chr=factor(.$chr, levels=.$chr %>% gtools::mixedsort())) %>%
#     arrange(chr) %>%
#     column_to_rownames('chr')

#   return(noise)
# }

# dim <- 50 # Dimension to use for normal subspace

# Y.combined <- readRDS(file=here(paste0('output/14_combineAndDivideTangentOutput/rescalingAfterTangent/dimension', dim), 'Y.RData')) %>%
#   as.data.frame()

# Y.combined.noise.list <- parallel::mclapply(Y.combined, calc.noise.each.chr, mc.cores=parallel::detectCores()-2)
# Y.combined.noise.df <- Y.combined.noise.list %>%
#   bind_cols() %>%
#   setNames(names(Y.combined.noise.list))
# saveRDS(Y.combined.noise.df, file=here('output/12_SignalNoise', 'Y.combined.noise.df.RData'), compress=FALSE)
# # Y.combined.noise.df <- readRDS(file=here('output/12_SignalNoise', 'Y.combined.noise.df.RData'))


# noise.df <- Y.combined.noise.df %>%
#   rownames_to_column('chr') %>%
#   pivot_longer(names_to='SampleID', values_to='noise', cols=-'chr') %>%
#   left_join(sif, by='SampleID') %>%
#   filter(gender!='na') %>%
#   filter(!(gender=='female' & chr=='Y')) %>%
#   mutate(chr=factor(.$chr, levels=.$chr %>% unique() %>% gtools::mixedsort()))

# g <- ggplot(noise.df, aes(x=chr, y=noise)) +
#   geom_violin(aes(col=gender), position = position_dodge(0.9)) +
#   geom_boxplot(aes(col=gender), outlier.shape=NA, width=0.1, position = position_dodge(0.9), show.legend=FALSE) +
#   scale_y_continuous(limits=c(0, 1.3), breaks=seq(0, 1.2, by=0.2)) +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_eachChr.png'), dpi=100, width=16, height=6)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_eachChr.pdf'), width=16, height=6)

# noise.df.x <- noise.df %>%
#   filter(chr=='X') %>%
#   group_by(project) %>%
#   arrange(project, noise) %>%
#   ungroup() %>%
#   mutate(SampleID=factor(.$SampleID, levels=.$SampleID))

# noise.threshold <- 0.3
# g <- ggplot(noise.df.x, aes(x=SampleID, y=noise)) +
#   geom_point(aes(col=project, shape=gender), size=3) +
#   geom_hline(yintercept=noise.threshold, linetype='dashed', col='red') +
#   # scale_x_discrete(expand=expansion(add=10)) +
#   scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2)) +
#   facet_wrap(~project, nrow=1, scales='free_x') +
#   coord_cartesian(clip='off') +
#   labs(title='Noise in chrX after ZSTangent') +
#   theme_bw(base_size=20) +
#   theme(panel.grid.major.x=element_blank(), axis.text.x=element_blank()) +
#   theme(strip.text.x=element_text(angle=90))
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX.png'), dpi=100, width=24, height=8)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX.pdf'), width=24, height=8)

# g <- ggplot(noise.df.x, aes(x=noise)) +
#   geom_density(aes(fill=project), alpha=0.2) +
#   geom_vline(xintercept=noise.threshold, linetype='dashed', col='red') +
#   scale_x_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2)) +
#   labs(title='Distribution of noise in chrX after ZSTangent') +
#   theme_bw(base_size=20)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX_density.png'), dpi=100, width=16, height=8)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX_density.pdf'), width=16, height=8)

# noise.autosomes <- noise.df %>%
#   filter(!chr %in% c('X', 'Y')) %>%
#   group_by(SampleID) %>%
#   summarize(noise.auto=median(noise))

# noise.x.auto <- noise.df.x %>%
#   left_join(noise.autosomes, by='SampleID')

# g <- ggplot(noise.x.auto, aes(x=noise, y=noise.auto)) +
#   geom_abline(intercept=0, slope=1, linetype='dashed') +
#   geom_point(aes(col=project, shape=gender)) +
#   scale_x_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2)) +
#   scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2)) +
#   coord_fixed(ratio=1) +
#   # facet_wrap(~project, nrow=1, scales='free_x') +
#   labs(title='Noise in chrX vs autosomes after ZSTangent', x='Noise in chrX', y='Noise in autosomes') +
#   theme_bw(base_size=20)
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrXvsAutosomes.png'), dpi=100, width=12, height=8)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrXvsAutosomes.pdf'), width=12, height=8)
  
# noise.x.pre.post <- noise.df.x %>%
#   left_join(X.signal.noise.df %>% select(SampleID, noise.chrx), by='SampleID')

# g <- ggplot(noise.x.pre.post, aes(x=noise.chrx, y=noise)) +
#   geom_abline(intercept=0, slope=1, linetype='dashed') +
#   geom_hline(yintercept=noise.threshold, linetype='dashed', col='red') +
#   geom_point(aes(col=project, shape=gender)) +
#   scale_x_continuous(limits=c(0, NA), breaks=seq(0, 1.4, by=0.2)) +
#   scale_y_continuous(limits=c(0, NA), breaks=seq(0, 1.0, by=0.2)) +
#   coord_fixed(ratio=1) +
#   # facet_wrap(~project, nrow=1, scales='free_x') +
#   labs(title='Noise in chrX (Pre-ZSTangent vs Post-ZSTangent)', x='Pre-ZSTangent noise', y='Post-ZSTangent noise') +
#   theme_bw(base_size=20)
# # show(g)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX_PrePost.png'), dpi=100, width=12, height=8)
# ggsave(g, file=here('output/12_SignalNoise', 'Noise_chrX_PrePost.pdf'), width=12, height=8)


# ## Figures for paper
# dimensions <- c('Pre-Tangent', 10, 30, 50, 100, 200, 500, 5000, 10441)

# ## Noise (violin plot)
# g <- ggplot(signal.noise %>%
#       filter(gender!='na') %>%
#       filter(dim %in% dimensions),
#     aes(x=dim, y=noise)) +
#   geom_violin(fill='red') +
#   geom_boxplot(outlier.shape=NA, width=0.1) +
#   ylim(0, NA) +
#   labs(x='Number of dimensions in reference plane', y='Noise') +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave(g, file=here('output/12_SignalNoise/FiguresForPaper', 'Noise_violin.pdf'), width=9, height=6)

# ## Signal-to-noise (violin plot)
# g <- ggplot(signal.noise %>%
#       filter(gender!='na') %>%
#       filter(dim %in% dimensions),
#     aes(x=dim, y=sn)) +
#   geom_violin(fill='red') +
#   geom_boxplot(outlier.shape=NA, width=0.1) +
#   labs(x='Number of dimensions in reference plane', y='S/N') +
#   theme_bw(base_size=30) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave(g, file=here('output/12_SignalNoise/FiguresForPaper', 'SN_violin.pdf'), width=9, height=6)
