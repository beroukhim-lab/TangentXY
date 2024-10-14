# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Part of the TangentXY tool for normalizing and converting read 
#              depths from WES/SNP array data to log2 relative copy values.
#              Code to install required R packages prior to running TangentXY.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

packages.list <- c('optparse', 'tidyverse', 'here', 'scales', 'pracma', 'ggbeeswarm', 'ggh4x', 'gtools')

for (i in seq_along(packages.list)) {
  package.i <- packages.list[i]

  if (!requireNamespace(package.i, quietly=TRUE)) {
    print(paste0('Package "', package.i, '" not installed. Installing...'))
    install.packages(package.i, dependencies=TRUE)      
  }
}

installed.packages <- installed.packages()
not.installed.packages <- setdiff(packages.list, installed.packages)

if (length(not.installed.packages)==0) {
  print('All required packages are installed.')
} else {
  print(paste('Packages below were not installed. Try installing manually.', paste(not.installed.packages, collapse=', '), sep="\n"))
}
