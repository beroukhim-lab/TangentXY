#!/bin/bash

# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Code to perform linear transformation and singular value
#							 decomposition (SVD) on reference plane.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

DefaultWorkingDir=$(pwd)

usage() {
	echo "usage: 00_InstallRPackages.sh"
}

Rscript --slave ./module/00_InstallRPackages.R
