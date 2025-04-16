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
	echo "usage: 01_PreprocessNormals.sh -d <WorkingDir> -s <SIF> -n <NormalsMatrix> -t <TumorsMatrix>"

	echo "  WorkingDir: Working directory. 'output' directory is automatically generated under this directory."
	echo "  SIF: Sample information file."
	echo "  NormalsMatrix: Normal sample signal matrix file with the values in log2(Relative Copy Number) format."
	echo "  TumorsMatrix: Tumor sample signal matrix file with the values in log2(Relative Copy Number) format."
}

##### -------------------- Step 1 (Linear transformation of male chrX signal) -------------------- #####
while getopts :d:s:n:t: option; do
	case "${option}" in
		d) WorkingDir=${OPTARG};;
		s) SIF=${OPTARG};;
		n) NormalsMatrix=${OPTARG};;
		t) TumorsMatrix=${OPTARG};;
		\?)
				echo "Unknown options:"
				usage
				exit 1;;
		:)
				echo "Missing required options:"
				usage
				exit 1;;
		h|*)
				usage
				exit 1;;
	esac
done

if [ -z ${SIF} ] || [ -z ${NormalsMatrix} ] || [ -z ${TumorsMatrix} ]; then
	echo "Key arguments (Sample information file, Normal samples matrix, Tumor samples matrix) not present."
	echo "Exiting..."
	exit 1
else
	SIF_Absolute=$(realpath ${SIF})
	NormalsMatrix_Absolute=$(realpath ${NormalsMatrix})
	TumorsMatrix_Absolute=$(realpath ${TumorsMatrix})
fi

if [ -z ${WorkingDir} ]; then
	WorkingDir=${DefaultWorkingDir}
	WorkingDir_Absolute=$(realpath ${WorkingDir})
	echo "Working directory not set. Using "${WorkingDir_Absolute}" instead."
else
	WorkingDir_Absolute=$(realpath ${WorkingDir})
	echo "Working directory = "${WorkingDir_Absolute}
fi


echo -e "\nRunning linear transformation...\n"

Rscript --slave ./module/01_LinearTransformation.R \
	-d ${WorkingDir_Absolute} -s ${SIF_Absolute} -n ${NormalsMatrix_Absolute} -t ${TumorsMatrix_Absolute}


##### -------------------- Step 2 (SVD on autosomes & chrX) -------------------- #####
echo -e "\nRunning SVD...\n"

Rscript --slave ./module/02_SVD.R \
	-d ${WorkingDir_Absolute} -s ${SIF_Absolute}




