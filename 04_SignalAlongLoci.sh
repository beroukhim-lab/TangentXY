#!/bin/bash

# ------------------------------------------------------------------------------
# Author: Kei Enomoto
# Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
# Description: Code to visualize genomic copy number signal intensity along
#              genomic loci.
# License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
# ------------------------------------------------------------------------------

DefaultWorkingDir=$(pwd)

usage() {
	echo "usage: 04_SignalNoise.sh -d <WorkingDir> -r <PIF> -t <TumorMatrix> -i <SampleOfInterest> -l <LatentFactorNum>"

	echo "  WorkingDir: Working directory. 'output' directory is automatically generated under this directory."
	echo "  PIF: Probe information file."
	echo "  TumorMatrix: Tumor sample signal matrix file with the values in log2(Relative Copy Number) format."
	echo "  SampleOfInterest: Sample of interest for which you want to visualize signals."
	echo "  LatentFactorNum: The number of latent factors to reconstruct a normal subspace. Need to be less than or equal to the number of normal samples in the normal sample signal matrix."
}

while getopts :d:r:t:i:l: option; do
	case "${option}" in
		d) WorkingDir=${OPTARG};;
		r) PIF=${OPTARG};;
		t) TumorsMatrix=${OPTARG};;
		i) SampleOfInterest=${OPTARG};;
		l) LatentFactorNum=${OPTARG};;
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

if [ -z ${PIF} ] || [ -z ${TumorsMatrix} ] || [ -z ${SampleOfInterest} ] || [ -z ${LatentFactorNum} ]; then
	echo "Key arguments (Sample information file, Tumor samples matrix, Number of latentfactors) not present."
	usage
	exit 1
else
	PIF_Absolute=$(realpath ${PIF})
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

echo -e "\nVisualizing signal of ${SampleOfInterest}...\n"

Rscript --slave ./module/05_SignalAlongLoci.R \
	-d ${WorkingDir_Absolute}  -r ${PIF_Absolute} -t ${TumorsMatrix_Absolute} -i ${SampleOfInterest} -l ${LatentFactorNum}
