#!/bin/bash

DefaultWorkingDir=$(pwd)

usage() {
	echo "usage: 03_SignalNoise.sh -d <WorkingDir> -s <SIF> -r <PIF> -t <TumorMatrix> -l <LatentFactorNum>"

	echo "  WorkingDir: Working directory. 'output' directory is automatically generated under this directory."
	echo "  SIF: Sample information file."
	echo "  PIF: Probe information file."
	echo "  TumorMatrix: Tumor sample signal matrix file with the values in log2(Relative Copy Number) format."
	echo "  LatentFactorNum: The number of latent factors to reconstruct a normal subspace. Need to be less than or equal to the number of normal samples in the normal sample signal matrix."
}

while getopts :d:s:r:t:l: option; do
	case "${option}" in
		d) WorkingDir=${OPTARG};;
		s) SIF=${OPTARG};;
		r) PIF=${OPTARG};;
		t) TumorsMatrix=${OPTARG};;
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

if [ -z ${SIF} ] || [ -z ${PIF} ] || [ -z ${TumorsMatrix} ] || [ -z ${LatentFactorNum} ]; then
	echo "Key arguments (Sample information file, Tumor samples matrix, Number of latentfactors) not present."
	usage
	exit 1
else
	SIF_Absolute=$(realpath ${SIF})
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


Rscript --slave ./module/04_SignalNoise.R \
	-d ${WorkingDir_Absolute} -s ${SIF_Absolute} -r ${PIF_Absolute} -t ${TumorsMatrix_Absolute} -l ${LatentFactorNum}

