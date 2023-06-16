#!/bin/bash

DefaultWorkingDir=$(pwd)

usage() {
	echo "usage: 01_PreprocessNormals.sh -d <WorkingDir> -s <SIF> -n <NormalsMatrix>"
}

##### -------------------- Step 1 (Linear transformation of male chrX signal) -------------------- #####
while getopts :d:s:n: option; do
	case "${option}" in
		d) WorkingDir=${OPTARG};;
		s) SIF=${OPTARG};;
		n) NormalsMatrix=${OPTARG};;
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

if [ -z ${SIF} ] || [ -z ${NormalsMatrix} ]; then
	echo "Key arguments (Sample information file, ) not present."
	echo "Exiting..."
	exit 1
else
	SIF_Absolute=$(realpath ${SIF})
	NormalsMatrix_Absolute=$(realpath ${NormalsMatrix})
fi

if [ -z ${WorkingDir} ]; then
	echo "Working directory not set. Using "${DefaultWorkingDir}" instead."
	WorkingDir=${DefaultWorkingDir}
else
	WorkingDir_Absolute=$(realpath ${WorkingDir})
	echo "Working directory = "${WorkingDir_Absolute}
fi


echo -e "\nRunning linear transformation...\n"

Rscript --slave ./module/01_LinearTransformation.R \
	-d ${WorkingDir_Absolute} -s ${SIF_Absolute} -n ${NormalsMatrix_Absolute}


##### -------------------- Step 2 (SVD on autosomes & chrX) -------------------- #####
echo -e "\nRunning SVD...\n"

Rscript --slave ./module/02_SVD.R \
	-d ${WorkingDir_Absolute} -s ${SIF_Absolute}




