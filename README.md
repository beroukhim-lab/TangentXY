# TangentXY

This repository contains the code for running TangentXY, a sex-informed extension of Tangent copy number inference pipeline.

## To run:
### 1. Download this repository to your local drive:
Use the "Download ZIP" button and decompress the downloaded file.

### 2. Prepare input files:
The sample input files are provided only for formatting references.
Prepare your own data files according to these formats.  
- Sample information file: ./sampledata/sample_information.txt
- Probe information file: ./sampledata/probe_annotation.txt
- Normal sample signal matrix file: ./sampledata/normal_log2RCN.txt
- Tumor sample signal matrix file: ./sampledata/tumor_log2RCN.txt
Note: Remove prefix "chr" from chromosome names in probe information file, normal sample signal matrix file, and tumor sample signal matrix file (e.g. chr1 -> 1, chrX -> X).

Note that in order to obtain good normalization results, normal and tumor sample signal matrix need to be quality-controlled based on, but not limited to, the following criteria prior to input into TangentXY:
- Remove likely failed samples (e.g. samples with too many zeros/NAs)
- Remove likely failed loci (e.g. loci with too many zeros/NAs)
- Remove samples likely to be mislabeled by sex (female <- male)
- Remove samples likely to be mislabeled by sample type (normal <-> tumor)
- Replace obvious outliers with a representative values (e.g. marginal mean/median)
- Remove loci corresponding to common germline CNvs  
The scripts for the quality-control are not provided in this repository, however, Tangent paper describes some of those steps.  
Gao, G. F. et al. Tangent normalization for somatic copy-number inference in cancer genome analysis. Bioinformatics (2022) doi:10.1093/bioinformatics/btac586.
  


### 2. Run:
Follow the steps below.
#### 1. Preprocessing of normal sample signals: Linear transformation of male chrX signal and SVD on autosomes and chrX
Run this tool first to preprocess normal sample signals.
```
## Change directory to the downloaded and decompressed directory
cd ~/Downloads/TangentXY-main
./01_PreprocessNormals.sh \
  -d ./ # Directory to store the outputs \  (Optional. Default: current directory)
  -s ./sampledata/sample_information.txt \  (Mandatory. Sample information file)
  -n ./sampledata/normal_log2RCN.txt \      (Mandatory. Normal sample signal matrix)
```

#### 2. TangentXY with a reconstructed normal subspace with a desired number of latent factors
Run Tangent using the preprocessed normal sample signals and the number of latent factors of your choice.  
A plot generated in the preprocessing might help you select a good number for latent factors.  
((working directory)/output/SVD/LatentFactors_Importance.pdf)
```
./02_TangentXY.sh \
  -d ./ # Directory to store the outputs \  (Optional. Default: current directory)
  -s ./sampledata/sample_information.txt \  (Mandatory. Sample information file)
  -n ./sampledata/normal_log2RCN.txt \      (Mandatory. Normal sample signal matrix)
  -p ./output/SVD/n.autox.svd.RData \       (Mandatory. Output of 01_PreProcessNormals.sh)
  -t ./sampledata/tumor_log2RCN.txt \       (Mandatory. Tumor sample signals to be normalizaed)
  -l 5                                      (Mandatory. The number of latent factors used to reconstruct a normal subspace)
```

#### 3. (Optional) Check the normalization performance by looking at signal level, noise level, and signal-to-noise ratio
You may want to run TangentXY with several different number of latent factors and decide the best number based on the normalization performance. Use this tool for the comparison.
```
./03_SignalNoise.sh \
  -d ./ # Directory to store the outputs \  (Optional. Default: current directory)
  -s ./sampledata/sample_information.txt \  (Mandatory. Sample information file)
  -r ./sampledata/probe_annotation.txt \    (Mandatory. Probe information file)
  -t ./sampledata/tumor_log2RCN.txt \       (Mandatory. Tumor sample signals to be normalizaed)
  -l 1,5,10                                 (Mandatory. The number of latent factors used to reconstruct a normal subspace)
```

#### 4. (Optional) Visualize the signal intensity of sample of interest along genomic position
You may also want to visually check the signal intensity of pre-/post-TangentXY data. This script generates a visualization of signal in a sample of your choice.
```
./04_SignalAlongLoci.sh
  -d ./ # Directory to store the outputs \  (Optional. Default: current directory)
  -r ./sampledata/probe_annotation.txt \    (Mandatory. Probe information file)
  -t ./sampledata/tumor_log2RCN.txt \       (Mandatory. Tumor sample signals to be normalizaed)
  -i tumor.female1 \                        (Mandatory. The name of sample you want to visualize)
  -l 1,5,10                                 (Mandatory. The number of latent factors used to reconstruct a normal subspace)
```

### 3. Output
```
(Working directory)
└── output
    ├── LinearTransformation
    │   ├── NormalSamples_ChrX_ZscoreDistribution.pdf
    │   └── normal_log2RCN_linearTrans.txt
    ├── SVD
    │   ├── LatentFactors_Importance.pdf
    │   └── n.autox.svd.RData
    ├── SignalAlongLoci
    │   └── tumor.female1_signal.pdf
    ├── SignalNoise
    │   ├── SignalNoise.png
    │   └── SignalNoise.txt
    └── TangentXY
        ├── TangentXYnormalized_tumor_log2RCN_10latentFactors.txt
        ├── TangentXYnormalized_tumor_log2RCN_1latentFactors.txt
        └── TangentXYnormalized_tumor_log2RCN_5latentFactors.txt
```
