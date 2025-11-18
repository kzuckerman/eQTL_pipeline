### Introduction
Pipeline for performing an eQTL analysis on single cell + genotype data, with 
the option to perform summary-data-based mendelian randomization using SMR

### Installation

Install and activate conda environment from env.yml

```
conda env create -f env.yml
conda activate eqtl
```
Note: if using MacOS, remove plink2 from the env.yml file and install PLINK2 otherwise

Download + unpack SMR and OSCA for your OS from 
https://yanglab.westlake.edu.cn/software/smr/#Overview, add both to PATH 


### Usage
```
./bash_eqtl.sh PLINK_OBJECT SCRNA_OBJECT CLUSTER_NAME COVARIATES WORKDIR [GWAS_SUM_STATS] [NUM_CPUS]
```
PLINK_OBJECT: Path to PLINK genotype object (omitting any file extension)

SCRNA_OBJECT: Path to h5ad object with gene expression data (QCed ideally)

CLUSTER_NAME: Name of the column in obs of SCRNA_OBJECT to use as cluster labels, used as covariates downstream

COVARIATES: Path to CSV containing covariates. Should contain at minimum "Sample", "Age", "Sex" columns

WORKDIR: Path to directory to save all files to

[GWAS_SUM_STATS]: Optional. Path to GWAS Summary statistics file (in COJO format) to analyze with eQTL data using SMR

[NUM_CPUS]: Optional. Number of CPUs to use. Default 1


### Notes
- All paths can be gcloud paths or local paths. If GCloud paths are passed, files will be copied to the WORKDIR before proceeding
- Make sure Sample names align between PLINK_OBJECT, SCRNA_OBJECT, and COVARIATES files. Nothing will work if they do not
- If running SMR, note that GWAS summary statistics file MUST be in  GCTA-COJO format (although n column can be NA as it is not access in SMR). Use finngen_to_COJO.sh for converting Finngen data, or see SMR docs for file format
