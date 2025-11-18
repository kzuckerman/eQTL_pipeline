#!/bin/bash
# $1 path to plink genotype object, gcloud (will cp) or local
# $2 path to h5ad file to analyze, gcloud (will cp) or local
# $3 name of anndata .obs cluster column to use
# $4 path to covariates file by sample 
# $5 directory to work in
# $6 [optional] GWAS Summary statistics file path, gcloud (will cp) or local, must be in GCTA-COJO format
# $7 [optional] number of CPUs available

set -euo pipefail

dir=$5
mkdir -p $dir

run_dir=$(basename $(realpath "$0"))

num_cpus=${7:-1}

echo "Loading files..."
if [[ $1 == gs* ]]; then
    gcloud storage cp ${1}'*' $dir
    plink_pref=$dir/$(basename "$1")
else
    plink_pref=$1
fi

if [[ $2 == gs* ]]; then
    gcloud storage cp ${2}'*' $dir
    h5ad_path=$dir/$(basename "$2")
else
    h5ad_path=$2
fi

if [[ $4 == gs* ]]; then
    gcloud storage cp ${4}'*' $dir
    cov_file=$dir/$(basename "$4")
else
    cov_file=$4
fi


count=$(wc -l ${plink_pref}.fam | awk '{print $1}')

echo "Getting PCAs..."
if [[ $count -gt 5000 ]]; then
    plink2 --bfile ${plink_pref} \
            --pca approx 5 \
            --out ${plink_pref} > /dev/null
else
    plink2 --bfile ${plink_pref} \
            --pca 5 \
            --out ${plink_pref} > /dev/null
fi


echo "Calculating expression and composition matrices..."
python3 $run_dir/expr_comp_matrices.py $dir ${h5ad_path} ${3} > /dev/null


echo "Formatting files for OSCA..."
Rscript $run_dir/OSCA_format.r $plink_pref $cov_file $dir $run_dir > /dev/null


# Log file names with prefix from final output
progress_log="$dir/osca_progress.log"

h5ad_pref=$(basename -s .h5ad $h5ad_path)
befile_prefix="$dir/${h5ad_pref}_osca"
final_output="$dir/${h5ad_pref}_final"

# Step 1: Create BOD file
echo "OSCA: Running gene-expression --make-bod..."
osca --efile $dir/Phenotype_${h5ad_pref}_ocsa.txt --gene-expression --make-bod --out $befile_prefix >> ${progress_log} 2>&1

# Step 2: Update OPI with user-specified input
echo "OSCA: Updating OPI..."
osca --befile $befile_prefix --update-opi $dir/Upprobe_${h5ad_pref}.opi >> ${progress_log} 2>&1

# Step 3: Run eQTL analysis
echo "OSCA: Running eQTL analysis..."
osca --eqtl --bfile $plink_pref --befile $befile_prefix --cis --cis-wind 2000 \
    --covar $dir/cov1_${h5ad_pref}.txt --qcovar $dir/cov2_${h5ad_pref}.txt --to-smr --thread-num 4 \
    --out $dir/tempeqtl_$h5ad_pref >> ${progress_log} 2>&1

# Step 4: Query eQTL summary with user-specified output
echo "OSCA: Querying eQTL summary..."
osca --beqtl-summary $dir/tempeqtl_$h5ad_pref  --query 1 --out $final_output >> ${progress_log} 2>&1


if [[ -n "${6:-}" ]]; then
    echo "Running SMR..."
    smr --bfile $plink_pref --gwas-summary $6 --beqtl-summary $dir/tempeqtl_$h5ad_pref --out ${final_output}_smr --thread-num $num_cpus >> $dir/smr_progress.log 2>&1
fi