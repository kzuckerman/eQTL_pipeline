#!/bin/bash
# $1 path to plink genotype object NO EXTENSION, gcloud (will cp) or local
# $2 path to h5ad file to analyze, gcloud (will cp) or local
# $3 name of anndata .obs cluster column to use
# $4 path to covariates file by sample 
# $5 directory to work in
# $6 [optional] GWAS Summary statistics file path, gcloud (will cp) or local, must be in GCTA-COJO format
# $7 [optional] number of CPUs available

set -euo pipefail

dir=$5
mkdir -p $dir $dir/processed_matrices $dir/osca_input $dir/osca_intermediate

run_dir=$(dirname $(realpath "$0"))

num_cpus=${7:-1}

echo "Loading files..."
if [[ $1 == gs* ]]; then
    gcloud storage cp ${1}'*' $dir
    plink_path=$dir/$(basename "$1")
else
    plink_path=$1
fi

plink_pref=$(basename "$1")

if [[ $2 == gs* ]]; then
    gcloud storage cp ${2}'*' $dir
    h5ad_path=$dir/$(basename "$2")
else
    h5ad_path=$2
fi
h5ad_pref=$(basename -s .h5ad $h5ad_path)

if [[ $4 == gs* ]]; then
    gcloud storage cp ${4}'*' $dir
    cov_file=$dir/$(basename "$4")
else
    cov_file=$4
fi


count=$(wc -l ${plink_path}.fam | awk '{print $1}')

echo "Getting PCAs..."
if [[ $count -gt 5000 ]]; then
    plink2 --bfile ${plink_path} \
            --pca approx 5 \
            --threads $num_cpus \
            --out ${dir}/$plink_pref > /dev/null
else
    plink2 --bfile ${plink_path} \
            --pca 5 \
            --threads $num_cpus \
            --out ${dir}/$plink_pref > /dev/null
fi


echo "Calculating expression and composition matrices..."
python3 $run_dir/expr_comp_matrices.py $dir ${h5ad_path} ${3} > /dev/null


echo "Formatting files for OSCA..."
Rscript $run_dir/OSCA_format.r $plink_path $cov_file $dir $run_dir $h5ad_pref > /dev/null


# Log file names with prefix from final output
progress_log="$dir/osca_intermediate/osca_progress.log"

befile_prefix="$dir/osca_intermediate/${h5ad_pref}_osca"
final_output="$dir/${h5ad_pref}_"

# Step 1: Create BOD file
echo "OSCA: Running gene-expression --make-bod..."
osca --efile $dir/osca_input/Phenotype_${h5ad_pref}_ocsa.txt --gene-expression --make-bod --out $befile_prefix >> ${progress_log} 

# Step 2: Update OPI with user-specified input
echo "OSCA: Updating OPI..."
osca --befile $befile_prefix --update-opi $dir/osca_input/Upprobe_${h5ad_pref}.opi >> ${progress_log} 

# Step 3: Run eQTL analysis
echo "OSCA: Running eQTL analysis..."
osca --eqtl --bfile $plink_path --befile $befile_prefix --cis --cis-wind 2000 \
    --covar $dir/osca_input/cov1_${h5ad_pref}.txt --qcovar $dir/osca_input/cov2_${h5ad_pref}.txt --to-smr --thread-num $num_cpus \
    --out $dir/osca_intermediate/tempeqtl_$h5ad_pref >> ${progress_log} 

# Step 4: Query eQTL summary with user-specified output
echo "OSCA: Querying eQTL summary..."
osca --beqtl-summary $dir/osca_intermediate/tempeqtl_$h5ad_pref  --query 1 --out $dir/${h5ad_pref}_all_eqtl.tsv >> ${progress_log} 


if [[ -n "${6:-}" ]]; then
    echo "Running SMR.."
    if [[ $6 == gs* ]]; then
        gcloud storage cp $6 $dir
        gwas=$dir/$(basename "$6")
    else
        gwas=$6
    fi
    smr --bfile $plink_path --gwas-summary $gwas --beqtl-summary $dir/osca_intermediate/tempeqtl_$h5ad_pref --out $dir/${h5ad_pref} --thread-num $num_cpus >> $dir/osca_intermediate/smr_progress.log 
fi

echo ""
echo ""