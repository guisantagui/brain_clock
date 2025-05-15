#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=perts_preproc_sub_alt
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_perts_preproc_sub_alt.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_perts_preproc_sub_alt.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# In this pre-processing pipeline of the merged perturbation data, we first apply the distribution
# of the reference dataset, and then perform SVA on a subsample of the samples to infer SVs in the rest.

# Variables for the pipeline
########################################################################################################################
input="/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/merged_perts_exprsn_qnorm.csv"
input="../../results/parsed/merged_perts/merged_perts_exprsn.csv"
clin="../../results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
train_test="../../results/parsed/merged/train_test.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/merged_perts_metdat.csv"
nSV_method="leek"
subsamp=0.2
nPCs=20
outDir="../../results/preproc/perts_preproc_subsamp_alt/"

mkdir -p $outDir

# Apply the distribution from the training set to the perturbation data
Rscript norm_perts_data.R $input --clin_dat $clin --train_test $train_test

norm_input=$(echo "$input" | sed 's/.csv$//')
norm_input="${norm_input}_clinRef_qnorm.csv"

# Do PCA before the SVA
Rscript bigPCA_exe.R $norm_input --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$norm_input" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Do SVA
Rscript create_batchRemObjkts.R $norm_input \
    --metDat $metDat \
    --whatData "pert" \
    --subSampSVAmods $subsamp \
    --outDir $outDir

mod=$(echo "$norm_input" | sed 's/.csv$//')
mod="${mod}_svaMod.rds"
mod=$(basename $mod)
mod0=$(echo "$norm_input" | sed 's/.csv$//')
mod0="${mod0}_svaMod0.rds"
mod0=$(basename $mod0)

Rscript sva_fast_exe.R $norm_input \
    --mod "${outDir}${mod}" \
    --mod0 "${outDir}${mod0}" \
    --nSV_method $nSV_method \
    --saveSVrem \
    --outDir $outDir

svaAdj=$(basename "$norm_input")
svaAdj=$(echo "$svaAdj" | sed 's/.csv$//')
svaAdj="${outDir}${svaAdj}_svaMod_svaAdj.csv"

# Do PCA after the SVA
Rscript bigPCA_exe.R $svaAdj --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$svaAdj" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Apply the distribution from the training set again, to adjust
# whatever changes SVA might have done on it
Rscript norm_perts_data.R $svaAdj --clin_dat $clin --train_test $train_test

svaAdj_clinNorm=$(basename "$svaAdj")
svaAdj_clinNorm=$(echo "$svaAdj_clinNorm" | sed 's/.csv$//')
svaAdj_clinNorm="${outDir}${svaAdj_clinNorm}_clinRef_qnorm.csv"

# Do PCA after the qnorm with training reference
Rscript bigPCA_exe.R $svaAdj_clinNorm --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$svaAdj_clinNorm" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir