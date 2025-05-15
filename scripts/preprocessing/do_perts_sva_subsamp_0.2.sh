#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=perts_preproc_sub0.2
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_perts_preproc_sub0.2.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_perts_preproc_sub0.2.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
input="../../results/parsed/merged_perts/merged_perts_exprsn.csv"
clin="../../results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
train_test="../../results/parsed/merged/train_test.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/merged_perts_metdat.csv"
nSV_method="leek"
subsamp=0.2
nPCs=20
outDir="../../results/preproc/perts_preproc_subsamp0.2/"

mkdir -p $outDir

# Normalize the merged dataset
Rscript norm_perts_data.R $input --metDat $metDat
norm_input="${input/.csv/_log2111_qnorm.csv}"

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

# Apply the distribution from the controls to the perturbation data
Rscript norm_perts_data.R $svaAdj --clin_dat $clin --train_test $train_test