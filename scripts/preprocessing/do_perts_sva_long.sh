#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=perts_preproc_lng
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500GB
#SBATCH -c 8
#SBATCH --time=14-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_perts_preproc_lng.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_perts_preproc_lng.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=long

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
input="/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/merged_perts_exprsn_qnorm.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/merged_perts_metdat.csv"
nSV_method="leek"
nPCs=20
outDir="../../results/preproc/perts_preproc_lng/"

# Do PCA before the SVA
Rscript bigPCA_exe.R $input --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$input" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Do SVA
Rscript create_batchRemObjkts.R $input \
    --metDat $metDat \
    --whatData "pert" \
    --outDir $outDir

mod=$(echo "$input" | sed 's/.csv$//')
mod="${mod}_svaMod.rds"
mod=$(basename $mod)
mod0=$(echo "$input" | sed 's/.csv$//')
mod0="${mod0}_svaMod0.rds"
mod0=$(basename $mod0)

Rscript sva_fast_exe.R $input \
    --mod "${outDir}${mod}" \
    --mod0 "${outDir}${mod0}" \
    --nSV_method $nSV_method \
    --saveSVrem \
    --outDir $outDir

# Do PCA after the SVA
Rscript bigPCA_exe.R $input --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$dat" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir