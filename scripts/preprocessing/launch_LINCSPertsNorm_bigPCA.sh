#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=lincsPrtsPCA
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_lincsPrtsPCA.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_lincsPrtsPCA.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

dat="../../results/parsed/merged_perts/merged_perts_counts_quantNorm.csv"
metDat="../../results/parsed/merged_perts/merged_perts_metdat.csv"
outDir="../../results/preproc/lincs_pert_norm_pca/"
nPCs=20

Rscript bigPCA_exe.R $dat --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$dat" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

conda deactivate