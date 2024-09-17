#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=pcaPredSVs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_pcaPredSVs.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_pcaPredSVs.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=long

conda activate r-4.3.1

# This script runs PCA of dataset built by merging SVAed dataset constituting everything but LINCS, and LINCS
# with batch effect correction based on predicted SVs, obtained with SV predictive model.

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integRealSVA_LINCSpredSVA/integRealSVA_LINCSpredSVA.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"
filtDF="none" # "none" for not filtering
#propZerosRem=0.8
nPCs=20
#tiss2rem="cerebellum,cerebellar hemisphere"
#outTag="noCerebell"
excludeSubstudy="none"
#nSV_method="leek"
#rinFilt=6
outDir="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integRealSVA_LINCSpredSVA/"

# Create output directory if it doesn't exist
if [ ! -d "$outDir" ]; then
    mkdir -p "$outDir"
fi

# Filter input dataframe, if specified
if [[ "$filtDF" != "none" || "$excludeSubstudy" != "none" ]]; then
    filtOutName=$(dirname $input)
    filtOutName="$filtOutName/"
    Rscript filt_genes.R $input --filtFile $filtDF --metDat $metDat --excludeSubstudy $excludeSubstudy --outDir $filtOutName
    filtDF_bn=$(basename "$filtDF")
    filtInput=$(echo "$input" | sed 's/.csv$//')
    filtInput="${filtInput}_${filtDF_bn}"
else
    filtInput="$input"
fi

Rscript bigPCA_exe.R $input --nPCs $nPCs --stand --outDir $outDir
pcaInFile=$(echo "$input" | sed 's/.csv$//')
pcaInFile="${pcaInFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaInFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

conda deactivate