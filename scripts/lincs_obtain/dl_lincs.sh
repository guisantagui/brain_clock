#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=dl_lincs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 8
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_dl_lincs.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_dl_lincs.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

wget -P /mnt/lscratch/users/gsantamaria/test_large_files/ https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/cp_predicted_RNAseq_profiles.gctx &
wget -P /mnt/lscratch/users/gsantamaria/test_large_files/ https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/ctl_predicted_RNAseq_profiles.gctx &
wait



# All this added to do it all at once
conda activate r-4.3.1
# Variables for the pipeline
########################################################################################################################
gctxFile_exp="/mnt/lscratch/users/gsantamaria/test_large_files/cp_predicted_RNAseq_profiles.gctx"
gctxFile_ctl="/mnt/lscratch/users/gsantamaria/test_large_files/ctl_predicted_RNAseq_profiles.gctx"
cellType="NPC"
batchSize=8000
gctxDir="/mnt/lscratch/users/gsantamaria/test_large_files/NPC/gctx/"

# Run
########################################################################################################################
Rscript lincs_filt.R $gctxFile_exp --cellType $cellType --batchSize $batchSize --outDir $gctxDir &
Rscript lincs_filt.R $gctxFile_ctl --cellType $cellType --batchSize $batchSize --outDir $gctxDir &

wait

# Parsing part
########################################################################################################################
outDir="/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/"
files_with_cType=$(find $gctxDir -type f -exec grep -l $cellType {} +)

# Optionally, you can iterate over each file
for file in $files_with_cType; do
    echo "Processing file: $file" &
    Rscript lincs_parse.R $file --outDir $outDir &
done

wait

conda deactivate