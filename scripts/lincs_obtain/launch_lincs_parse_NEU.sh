#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=lincsParseNEU
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 8
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/output_lincs_parse_NEU.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/error_lincs_parse_NEU.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
gctxDir="/mnt/lscratch/users/gsantamaria/test_large_files/NEU/gctx/"
cellType="NEU"



#gctxFile="/mnt/lscratch/users/gsantamaria/test_large_files/ctl_predicted_RNAseq_profiles.gctx"
outDir="/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/"

# Run
########################################################################################################################
files_with_cType=$(find $gctxDir -type f -exec grep -l $cellType {} +)

# Optionally, you can iterate over each file
for file in $files_with_cType; do
    echo "Processing file: $file" &
    Rscript lincs_parse.R $file --outDir $outDir &
done

wait

