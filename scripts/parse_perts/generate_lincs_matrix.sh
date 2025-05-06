#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=gen_lincs_mat
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 8
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/parse_perts/output_gen_lincs_mat.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/parse_perts/error_gen_lincs_mat.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

# Download LINCS L1000 level 3 GCTX files
########################################################################################################################
out_dir="../../data/perturbation/lincs/"

mkdir -p $out_dir

wget -P "$out_dir" https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/cp_predicted_RNAseq_profiles.gctx &
wget -P "$out_dir" https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/ctl_predicted_RNAseq_profiles.gctx &
wait


# All this added to do it all at once
conda activate r-4.3.1
# Variables for the pipeline
########################################################################################################################
gctxFile_exp="${out_dir}cp_predicted_RNAseq_profiles.gctx"
gctxFile_ctl="${out_dir}ctl_predicted_RNAseq_profiles.gctx"
cellTypes="NPC NEU MICROGLIA"
batchSize=8000

# Generate the files for each cell type
########################################################################################################################
for cell in $cellTypes; do
    # Create smaller GCTX with only the target cell type for both controls and perturbation in batches of $batchSize
    gctxDir="${oud_dir}${cell}/gctx/"
    Rscript lincs_filt.R $gctxFile_exp --cellType $cellType --batchSize $batchSize --outDir $gctxDir &
    Rscript lincs_filt.R $gctxFile_ctl --cellType $cellType --batchSize $batchSize --outDir $gctxDir &
    wait
    # Extract CSV files for metadata and expression from each GCTX file
    out_dir_pars="${out_dir}${cell}/parsed_mats/"
    files_with_cType=$(find $gctxDir -type f -exec grep -l $cell {} +)
    for file in $files_with_cType; do
        echo "Processing file: $file" &
        Rscript lincs_parse.R $file --outDir $outDir &
    done
    wait
    # Concatenate the different CSV files of each cell type into a single file
    Rscript concat_parsed_lincs_files.R $out_dir_pars --cellType $cell --whichDF "expMat" --outDir $out_dir_pars &
    Rscript concat_parsed_lincs_files.R $out_dir_pars --cellType $cell --whichDF "metDat" --outDir $out_dir_pars &
    wait
done


conda deactivate