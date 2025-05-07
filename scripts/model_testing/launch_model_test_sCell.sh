#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTstSCell
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output=output_model_test_sCell.txt
#Define sderr path
#SBATCH --error=error_model_test_sCell.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
datFile="../../results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
modDir="../../results/models/secnd_round/mod_alpha1/"
metDat="../../results/parsed/merged/merged_metdat.csv"
modCoefsFile="../../results/models/secnd_round/mod_alpha1_coefs.csv"
batchSize=8000
whatSampsTest="single_cell"
mem="50G"
youngThrshld=30
oldThrshld=70
outDir="../../results/model_test_sCell/"

# Run the simulations
########################################################################################################################

Rscript model_test.R $datFile \
        --modDir $modDir \
        --metDat $metDat \
        --sizeBatch $batchSize \
        --whatSampsTest $whatSampsTest \
        --mem $mem \
        --outDir $outDir

sCell_stats_in="${outDir}pred_ages.csv"

Rscript model_sCell_stats.R $sCell_stats_in \
        --metDat $metDat \
        --excludeYoung

Rscript get_sc_gene_dir_changes_prop.R $datFile \
        --metDat $metDat \
        --modCoefs $modCoefsFile \
        --youngThrshld $youngThrshld \
        --oldThrshld $oldThrshld \
        --outDir $outDir

conda deactivate