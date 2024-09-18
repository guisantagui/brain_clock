#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTrainAllLINCS
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/output_model_training_allLINCS_chronAge.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/error_model_training_allLINCS_chronAge.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# This version prefilters the genes before the fitting to keep the ones that were significant in the SC model
# (obtained with mod_train_sc.R) and the ones that were significant with the bulk data after SVAing the whole matrix
# with all LINCS brain cell types, etc.

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
data="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_allLINCSBrain/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC_noCerebell_onlyAge_svaAdj.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"
respVar="age_chron"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
alpha="0"
mem="24G"
preFiltGenes="/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/fullMod_andSCMod_genes.csv" # Combined genes with non-zero coefficients of SC model and full model
braakThrshld="4"
outDir="/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_chron_age_scSignAndFullModSignPreFilt/"

# Train the model
########################################################################################################################
Rscript mod_train_and_test.R $data --metDat $metDat --respVar $respVar --ageTransPars $ageTransPars --alpha $alpha --mem $mem --preFiltGenes $preFiltGenes --braakThrshld $braakThrshld --outDir $outDir

modFile="${outDir}modFuncsAlpha${alpha}.rds"

Rscript model_getCoefs.R $modFile --outDir $outDir

conda deactivate