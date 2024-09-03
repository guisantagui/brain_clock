#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTstPrts
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:40:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_test_perts/output_model_test_perts.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_test_perts/error_model_test_perts.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
datFile="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_noCerebell_onlyAge_svaAdj.csv"
modFile="/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWLincs_and_sc_sva_chron_age/modFuncsAlpha1/GLM_model_R_1724154505295_1"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"
respVar="age_chron"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
batchSize=8000
whatSampsTest="single_cell"
mem="50G"
outDir="/home/users/gsantamaria/projects/brain_clock/results/model_test_sCell/modAllGenes_integWLincs_and_sc_oAge_chronAge_alph1/"

# Run the simulations
########################################################################################################################

Rscript model_test.R $datFile --modFile $modFile --metDat $metDat --respVar $respVar --ageTransPars $ageTransPars --sizeBatch $batchSize --whatSampsTest $whatSampsTest --mem $mem --outDir $outDir

conda deactivate