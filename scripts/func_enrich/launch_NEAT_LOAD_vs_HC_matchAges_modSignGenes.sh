#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=NEAT_yooMAges
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/func_enrich/output_NEAT_yooMAges.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/func_enrich/error_NEAT_yooMAges.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# Run NEAT enrichment using the genes with non-zero coefficients in  the GLM obtained by fitting 
# integrated dataset + LINCS NPCs + SC data with chronological age as response variable.


conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
DE_file="/home/users/gsantamaria/projects/brain_clock/results/yoo_science_results/LOAD_vs_HC_matchAges_predSign.csv"
outNameNEAT="LOAD_vs_HC_matchAges_modSignGenes.csv"
extGeneSets=NULL
FCVersion=5
conv2ENSEMBL=false
FCFile="/home/users/gsantamaria/projects/cureMILS/proteomics/data/FunCoup/FC5.0_H.sapiens_full"
outDir="/home/users/gsantamaria/projects/brain_clock/results/enrich_NEAT/yoo_science/"

# If the datasets and dataformat are not installed this script will install them
bash install_ncbi_CLTools.sh

# Create the dataframe with ENSEMBL IDs
if [ $conv2ENSEMBL = true ]
then
    Rscript symbs_2_IDs.R $DE_file --whatOutID ENSEMBL --remResDups
    csvExt=".csv"
    DE4NEAT="${DE_file/$csvExt/}"
    DE4NEAT=$DE4NEAT"_ENSEMBL_IDs.csv"
else
    DE4NEAT=$DE_file
fi

# Run NEAT
Rscript NEAT_enrich.R --DE $DE4NEAT \
                      --outName $outNameNEAT \
                      --extGeneSets $extGeneSets \
                      --FCFile $FCFile \
                      --FCVersion $FCVersion \
                      --outDir $outDir

conda deactivate