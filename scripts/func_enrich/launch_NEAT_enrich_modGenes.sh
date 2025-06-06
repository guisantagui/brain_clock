#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=NEAT_enr_modGns
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output=output_NEAT_enrich_modGenes.txt
#Define sderr path
#SBATCH --error=error_NEAT_enrich_modGenes.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# Run NEAT enrichment using the genes with non-zero coefficients in the GLM obtained
# in the second fitting round.


conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
DE_file="../../results/models/secnd_round/mod_alpha1_coefs.csv"
outNameNEAT="modGenes_neat.csv"
extGeneSets=NULL
FCVersion=5
conv2ENSEMBL=false
FCFile="~/projects/cureMILS/proteomics/data/FunCoup/FC5.0_H.sapiens_full"
outDir="../../results/enrich_NEAT/mod_second_round/"

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