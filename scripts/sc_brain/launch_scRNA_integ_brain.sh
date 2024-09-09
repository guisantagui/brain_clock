#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=scRNAintBrain
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 1
#SBATCH --time=00-12:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/sc_brain/output_scRNAseq_integ_brain.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/sc_brain/error_scRNAseq_integ_brain.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

module load devel/CMake/3.18.4-GCCcore-10.2.0
module load math/GMP/6.2.0-GCCcore-10.2.0
module load vis/cairo/1.16.0-GCCcore-10.2.0

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

preProcDir="/home/users/gsantamaria/projects/aging/results/preprocessing/brain"
cellCycleGenes="/home/users/gsantamaria/projects/brain_clock/data/scRNAseq/CellCycleGenes_Human.csv"
method="rpca"

Rscript integration_exec.R --preProcDir $preProcDir --method $method


conda deactivate