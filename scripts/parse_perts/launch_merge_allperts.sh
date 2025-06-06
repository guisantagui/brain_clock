#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=merge_allperts
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 8
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/output_merge_allperts.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/error_merge_allperts.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

Rscript merge_lincs_pert111.R

conda deactivate