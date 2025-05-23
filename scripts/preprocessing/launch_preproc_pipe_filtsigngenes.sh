#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=preproc_secnd
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=07-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_preproc_secnd.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_preproc_secnd.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=long

conda activate r-4.3.1

# This version of the preproc pipeline uses as input the clinical dataset counts,
# consisting of synapse.org datasets, TBI, brainSeq Phase 1 and 2 and SC data from ageAnno. This dataset has been already
# filtered in terms of removal of samples with low RIN, variables, samples with high proportion of zeros, etc. 
# Previous to doing sva filters the genes to keep only the ones that had non-zero coefficients in the fitting of
# the GLM with all the genes using chronological age as response variable.

# Variables for the pipeline
########################################################################################################################

input="../../results/parsed/merged/merged_counts.csv"
metDat="../../results/parsed/merged/merged_metdat.csv"
filtDF="../../results/models/first_round_sva_strat_wBSPIext/mod_alpha1_coefs.csv" # "none" for not filtering
nPCs=20
train_test="../../brain_clock/results/parsed/merged/train_test.csv"
tiss2rem="cerebellum,cerebellar hemisphere"
outTag="noCerebell"
excludeSubstudy="none" # "none" for not excluding any substudy before the preprocessing
nSV_method="leek"
outDir="../../results/preproc/second_round/"

# Create output directory if it doesn't exist
if [ ! -d "$outDir" ]; then
    mkdir -p "$outDir"
fi

# Filter input dataframe, if specified
if [ "$filtDF" != "none" ]; then
    filtOutName=$(dirname $input)
    filtOutName="$filtOutName/"
    Rscript filt_genes.R $input --filtFile $filtDF --outDir $filtOutName
    filtDF_bn=$(basename "$filtDF")
    filtInput=$(echo "$input" | sed 's/.csv$//')
    filtInput="${filtInput}_${filtDF_bn}"
else
    filtInput="$input"
fi

# Normalize the dataset using as reference the training samples defined
# during the first fitting round
Rscript normalize_dataset.R $filtInput --train_test $train_test
norm_input="${filtInput/.csv/_log2_qnorm.csv}"

# Do PCA of data.csv file and plot it.
Rscript bigPCA_exe.R $norm_input --nPCs $nPCs --stand --outDir $outDir
pcaFile=$(echo "$norm_input" | sed 's/.csv$//')
pcaFile=$(basename $pcaFile)
pcaFile="${outDir}${pcaFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Remove cerebellum samples and run PCA again
Rscript tpms_tissOutFilt.R $norm_input --metDat $metDat --tiss2rem "$tiss2rem" --outTag $outTag --outDir $outDir

noCerebFile=$(echo "$norm_input" | sed 's/.csv$//')
noCerebFile=$(basename $noCerebFile)
noCerebFile="${outDir}${noCerebFile}_${outTag}.csv"

Rscript bigPCA_exe.R $noCerebFile --nPCs $nPCs --stand --outDir $outDir
pcaNoCerebFile=$(echo "$noCerebFile" | sed 's/.csv$//')
pcaNoCerebFile="${pcaNoCerebFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaNoCerebFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Create objects necessary for batch removal steps
Rscript create_batchRemObjkts.R $noCerebFile --metDat $metDat  --whatData "clin" --outDir $outDir

# Run SVA with no covariates and with covariates, and ComBat
mod_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
mod_onlyAge="${mod_onlyAge}_svaMod_onlyAge.rds"
mod0_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
mod0_onlyAge="${mod0_onlyAge}_svaMod0_onlyAge.rds"

mod_all=$(echo "$noCerebFile" | sed 's/.csv$//')
mod_all="${mod_all}_svaMod_all.rds"
mod0_all=$(echo "$noCerebFile" | sed 's/.csv$//')
mod0_all="${mod0_all}_svaMod0_all.rds"

batch=$(echo "$noCerebFile" | sed 's/.csv$//')
batch="${batch}_batches.rds"

modCombat=$(echo "$noCerebFile" | sed 's/.csv$//')
modCombat="${modCombat}_combatMod.rds"

Rscript sva_fast_exe.R $noCerebFile --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --saveSVrem --outDir $outDir

svaAdj_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
svaAdj_onlyAge="${svaAdj_onlyAge}_onlyAge_svaAdj.csv"

# Run PCAs of the results
Rscript bigPCA_exe.R $svaAdj_onlyAge --nPCs $nPCs --stand --outDir $outDir

# Run SVA on the SVA-removed dataset and plot surrogate
# variables to see if batch effect is still present (without saving df
# with regressed-out SVs)
Rscript sva_fast_exe.R $svaAdj_onlyAge --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --outDir $outDir

sva_after_sva_onlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
sva_after_sva_onlyAge="${sva_after_sva_onlyAge}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj_onlyAge --sva $sva_after_sva_onlyAge --metDat $metDat --outDir $outDir

# Plot PCAs
pcaSVAOnlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
pcaSVAOnlyAge="${pcaSVAOnlyAge}_pca.rds"

Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir

conda deactivate