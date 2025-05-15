#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=preproc_first
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/output_preproc_first.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/preprocessing/error_preproc_first.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# This version of the preproc pipeline uses as input the quantNorm(log(counts + 1)) of the integrated dataset,
# consisting of synapse.org datasets, TBI, 111 compilation of GEO obtained by Sascha and Javier for brain specific
# cell types and NPC, NEU and MIC samples of LINCS1000 (level 3), and SC data from ageAnno. This dataset has been
# already preprocessed in terms of removal of samples with low RIN, variables, samples with high proportion of zeros,
# etc. 

# Variables for the pipeline
########################################################################################################################
input="../../results/parsed/merged/merged_counts.csv"
metDat="../../results/parsed/merged/merged_metdat.csv"
filtDF="none" # "none" for not filtering
nPCs=20
trainProp=0.66
tiss2rem="cerebellum,cerebellar hemisphere"
outTag="noCerebell"
excludeSubstudy="none" # "none" for not excluding any substudy before the preprocessing
nSV_method="leek"
outDir="../../results/preproc/first_round/"

# Create output directory if it doesn't exist
if [ ! -d "$outDir" ]; then
    mkdir -p "$outDir"
fi

# Filter input dataframe, if specified
if [ "$filtDF" != "none" ]; then
    filtOutName=$(dirname $input)
    filtOutName="$filtOutName/"
    Rscript filt_genes.R $input --filtFile $filtDF --metDat $metDat --excludeSubstudy $excludeSubstudy --outDir $filtOutName
    filtDF_bn=$(basename "$filtDF")
    filtInput=$(echo "$input" | sed 's/.csv$//')
    filtInput="${filtInput}_${filtDF_bn}"
else
    filtInput="$input"
fi

# Normalize the dataset using as reference the entire input dataset: this normalization is only
# for visualizing the data in PCAs.
Rscript normalize_dataset.R $filtInput
norm_input="${filtInput/.csv/_log2_qnorm.csv}"

# Do PCA of preprocessed file and plot it.
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

# Do train-test splitting here, as we will use it for normalization to prevent
# data leakage, and we have already removed cerebellum samples.
train_test_outdir=$(echo $(dirname $input))
Rscript train_test_split.R $noCerebFile \
    --metDat $metDat \
    --trainProp $trainProp \
    --exclude_substudy "brainSeq_pI" \
    --outdir $train_test_outdir

# Normalize raw data again, but now using only train set as reference. Will
# overwrite normalized file in $train_test_outdir
Rscript normalize_dataset.R $filtInput  \
    --train_test "${train_test_outdir}/train_test.csv"

# Remove cerebellum samples again. Will overwrite the noCerebell file previously
# generated in $outDir
Rscript tpms_tissOutFilt.R $norm_input \
    --metDat $metDat \
    --tiss2rem "$tiss2rem" \
    --outTag $outTag \
    --outDir $outDir

# Create objects necessary for batch removal steps
Rscript create_batchRemObjkts.R $noCerebFile --metDat $metDat --whatData "clin" --outDir $outDir

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

# Remove batch effect: ComBat --> SVA

# These lines are for running ComBat
Rscript combat_exe.R $noCerebFile --batch $batch --combatMod $modCombat --outDir $outDir
combatAdj=$(echo "$noCerebFile" | sed 's/.csv$//')
combatAdj="${combatAdj}_combat.csv"

# These lines are for running only SVA
Rscript sva_fast_exe.R $noCerebFile --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --saveSVrem --outDir $outDir
svaAdj_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
svaAdj_onlyAge="${svaAdj_onlyAge}_onlyAge_svaAdj.csv"

# Run SVA on ComBat-ed data
Rscript sva_fast_exe.R $combatAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --saveSVrem --outDir $outDir
cmbt_svaAdj_onlyAge=$(echo "$combatAdj" | sed 's/.csv$//')
cmbt_svaAdj_onlyAge="${cmbt_svaAdj_onlyAge}_onlyAge_svaAdj.csv"

# Run PCAs of the results
Rscript bigPCA_exe.R $combatAdj --nPCs $nPCs --stand --outDir $outDir
Rscript bigPCA_exe.R $svaAdj_onlyAge --nPCs $nPCs --stand --outDir $outDir
Rscript bigPCA_exe.R $cmbt_svaAdj_onlyAge --nPCs $nPCs --stand --outDir $outDir

# Run SVA on the SVA-removed dataset and plot surrogate
# variables to see if batch effect is still present (without saving df
# with regressed-out SVs)
Rscript sva_fast_exe.R $svaAdj_onlyAge --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --outDir $outDir
sva_after_sva_onlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
sva_after_sva_onlyAge="${sva_after_sva_onlyAge}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj_onlyAge --sva $sva_after_sva_onlyAge --metDat $metDat --outDir $outDir

Rscript sva_fast_exe.R $combatAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --outDir $outDir
sva_after_combat=$(echo "$combatAdj" | sed 's/.csv$//')
sva_after_combat="${sva_after_combat}_onlyAge_svobj.rds"
Rscript sva_plot.R $combatAdj --sva $sva_after_combat --metDat $metDat --outDir $outDir

Rscript sva_fast_exe.R $cmbt_svaAdj_onlyAge --mod $mod_onlyAge --mod0 $mod0_onlyAge --nSV_method $nSV_method --outDir $outDir
sva_after_combat_and_sva=$(echo "$cmbt_svaAdj_onlyAge" | sed 's/.csv$//')
sva_after_combat_and_sva="${sva_after_combat_and_sva}_onlyAge_svobj.rds"
Rscript sva_plot.R $cmbt_svaAdj_onlyAge --sva $sva_after_combat_and_sva --metDat $metDat --outDir $outDir

# Plot PCAs
pcaSVAOnlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
pcaSVAOnlyAge="${pcaSVAOnlyAge}_pca.rds"
Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir

pcaComBat=$(echo "$combatAdj" | sed 's/.csv$//')
pcaComBat="${pcaComBat}_pca.rds"
Rscript plotBigPCA_exe.R $pcaComBat --metDat $metDat --x PC1 --y PC2 --outDir $outDir

pcaSVAafterComBat=$(echo "$cmbt_svaAdj_onlyAge" | sed 's/.csv$//')
pcaSVAafterComBat="${pcaSVAafterComBat}_pca.rds"
Rscript plotBigPCA_exe.R $pcaSVAafterComBat --metDat $metDat --x PC1 --y PC2 --outDir $outDir

conda deactivate