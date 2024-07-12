#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=svaFastAll
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 8
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_svaFastAll_lincsLndmrk.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_svaFastAll_lincsLndmrk.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# This version of the preproc pipeline uses as input the quantNorm(log(counts + 1)) of the integrated dataset,
# consisting of synapse.org datasets, TBI, 111 compilation of GEO obtained by Sascha and Javier for brain specific
# cell types and NPC samples of LINCS1000 (level 3). This dataset has been already preprocessed in terms of removal
# of samples with low RIN, variables, samples with high proportion of zeros, etc. 

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_counts_wTBI_wPert111_log2_quantNorm_preproc_wLINCS.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wLINCS.csv"
filtDF="none" # "none" for not filtering
#propZerosRem=0.8
nPCs=20
tiss2rem="cerebellum,cerebellar hemisphere"
outTag="noCerebell"
#rinFilt=6
outDir="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_all_sva_fast/"

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

# Preprocess the TPMs file
#Rscript tpms_preproc.R $input --metDat $metDat --propZerosRem $propZerosRem --rinFilt $rinFilt --outDir $outDir

#logSuff="_preproc_log10.csv"

#inPCALog=$(basename "$input")
#inPCALog=$(echo "$inPCALog" | sed 's/.csv$//')
#inPCALog="${outDir}${inPCALog}${logSuff}"

# Do PCA of preprocessed file and plot it.
#Rscript bigPCA_exe.R $filtInput --nPCs $nPCs --stand --outDir $outDir
#pcaFile=$(echo "$filtInput" | sed 's/.csv$//')
#pcaFile=$(basename $pcaFile)
#pcaFile="${outDir}${pcaFile}_pca.rds"
#Rscript plotBigPCA_exe.R $pcaFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Remove cerebellum samples and run PCA again
Rscript tpms_tissOutFilt.R $filtInput --metDat $metDat --tiss2rem "$tiss2rem" --outTag $outTag --outDir $outDir

noCerebFile=$(echo "$filtInput" | sed 's/.csv$//')
noCerebFile=$(basename $noCerebFile)
noCerebFile="${outDir}${noCerebFile}_${outTag}.csv"

#Rscript bigPCA_exe.R $noCerebFile --nPCs $nPCs --stand --outDir $outDir
#pcaNoCerebFile=$(echo "$noCerebFile" | sed 's/.csv$//')
#pcaNoCerebFile="${pcaNoCerebFile}_pca.rds"
#Rscript plotBigPCA_exe.R $pcaNoCerebFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Create objects necessary for batch removal steps
Rscript create_batchRemObjkts.R $noCerebFile --metDat $metDat --outDir $outDir

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

Rscript sva_fast_exe.R $noCerebFile --mod $mod_onlyAge --mod0 $mod0_onlyAge --saveSVrem --outDir $outDir &
#Rscript sva_exe.R $noCerebFile --mod $mod_all --mod0 $mod0_all --saveSVrem --outDir $outDir &
#Rscript combat_exe.R $noCerebFile --batch $batch --combatMod $modCombat --outDir $outDir &
wait

svaAdj_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
svaAdj_onlyAge="${svaAdj_onlyAge}_onlyAge_svaAdj.csv"

#svaAdj_all=$(echo "$noCerebFile" | sed 's/.csv$//')
#svaAdj_all="${svaAdj_all}_all_svaAdj.csv"

#combatAdj=$(echo "$noCerebFile" | sed 's/.csv$//')
#combatAdj="${combatAdj}_combat.csv"

# Run PCAs of the results
Rscript bigPCA_exe.R $svaAdj_onlyAge --nPCs $nPCs --stand --outDir $outDir &
#PID_sva_onlyAge=$!

#Rscript bigPCA_exe.R $svaAdj_all --nPCs $nPCs --stand --outDir $outDir &
#PID_sva_all=$!

#Rscript bigPCA_exe.R $combatAdj --nPCs $nPCs --stand --outDir $outDir &
#PID_combat=$!
wait

# Run SVA of the three batch effect removal approaches
# and plot surrogate variables to see if batch effect is 
# still present (without saving df with regressed-out SVs)
Rscript sva_fast_exe.R $svaAdj_onlyAge --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
#Rscript sva_exe.R $svaAdj_all --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
#Rscript sva_exe.R $combatAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
wait

sva_after_sva_onlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
sva_after_sva_onlyAge="${sva_after_sva_onlyAge}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj_onlyAge --sva $sva_after_sva_onlyAge --metDat $metDat --outDir $outDir

#sva_after_sva_all=$(echo "$svaAdj_all" | sed 's/.csv$//')
#sva_after_sva_all="${sva_after_sva_all}_onlyAge_svobj.rds"
#Rscript sva_plot.R $svaAdj_all --sva $sva_after_sva_all --metDat $metDat --outDir $outDir

#sva_after_combat=$(echo "$combatAdj" | sed 's/.csv$//')
#sva_after_combat="${sva_after_combat}_onlyAge_svobj.rds"
#Rscript sva_plot.R $combatAdj --sva $sva_after_combat --metDat $metDat --outDir $outDir

# Plot PCAs
pcaSVAOnlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
pcaSVAOnlyAge="${pcaSVAOnlyAge}_pca.rds"

#pcaSVAAll=$(echo "$svaAdj_all" | sed 's/.csv$//')
#pcaSVAAll="${pcaSVAAll}_pca.rds"

#pcaCombat=$(echo "$combatAdj" | sed 's/.csv$//')
#pcaCombat="${pcaCombat}_pca.rds"

#(wait $PID_sva_onlyAge && Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
#(wait $PID_sva_all && Rscript plotBigPCA_exe.R $pcaSVAAll --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
#(wait $PID_combat && Rscript plotBigPCA_exe.R $pcaCombat --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir
#Rscript plotBigPCA_exe.R $pcaSVAAll --metDat $metDat --x PC1 --y PC2 --outDir $outDir
#Rscript plotBigPCA_exe.R $pcaCombat --metDat $metDat --x PC1 --y PC2 --outDir $outDir

conda deactivate