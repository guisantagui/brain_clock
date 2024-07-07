#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=preproc
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_preproc_pipe.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_preproc_pipe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/combined_TPMs_wTBI_lincs.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/combined_metDat_wTBI.csv"
propZerosRem=0.8
nPCs=20
tiss2rem="cerebellum,cerebellar hemisphere"
outTag="noCerebell"
rinFilt=6
outDir="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/wTPI_lincs/"

# Preprocess the TPMs file
Rscript tpms_preproc.R $input --metDat $metDat --propZerosRem $propZerosRem --rinFilt $rinFilt --outDir $outDir

logSuff="_preproc_log10.csv"

inPCALog=$(basename "$input")
inPCALog=$(echo "$inPCALog" | sed 's/.csv$//')
inPCALog="${outDir}${inPCALog}${logSuff}"

# Do PCA of preprocessed file and plot it.
Rscript bigPCA_exe.R $inPCALog --nPCs $nPCs --stand --outDir $outDir
pcaLogFile=$(echo "$inPCALog" | sed 's/.csv$//')
pcaLogFile="${pcaLogFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaLogFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Remove cerebellum samples and run PCA again
Rscript tpms_tissOutFilt.R $inPCALog --metDat $metDat --tiss2rem "$tiss2rem" --outTag $outTag --outDir $outDir

noCerebFile=$(echo "$inPCALog" | sed 's/.csv$//')
noCerebFile="${noCerebFile}_${outTag}.csv"

Rscript bigPCA_exe.R $noCerebFile --nPCs $nPCs --stand --outDir $outDir
pcaNoCerebFile=$(echo "$noCerebFile" | sed 's/.csv$//')
pcaNoCerebFile="${pcaNoCerebFile}_pca.rds"
Rscript plotBigPCA_exe.R $pcaNoCerebFile --metDat $metDat --x PC1 --y PC2 --outDir $outDir

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

Rscript sva_exe.R $noCerebFile --mod $mod_onlyAge --mod0 $mod0_onlyAge --saveSVrem --outDir $outDir &
Rscript sva_exe.R $noCerebFile --mod $mod_all --mod0 $mod0_all --saveSVrem --outDir $outDir &
Rscript combat_exe.R $noCerebFile --batch $batch --combatMod $modCombat --outDir $outDir &
wait

svaAdj_onlyAge=$(echo "$noCerebFile" | sed 's/.csv$//')
svaAdj_onlyAge="${svaAdj_onlyAge}_onlyAge_svaAdj.csv"

svaAdj_all=$(echo "$noCerebFile" | sed 's/.csv$//')
svaAdj_all="${svaAdj_all}_all_svaAdj.csv"

combatAdj=$(echo "$noCerebFile" | sed 's/.csv$//')
combatAdj="${combatAdj}_combat.csv"

# Run PCAs of the results
Rscript bigPCA_exe.R $svaAdj_onlyAge --nPCs $nPCs --stand --outDir $outDir &
#PID_sva_onlyAge=$!

Rscript bigPCA_exe.R $svaAdj_all --nPCs $nPCs --stand --outDir $outDir &
#PID_sva_all=$!

Rscript bigPCA_exe.R $combatAdj --nPCs $nPCs --stand --outDir $outDir &
#PID_combat=$!
wait

# Run SVA of the three batch effect removal approaches
# and plot surrogate variables to see if batch effect is 
# still present (without saving df with regressed-out SVs)
Rscript sva_exe.R $svaAdj_onlyAge --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
Rscript sva_exe.R $svaAdj_all --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
Rscript sva_exe.R $combatAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir &
wait

sva_after_sva_onlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
sva_after_sva_onlyAge="${sva_after_sva_onlyAge}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj_onlyAge --sva $sva_after_sva_onlyAge --metDat $metDat --outDir $outDir

sva_after_sva_all=$(echo "$svaAdj_all" | sed 's/.csv$//')
sva_after_sva_all="${sva_after_sva_all}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj_all --sva $sva_after_sva_all --metDat $metDat --outDir $outDir

sva_after_combat=$(echo "$combatAdj" | sed 's/.csv$//')
sva_after_combat="${sva_after_combat}_onlyAge_svobj.rds"
Rscript sva_plot.R $combatAdj --sva $sva_after_combat --metDat $metDat --outDir $outDir

# Plot PCAs
pcaSVAOnlyAge=$(echo "$svaAdj_onlyAge" | sed 's/.csv$//')
pcaSVAOnlyAge="${pcaSVAOnlyAge}_pca.rds"

pcaSVAAll=$(echo "$svaAdj_all" | sed 's/.csv$//')
pcaSVAAll="${pcaSVAAll}_pca.rds"

pcaCombat=$(echo "$combatAdj" | sed 's/.csv$//')
pcaCombat="${pcaCombat}_pca.rds"

#(wait $PID_sva_onlyAge && Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
#(wait $PID_sva_all && Rscript plotBigPCA_exe.R $pcaSVAAll --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
#(wait $PID_combat && Rscript plotBigPCA_exe.R $pcaCombat --metDat $metDat --x PC1 --y PC2 --outDir $outDir)
Rscript plotBigPCA_exe.R $pcaSVAOnlyAge --metDat $metDat --x PC1 --y PC2 --outDir $outDir
Rscript plotBigPCA_exe.R $pcaSVAAll --metDat $metDat --x PC1 --y PC2 --outDir $outDir
Rscript plotBigPCA_exe.R $pcaCombat --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Remove effect of age and PMI from the ComBat dataset, as it is the only one that is properly
# removing effect from substudy
Rscript remAgePMI_exe.R $combatAdj --metDat $metDat --outDir $outDir

# Run PCA of the result
combatAgePMIAdj=$(echo "$combatAdj" | sed 's/.csv$//')
combatAgePMIAdj="${combatAgePMIAdj}_pmiAgeAdj.csv"
Rscript bigPCA_exe.R $combatAgePMIAdj --nPCs $nPCs --stand --outDir $outDir

pcaCombatAgePMIAdj=$(echo "$combatAgePMIAdj" | sed 's/.csv$//')
pcaCombatAgePMIAdj="${pcaCombatAgePMIAdj}_pca.rds"
Rscript plotBigPCA_exe.R $pcaCombatAgePMIAdj --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Run SVA again (only giving age info) to assess if there are any remaining batch effects
Rscript sva_exe.R $combatAgePMIAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --saveSVrem --outDir $outDir
svaFile=$(echo "$combatAgePMIAdj" | sed 's/.csv$//')
svaFile="${svaFile}_onlyAge_svobj.rds"
Rscript sva_plot.R $combatAgePMIAdj --sva $svaFile --metDat $metDat --outDir $outDir

# Run a PCA of the resulting dataset with the SVA regressed out.
svaAdj=$(echo "$combatAgePMIAdj" | sed 's/.csv$//')
svaAdj="${svaAdj}_onlyAge_svaAdj.csv"

Rscript bigPCA_exe.R $svaAdj --nPCs $nPCs --stand --outDir $outDir

pcaSvaAdj=$(echo "$svaAdj" | sed 's/.csv$//')
pcaSvaAdj="${pcaSvaAdj}_pca.rds"

Rscript plotBigPCA_exe.R $pcaSvaAdj --metDat $metDat --x PC1 --y PC2 --outDir $outDir

# Run SVA again on this dataset (without saving the dataset with SVs regressed out)
# to assess if the batch effect has been completely removed.
Rscript sva_exe.R $svaAdj --mod $mod_onlyAge --mod0 $mod0_onlyAge --outDir $outDir
svaFileFinal=$(echo "$svaAdj" | sed 's/.csv$//')
svaFileFinal="${svaFileFinal}_onlyAge_svobj.rds"
Rscript sva_plot.R $svaAdj --sva $svaFileFinal --metDat $metDat --outDir $outDir

conda deactivate