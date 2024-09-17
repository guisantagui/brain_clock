################################################################################
# Brain clock: predict surrogate variables using model built with              #
# sva_model_fit.R.                                                             #
################################################################################
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
library(limma)
# Functions
################################################################################

# Read csv faster
readCsvFast <- function(f){
        df <- data.frame(data.table::fread(f))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# write.csv, but faster
writeCsvFst <- function(df, file, rowNames = T, colNames = T){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames)
}

# Directory stuff
################################################################################

datFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC.csv"
metDatFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"
modFile <- "/home/users/gsantamaria/projects/brain_clock/results/svaMod/allButLINCS_signChron/sva_pred_mod.rds"
outName <- "/home/users/gsantamaria/projects/brain_clock/results/svaMod/allButLINCS_signChron/LINCS_svaPredCorr.csv"

substudFilt <- "LINCS"

# Load data
################################################################################

df <- readCsvFast(datFile)
if(substudFilt != "none"){
        metDat <- readCsvFast(metDatFile)
        substudSamps <- metDat$specimenID[metDat$substudy == substudFilt]
        df <- df[make.names(rownames(df)) %in% make.names(substudSamps), ]
}

mod <- readRDS(modFile)

genesMod <- rownames(mod$coefficients)
genesMod <- genesMod[!grepl("Intercept", genesMod)]

df <- df[, match(genesMod, colnames(df))]


# Predict SVs
################################################################################
preds <- predict(mod, df)

# Regress out SVs and save adjusted matrix
################################################################################
df_adj <- removeBatchEffect(t(df), covariates = preds)

df_adj <- t(df_adj)

writeCsvFst(df_adj, file = outName)