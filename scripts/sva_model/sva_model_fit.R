################################################################################
# Brain clock: build predictive model of surrogate variables from integrated.  #
# database, with the aim of using it with new samples user might give, to      #
# estimate surrogate variables in new data and be able to adjust for new batch #
# effects.                                                                     #
################################################################################

if(!require("BiocManager", quietly = T)){
        install.packages("BiocManager",
        repos = "https://pbil.univ-lyon1.fr/CRAN/")
}

if(!require("limma")) BiocManager::install("limma")

library(limma)
library(dplyr)

svObjFile <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_noCerebell_onlyAge_svobj.rds"
inFile <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_noCerebell.csv"

readCsvFst <- function(path){
        df <- data.frame(data.table::fread(path))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}


# Load data
################################################################################
svobj <- readRDS(svObjFile)
dat <- readCsvFst(inFile)

# Fit model
################################################################################
sv_mat <- svobj$sv
colnames(sv_mat) <- paste("SV", 1:ncol(sv_mat), sep = "_")

dat <- dat[match(rownames(sv_mat), rownames(dat)), ]

surrogate_model <- lm(sv_mat ~ ., data = as.data.frame(dat))