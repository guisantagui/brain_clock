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