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
if(!require("argparser")){
        install.packages("argparser",
                         repos = 'https://pbil.univ-lyon1.fr/CRAN/')
}

library(limma)
library(dplyr)
library(argparser)

# Terminal argument parser
################################################################################

parser <- arg_parser("scRNAseq preprocessing")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--svobj",
                               "--outDir"),
                       help = c("Non sva-treated expression matrix (genes in columns, samples in rows).",
                                "SV RDS object (generated with sva_fast_exe.R or sva_exe.R).",
                                "Output directory"),
                       flag = c(F, F, F))

parsed <- parse_args(parser)


# Functions
################################################################################
readCsvFst <- function(path){
        df <- data.frame(data.table::fread(path))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Directory stuff
################################################################################
inFile <- parsed$input
svObjFile <- parsed$svobj
outDir <- parsed$outDir

outDir <- addSlashIfNot(outDir)
createIfNot(outDir)

outName <- sprintf("%ssva_pred_mod.rds", outDir)

# Load data
################################################################################
print("Loading the data...")
svobj <- readRDS(svObjFile)
dat <- readCsvFst(inFile)

# Fit model
################################################################################
sv_mat <- svobj$sv
colnames(sv_mat) <- paste("SV", 1:ncol(sv_mat), sep = "_")

dat <- dat[match(rownames(sv_mat), rownames(dat)), ]

print("Fitting the model...")
surrogate_model <- lm(sv_mat ~ ., data = as.data.frame(dat))

saveRDS(surrogate_model, file = outName)

print(sprintf("%s saved at %s", basename(outName), dirname(outName)))