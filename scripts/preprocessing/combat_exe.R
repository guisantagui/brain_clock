################################################################################
# Script for running Combat to remove the influcence from known batches.       #
################################################################################

if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require(sva, quietly = T)) BiocManager::install("sva", update = F)
library(sva)

# Terminal argument parser
################################################################################
parser <- arg_parser("Run ComBat.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--batch",
                               "--combatMod",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "RDS object of a vector of the batch of each sample in the input matrix.",
                                "ComBat model object.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)

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
dat <- parsed$input
batch <- parsed$batch
combatMod <- parsed$combatMod
outDir <- parsed$outDir

outBaseName <- gsub(".csv", "", basename(parsed$input))

# Load data
################################################################################
batch <- readRDS(batch)
modcombat <- readRDS(combatMod)
dat <- readCsvFast(dat)

# Run ComBat
################################################################################
edata <- t(dat)

print("Running ComBat...")
combat_edata <- ComBat(dat = edata,
                       batch = batch,
                       mod = modcombat,
                       par.prior = TRUE,
                       prior.plots = FALSE)

# Return to the format of the input (vars in columns, samples in rows)
combat_edata <- t(combat_edata)

outName <- sprintf("%s%s_combat.csv", outDir, outBaseName)
writeCsvFst(combat_edata, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))