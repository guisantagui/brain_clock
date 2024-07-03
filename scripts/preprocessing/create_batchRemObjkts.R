################################################################################
# Create necessary files for batch effect removal                              #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("BiocManager", quietly = TRUE)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
library(limma)

# Parser
################################################################################
parser <- arg_parser("Create objects needed for running batch effect removal")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata CSV file.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
inFile <- parsed$input
metDatFile <- parsed$metDat
outDir <- parsed$outDir

outBaseName <- gsub(".csv", "", basename(inFile))

# Load data
################################################################################
metDat <- read.csv(metDatFile, row.names = 1)
input <- data.frame(data.table::fread(inFile))
rownames(input) <- input$V1
input <- input[, colnames(input) != "V1"]


substudNew <- paste("ROSMAP",
                    metDat$batch_seq[metDat$substudy == "ROSMAP"],
                    sep = "_")

metDat$substudy[metDat$substudy == "ROSMAP"] <- substudNew

pheno <- data.frame(age_death = metDat$ageDeath[match(rownames(input),
                                                      make.names(metDat$specimenID))],
                    pmi = metDat$pmi[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    sex = metDat$sex[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    substudy = metDat$substudy[match(rownames(input),
                                                     make.names(metDat$specimenID))])

# Create mod objects for running SVA
################################################################################

# Substutute NAs in PMI for the median value
pheno$pmi[is.na(pheno$pmi)] <- median(pheno$pmi[!is.na(pheno$pmi)])

# Create design matrixes for SVA accounting only for ages
mod_onlyAge = model.matrix(~age_death, data=pheno)
mod0_onlyAge = model.matrix(~1, data=pheno)

# Save as RDS files
saveRDS(mod_onlyAge, file = sprintf("%s%s_svaMod_onlyAge.rds",
                                    outDir,
                                    outBaseName))
saveRDS(mod0_onlyAge, file = sprintf("%s%s_svaMod0_onlyAge.rds",
                                     outDir,
                                     outBaseName))

# Create design matrixes for SVA accounting for ages and adjustment variables
# covariates (pmi, sex and substudy)
mod = model.matrix(~., data=pheno)
mod0 = model.matrix(~., data=pheno[, colnames(pheno) != "age_death"])

# Save as RDS files
saveRDS(mod, file = sprintf("%s%s_svaMod_all.rds",
                            outDir,
                            outBaseName))
saveRDS(mod0, file = sprintf("%s%s_svaMod0_all.rds",
                             outDir,
                             outBaseName))

# Create mod and batch objects for running ComBat
################################################################################
batch <- pheno$substudy
modCombat <- model.matrix(~1 + age_death, data = pheno)

saveRDS(modCombat, file = sprintf("%s%s_combatMod.rds", outDir, outBaseName))
saveRDS(batch, file = sprintf("%s%s_batches.rds", outDir, outBaseName))