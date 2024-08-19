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
                               "--subSampSVAmods",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata CSV file.",
                                "Proportion for subsampling the SVA mod objects. If not specified all samples will be used. Subsampling will be done to ensure that all the substudies are represented. This is done because if the dataset is very big it might take a lot of time to run SVA",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F),
                       default = list("input" = NA,
                                      "--metDat" = NA,
                                      "--subSampSVAmods" = 1,
                                      "--outDir" = NA))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
inFile <- parsed$input
metDatFile <- parsed$metDat
subSampSVAmods <- as.numeric(parsed$subSampSVAmods)
outDir <- parsed$outDir

outBaseName <- gsub(".csv", "", basename(inFile))

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

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

pheno <- data.frame(specimenID = rownames(input),
                    age_death = metDat$ageDeath[match(rownames(input),
                                                      make.names(metDat$specimenID))],
                    pmi = metDat$pmi[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    sex = metDat$sex[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    substudy = metDat$substudy[match(rownames(input),
                                                     make.names(metDat$specimenID))])

# Obtain vector of sampled specimenIDs for subsampling. If proportion is not 
# specified all samples will be used
################################################################################

if(subSampSVAmods < 1){
        sampledIDsAll <- c()
        for(sbstdy in unique(pheno$substudy)){
                pheno_sbstdy <- pheno[pheno$substudy == sbstdy, ]
                nSamp <- round(nrow(pheno_sbstdy) * subSampSVAmods)
                sampdIDs <- sample(pheno_sbstdy$specimenID, nSamp)
                sampledIDsAll <- c(sampledIDsAll, sampdIDs)
        }
}else{
        sampledIDsAll <- pheno$specimenID
}

# Create mod objects for running SVA
################################################################################

# Substutute NAs in PMI for the median value
pheno$pmi[is.na(pheno$pmi)] <- median(pheno$pmi[!is.na(pheno$pmi)])

# Create design matrixes for SVA accounting only for ages
mod_onlyAge = model.matrix(~age_death, data=pheno[match(sampledIDsAll,
                                                        pheno$specimenID), ])
mod0_onlyAge = model.matrix(~1, data=pheno[match(sampledIDsAll,
                                                 pheno$specimenID), ])

rownames(mod_onlyAge) <- sampledIDsAll
rownames(mod0_onlyAge) <- sampledIDsAll

# Save as RDS files
saveRDS(mod_onlyAge, file = sprintf("%s%s_svaMod_onlyAge.rds",
                                    outDir,
                                    outBaseName))
saveRDS(mod0_onlyAge, file = sprintf("%s%s_svaMod0_onlyAge.rds",
                                     outDir,
                                     outBaseName))

# Create design matrixes for SVA accounting for ages and adjustment variables
# covariates (pmi, sex and substudy)
mod = model.matrix(~., data=pheno[match(sampledIDsAll,
                                        pheno$specimenID),
                                  !colnames(pheno) %in% c("specimenID",
                                                          "batch_seq")])
mod0 = model.matrix(~., data=pheno[match(sampledIDsAll,
                                         pheno$specimenID),
                                   !colnames(pheno) %in% c("age_death",
                                                           "specimenID",
                                                           "batch_seq")])

rownames(mod) <- sampledIDsAll
rownames(mod0) <- sampledIDsAll

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