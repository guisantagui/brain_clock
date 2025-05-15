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
                               "--whatData",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata CSV file.",
                                "Proportion for subsampling the SVA mod objects. If not specified all samples will be used. Subsampling will be done to ensure that all the substudies are represented. This is done because if the dataset is very big it might take a lot of time to run SVA",
                                "The data the objects are for: 'clin' clinical (postmortem) or 'pert' perturbation.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F),
                       default = list("input" = NA,
                                      "--metDat" = NA,
                                      "--subSampSVAmods" = 1,
                                      "--whatData" = NA,
                                      "--outDir" = NA))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
inFile <- parsed$input
metDatFile <- parsed$metDat
subSampSVAmods <- as.numeric(parsed$subSampSVAmods)
whatData <- parsed$whatData
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

metDat$tissue <- gsub("NPCs", "NPC", metDat$tissue)
metDat$tissue <- gsub("Cultured cortical neurons,", "NEU", metDat$tissue)
metDat$tissue <- gsub("GABAergic neurons", "NEU", metDat$tissue)
metDat$tissue <- gsub("human neurons", "NEU", metDat$tissue)
metDat$tissue <- gsub("Neuron", "NEU", metDat$tissue)
metDat$tissue <- gsub("human motor neurons", "NEU", metDat$tissue)
metDat$tissue <- gsub("U5-NPC (G1 Phase)", "NPC", metDat$tissue)

pheno <- data.frame(specimenID = rownames(input),
                    age_death = metDat$ageDeath[match(rownames(input),
                                                      make.names(metDat$specimenID))],
                    pmi = metDat$pmi[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    sex = factor(metDat$sex[match(rownames(input),
                                                  make.names(metDat$specimenID))]),
                    substudy = factor(metDat$substudy[match(rownames(input),
                                                            make.names(metDat$specimenID))]),
                    exper_group = factor(metDat$exper_group[match(rownames(input),
                                                                  make.names(metDat$specimenID))]),
                    tissue = factor(metDat$tissue[match(rownames(input),
                                                        make.names(metDat$specimenID))]))

# Obtain vector of sampled specimenIDs for subsampling. If proportion is not 
# specified all samples will be used
################################################################################

if(subSampSVAmods < 1){
        sampledIDsAll <- c()
        seed <- 567
        for(sbstdy in unique(pheno$substudy)){
                pheno_sbstdy <- pheno[pheno$substudy == sbstdy, ]
                nSamp <- round(nrow(pheno_sbstdy) * subSampSVAmods)
                sampdIDs <- sample(pheno_sbstdy$specimenID, nSamp)
                sampledIDsAll <- c(sampledIDsAll, sampdIDs)
                seed <- seed + 1
        }
}else{
        sampledIDsAll <- pheno$specimenID
}

pheno <- pheno[match(sampledIDsAll, pheno$specimenID), ]
pheno$tissue <- droplevels(pheno$tissue)

# Create mod objects for running SVA
################################################################################

if (whatData == "clin"){
        # Substutute NAs in PMI for the median value
        pheno$pmi[is.na(pheno$pmi)] <- median(pheno$pmi[!is.na(pheno$pmi)])

        # Substutute NAs in ageDeath for the median value (they are the perts_111,
        # so they dont actually have an age).
        pheno$age_death[is.na(pheno$age_death)] <- median(pheno$age_death[!is.na(pheno$age_death)])

        # Create design matrixes for SVA accounting only for ages
        mod_onlyAge = model.matrix(~age_death, data=pheno)
        mod0_onlyAge = model.matrix(~1, data=pheno)

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
        mod = model.matrix(~ age_death + pmi + sex + substudy,
                           data=pheno[, !colnames(pheno) %in% c("specimenID",
                                                                "batch_seq")])
        mod0 = model.matrix(~., data=pheno[, !colnames(pheno) %in% c("age_death",
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
}else if (whatData == "pert"){
        mod = model.matrix(~ exper_group + tissue,
                           data=pheno)
        mod <- mod[, -1]
        mod0 = model.matrix(~1, data=pheno)
        rownames(mod) <- pheno$specimenID
        rownames(mod0) <- pheno$specimenID
        # Save as RDS files
        saveRDS(mod, file = sprintf("%s%s_svaMod.rds",
                                    outDir,
                                    outBaseName))
        saveRDS(mod0, file = sprintf("%s%s_svaMod0.rds",
                                     outDir,
                                     outBaseName))
}



