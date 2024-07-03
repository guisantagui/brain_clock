################################################################################
#                                                                              #
# Regress out the effect from PMI and sex from input dataset using limma.      #
#                                                                              #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
library(limma)

# Parser
################################################################################
parser <- arg_parser("Remove a set of tissues from input dataset.")

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
inputFile <- parsed$input
metDatFile <- parsed$metDat
outDir <- parsed$outDir

outBaseName <- gsub(".csv", "", basename(inputFile))

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Load data
################################################################################
metDat <- read.csv(metDatFile, row.names = 1)
dat <- data.frame(data.table::fread(inputFile))
rownames(dat) <- dat$V1
dat <- dat[, colnames(dat) != "V1"]


edata <- t(dat)

pheno <- data.frame(age_death = metDat$ageDeath[match(rownames(dat),
                                                      make.names(metDat$specimenID))],
                    pmi = metDat$pmi[match(rownames(dat),
                                           make.names(metDat$specimenID))],
                    sex = metDat$sex[match(rownames(dat),
                                           make.names(metDat$specimenID))])

# Substutute NAs in PMI for the median value
pheno$pmi[is.na(pheno$pmi)] <- median(pheno$pmi[!is.na(pheno$pmi)])

# Remove covariate influence (PMI and sex)
################################################################################

# Create design matrixes for SVA accounting only for ages
mod_age = model.matrix(~age_death, data=pheno)

# Create a matrix of the covariates to remove
covMat <- pheno[, colnames(pheno) != "age_death"]

covMod <- model.matrix(~., covMat)



edata_adj <- removeBatchEffect(edata,
                               covariates = covMod[, colnames(covMod) != "(Intercept)"],
                               design = mod_age)

edata_adj <- t(edata_adj)

outName <- sprintf("%s%s_pmiAgeAdj.csv", outDir, outBaseName)
write.csv(edata_adj, file = outName)

print(sprintf("%s saved at %s", basename(outName), dirname(outName)))