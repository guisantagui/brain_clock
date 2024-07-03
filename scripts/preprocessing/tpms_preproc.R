################################################################################
# Preprocess merged TPM file, for removing variables with high proportion of.  #
# zeros and log-transforming it.                                               #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Parser
################################################################################
parser <- arg_parser("Preprocess merged TPM file to deal with zeros, and return preprocessed and log10 transformed file.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--propZerosRem",
                               "--rinFilt",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata CSV file.",
                                "Threshold proportion of zeros within samples and features to drop it out of the dataset.",
                                "Minimum RIN allowed for the samples to be considered",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
tpmFile <- parsed$input
metDatFile <- parsed$metDat
propZeroRem <- as.numeric(parsed$propZerosRem)
rinFilt <- as.numeric(parsed$rinFilt)
outDir <- parsed$outDir

outName <- gsub(".csv", "", basename(tpmFile))


if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Load data
################################################################################
metDat <- read.csv(metDatFile, row.names = 1)
comb_TPMs <- data.frame(data.table::fread(tpmFile))
rownames(comb_TPMs) <- comb_TPMs$V1
comb_TPMs <- comb_TPMs[, colnames(comb_TPMs) != "V1"]

# Preprocessing
################################################################################
plot(density(as.matrix(comb_TPMs)))
plot(density(as.matrix(log10(comb_TPMs + 0.00001))))


# Check proportion of zeros
printZeroProp <- function(df){
        zeroProp <- sum(df == 0)/(nrow(df)*ncol(df))
        
        print(sprintf("Proportion of zeros: %s",
                      as.character(round(zeroProp,
                                         digits = 4))))
}

printZeroProp(comb_TPMs)

# There is a large proportion of values that are zeros (43%). So let's iterate
# through each age bin and identify the genes that do not have at least 1 TPM
# in all the samples for each age group. Once identified, we will remove from
# the dataset the genes that are common across all ages.
lowExpGenes_perAge <- list()
for(age in unique(metDat$ageDeath)){
        #age <- 67
        sampsAge <- metDat$specimenID[metDat$ageDeath == age]
        sampsAge <- make.names(sampsAge)
        sampsAge <- sampsAge[sampsAge %in% rownames(comb_TPMs)]
        tpms_age <- comb_TPMs[make.names(sampsAge), ]
        #make.names(sampsAge) %in% row.names(comb_TPMs)
        #any(is.na(tpms_age))
        #tpms_age[1:100, 1:10]
        lowExpGenes <- colnames(tpms_age)[colSums(tpms_age) < nrow(tpms_age)]
        lowExpGenes_perAge[[paste("age",
                                  as.character(age),
                                  sep = "_")]] <- lowExpGenes
}

# Check if there are NAs
names(lowExpGenes_perAge)[unlist(lapply(lowExpGenes_perAge, function(x) all(is.na(x))))]

# Look for the intersection across all ages.
recurs_intersect <- function(lst){
        if (length(lst) == 1) {
                return(lst[[1]])
        }
        
        # Recursive case: find the intersection of the first element with the
        # intersection of the rest of the list
        return(intersect(lst[[1]], recurs_intersect(lst[-1])))
}

lowExpGenes_intersect <- recurs_intersect(lowExpGenes_perAge)

# Remove these low-expressed genes across all ages from the TPMs DF
comb_TPMs <- comb_TPMs[, !colnames(comb_TPMs) %in% lowExpGenes_intersect]

printZeroProp(comb_TPMs)

plot(density(as.matrix(log10(comb_TPMs + 0.00001))))

# Still, data is quite sparse: 26% of values are zero

# Check variables that have 80% or more of the values zero
varsHighPropZero <- colnames(comb_TPMs)[apply(comb_TPMs,
                                              2,
                                              function(x) sum(x == 0)/length(x)) > propZeroRem]

sampsHighPropZero <- rownames(comb_TPMs)[apply(comb_TPMs,
                                               1,
                                               function(x) sum(x == 0)/length(x)) > propZeroRem] 


comb_TPMs <- comb_TPMs[!rownames(comb_TPMs) %in% sampsHighPropZero,
                       !colnames(comb_TPMs) %in% varsHighPropZero]
printZeroProp(comb_TPMs)
plot(density(as.matrix(log10(comb_TPMs + 0.00001))))

# Log-transform 
comb_TPMs_log <- log10(comb_TPMs + 0.00001)
plot(density(as.matrix(comb_TPMs_log)))

# Remove samples with a RIN below of the user-specified threshold
# Filter out samples with RIN < 6
sampsGoodRIN <- make.names(metDat$specimenID[metDat$RIN >= rinFilt])
comb_TPMs <- comb_TPMs[rownames(comb_TPMs) %in% sampsGoodRIN, ]
comb_TPMs_log <- comb_TPMs_log[rownames(comb_TPMs_log) %in% sampsGoodRIN, ]

# Save preprocessed datasets
################################################################################
prepOutName <- sprintf("%s%s_preproc.csv", outDir, outName)
write.csv(comb_TPMs, prepOutName)
print(sprintf("%s saved at %s", basename(prepOutName), dirname(prepOutName)))

logOutName <- sprintf("%s%s_preproc_log10.csv", outDir, outName)
write.csv(comb_TPMs_log, logOutName)
print(sprintf("%s saved at %s", basename(logOutName), dirname(logOutName)))