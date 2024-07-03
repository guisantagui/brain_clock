################################################################################
# Remove tissues from preprocessed TPM file (intended for cerebellum, as it    #
# is very different).                                                          #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)


# Parser
################################################################################
parser <- arg_parser("Remove a set of tissues from input dataset.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--tiss2rem",
                               "--outTag",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata CSV file.",
                                "Tissues to be removed, separated by commas.",
                                "Tag to add to the returned file name",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
preprocTPMsFile <- parsed$input
metDatFile <- parsed$metDat
tiss2rem <- parsed$tiss2rem
outTag <- parsed$outTag
outDir <- parsed$outDir
outName <- sprintf("%s%s_%s.csv", outDir, gsub(".csv", "", basename(preprocTPMsFile)), outTag)

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Load data
################################################################################
metDat <- read.csv(metDatFile, row.names = 1)
preprocTPMs <- data.frame(data.table::fread(preprocTPMsFile))
rownames(preprocTPMs) <- preprocTPMs$V1
preprocTPMs <- preprocTPMs[, colnames(preprocTPMs) != "V1"]

# Filter out samples belonging to queried tisues
################################################################################
tiss2rem <- strsplit(tiss2rem, split = ",")[[1]]


boolVec <- !metDat$tissue[match(rownames(preprocTPMs),
                                make.names(metDat$specimenID))] %in% tiss2rem

preprocTPMs <- preprocTPMs[boolVec, ]

write.csv(preprocTPMs, outName)

print(sprintf("%s removed from input dataset.",
              paste(tiss2rem, collapse = ", ")))

print(sprintf("%s saved in %s",
              basename(outName),
              dirname(outName)))