################################################################################
#                                                                              #
# Filters genes to keep only the ones in a subset dataset. Intended for        #
# considering only LINCS landmark genes, or protein coding genes               #
#                                                                              #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Run SVA.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--filtFile",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "File with ENSEMBL IDs of genes that are going to be kept, in a column called ensembl_gene_id_orig.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Terminal argument parser
################################################################################

inFile <- parsed$input
filtFile <- parsed$filtFile
outDir <- parsed$outDir

outName <- sprintf("%s%s_%s",
                   outDir,
                   gsub(".csv", "", basename(inFile)),
                   basename(filtFile))

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

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

# Load data
################################################################################
print(sprintf("Filtering %s to keep only genes included in %s",
              basename(inFile),
              basename(filtFile)))

df <- readCsvFast(inFile)
filtDF <- readCsvFast(filtFile)

# Filter the dataset
################################################################################

keepGenes <- intersect(colnames(df), filtDF$ensembl_gene_id)
df <- df[, keepGenes]

# Save the dataset
################################################################################
writeCsvFst(df, file = outName)
print(sprintf("%s saved at %s.",
              basename(outName),
              dirname(outName)))