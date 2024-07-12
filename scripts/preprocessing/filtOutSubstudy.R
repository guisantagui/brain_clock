################################################################################
#                                                                              #
# Script for filtering out samples belonging to a particular substudy.         #
#                                                                              #
################################################################################

# Terminal argument parser
################################################################################
parser <- arg_parser("Remove substudy from input data.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--substud2Rem"
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Metadata file.",
                                "Substudy to be filtered out.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

inFile <- parsed$input
metDatFile <- parsed$metDatFile
substud2Rem <- parsed$substud2Rem
outDir <- parsed$outDir

outName <- sprintf("%s%s_no%s.csv",
                   outDir,
                   gsub(".csv", "", basename(inFile)),
                   substud2Rem)

# Functions
################################################################################

# Read csv faster
readCsvFast <- function(f, header = T){
        df <- data.frame(data.table::fread(f, header = header))
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

dat <- readCsvFast(inFile)
metDat <- readCsvFast(metDatFile)

# Filter out samples belonging to the selected study
################################################################################
print(sprintf("Filtering out samples belinging to %s substudy...",
              substud2Rem))
sampsKeep <- make.names(metDat$specimenID[metDat$substudy != substud2Rem])

dat <- dat[sampsKeep, ]
writeCsvFst(dat, file = outName)
print(sprintf("%s saved at %s.",
              basename(outName),
              dirname(outName)))
