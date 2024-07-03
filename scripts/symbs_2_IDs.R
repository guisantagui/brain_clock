################################################################################
# CureMILS project:                                                            #
# TRANSCRIPTOMICS: This script takes as input any dataframe using gene Symbols #
# as either rownames or colnames (needs to be specified in an argument) and    #
# translates them to either to ENTREZ IDs or ENSEMBL. It makes use of the      #
# "datasets" NCBI command line tool to retrieve the IDs from NCBIs site of the #
# gene symbols that, due to being using synonym symbols not included in        #
# org.Hs.eg.db, are not being mapped correctly by mapIds function. This way    #
# the number of NAs after conversion from symbol to IDs is significantly       #
# reduced. Output can be either ENTREZID or ENSEMBL.                           #
# Dependencies: datasets tool, in PATH                                         #
################################################################################

library(argparser)

# Terminal argument parser
################################################################################

parser <- arg_parser("Transformation of gene SYMBOLS to ENTREZ or ENSEMBL, overcoming limitations of the use of synonym symbols")

parser <- add_argument(parser = parser,
                       arg = c("input", "--whatOutID", "--remResDups", "--symbsInCols"),
                       help = c("Input file with either colnames or rownames in gene symbols",
                                "What ID to include in the output file. ENTREZID or ENSEMBL",
                                "If duplicates that migh arise from the use of different synonyms symbols in the input data should be removed from the output",
                                "If the symbols are in the colnames of the input data. If flag not present it is considered that they are in the rows."),
                       flag = c(F, F, T, T))

parsed <- parse_args(parser)


# packrat stuff
################################################################################

#rootDir <- parsed$rootDir

#packRatDir <- paste0(rootDir, "packrat/")
# Activate packrat and install packages, if necessary
#packrat::on(project = packRatDir)
#packrat::on()
if(!require("readxl", quietly = T)){
        install.packages("readxl", repos = 'http://cran.us.r-project.org')
}
library(readxl)
if(!require("ndjson", quietly = T)){
        install.packages("ndjson", repos = 'http://cran.us.r-project.org')
}
if(!require("org.Hs.eg.db", quietly = T)) BiocManager::install("org.Hs.eg.db",
                                                               update = F)
library(org.Hs.eg.db)
if(!require("argparser", quietly = T)){
        install.packages("argparser",
                         repos = 'http://cran.us.r-project.org')
}
library(argparser)
library(org.Hs.eg.db)

# Directory stuff
################################################################################

inFile <- parsed$input
whatOutID <- parsed$whatOutID
remResDups <- parsed$remResDups
symbsInRows <- !parsed$symbsInCols

name4Outs <- basename(inFile)
if(grepl(".xlsx", name4Outs)){
        name4Outs <- gsub(".xlsx", "", name4Outs, fixed = T)
}else if(grepl(".csv", name4Outs)){
        name4Outs <- gsub(".csv", "", name4Outs, fixed = T)
}

rootDir <- gsub("results/DE", "", dirname(inFile))

outDir <- paste0(dirname(inFile), "/")

# Load the data
################################################################################

if(grepl(".xlsx", inFile)){
        DF <- as.data.frame(read_xlsx(inFile))
        rownames(DF) <- DF[, 1]
        DF <- DF[, 2:ncol(DF)]
}else if(grepl(".csv", inFile)){
        DF <- read.csv(file = inFile,
                       row.names = 1)
}
if(!symbsInRows){
        DF <- as.data.frame(t(DF))
}
DF <- DF[rownames(DF) != "Intercept", ]
# Parsing
################################################################################

# Get a DF of equivalence between SYMBOL, ENTREZIDs and ENSEMBL IDs

IDs_DF <- data.frame(SYMBOL = rownames(DF),
                     ENTREZID = mapIds(org.Hs.eg.db,
                                       keys = rownames(DF),
                                       keytype = "SYMBOL",
                                       column = "ENTREZID",
                                       fuzzy = T))


NA_symbs <- IDs_DF$SYMBOL[is.na(IDs_DF$ENTREZID)]

# Try to assign ENTREZ ID to the symbols that have a NA ENTREZ ID when mapping
# them using org.Hs.eg.db

# This function takes as input gene symbols and retrieves the ncbi IDs and 
# main symbol names using datasets ncbi command line tool. Needs to be installed
# and accessible in PATH.
symb_to_entrez_ncbi <- function(symbs){
        symbsCollps <- paste(symbs, collapse = ",")
        print(symbsCollps)
        tempDir <- paste0(rootDir, "/temp/")
        if(!dir.exists(tempDir)){
                dir.create(tempDir)
        }
        
        dwComm <- sprintf("datasets download gene symbol %s --taxon human --filename %sgenes.zip",
                          symbsCollps,
                          tempDir)
        system(dwComm)
        
        commUnzip <- sprintf("unzip -o %sgenes.zip -d %s",
                             tempDir,
                             tempDir)
        
        system(commUnzip)
        
        jsonFile <- paste0(tempDir, "ncbi_dataset/data/data_report.jsonl")
        
        jsonDF <- as.data.frame(ndjson::stream_in(jsonFile, cls = "tbl"))
        
        synonyms <- jsonDF[, grep("synonym", colnames(jsonDF))]
        
        synonyms <- apply(synonyms, 1,
                          function(x) paste(x[!is.na(x)],
                                            collapse = ","))
        jsonDF_red <- jsonDF[, c("geneId", "symbol")]
        jsonDF_red$synonyms <- synonyms
        newSymbsVec <- c()
        ncbiIDVec <- c()
        for(s in symbs){
                grepVec <- which(sapply(jsonDF_red$synonyms,
                                        function(x) s %in% strsplit(x,
                                                                    ",")[[1]]))[1]
                if(length(grepVec) == 0){
                        newSymb <- NA
                        ncbiID <- NA
                }else{
                        newSymb <- jsonDF_red$symbol[grepVec]
                        newNcbiID <- jsonDF_red$geneId[grepVec]
                }
                newSymbsVec <- c(newSymbsVec, newSymb)
                ncbiIDVec <- c(ncbiIDVec, newNcbiID)
        }
        commRmTmp <- sprintf("rm -rf %s",
                             tempDir)
        system(commRmTmp)
        out_DF <- data.frame(symbol = symbs,
                             newSymbol = newSymbsVec,
                             entrezid =ncbiIDVec)
        return(out_DF)
}

if(length(NA_symbs) > 0){
        genesNA_idMapped <- symb_to_entrez_ncbi(NA_symbs)
        IDs_DF$ENTREZID[is.na(IDs_DF$ENTREZID)] <- genesNA_idMapped$entrezid[match(IDs_DF$SYMBOL[is.na(IDs_DF$ENTREZID)],
                                                                                   genesNA_idMapped$symbol)]
}

# Add ENSEMBL IDs column
IDs_DF$ENSEMBL <- mapIds(org.Hs.eg.db,
                         keys = IDs_DF$ENTREZID,
                         keytype = "ENTREZID",
                         column = "ENSEMBL",
                         fuzzy = T)

ENSEMBL_ids <- c()
for(i in seq_along(IDs_DF$ENSEMBL)){
        item <- IDs_DF$ENSEMBL[[i]]
        if(is.null(item)){
                ENSEMBL_ids <- c(ENSEMBL_ids, NA)
        }else if(is.na(item)){
                ENSEMBL_ids <- c(ENSEMBL_ids, NA)
        }else{
                ENSEMBL_ids <- c(ENSEMBL_ids, item)
        }
}
IDs_DF$ENSEMBL <- ENSEMBL_ids

# Create the parsed DF with new rownames, previous removal of those features
# that didn't have a match to the target ID. If indicated in args, remove 
# replicates that might result after the mapping (there are some gene symbols)
# that are considered synonims by NCBI but that appear in the datasets as 
# differet features.
DF_parsed <- DF
newRowNames <- IDs_DF[match(rownames(DF), IDs_DF$SYMBOL), whatOutID]
DF_parsed <- DF_parsed[!is.na(newRowNames), ]
newRowNames <- newRowNames[!is.na(newRowNames)]

if(remResDups){
        DF_parsed <- DF_parsed[!duplicated(newRowNames), ]
        newRowNames <- newRowNames[!duplicated(newRowNames)]
        rownames(DF_parsed) <- newRowNames
}else{
        rownames(DF_parsed) <- make.unique(newRowNames)
}

if(!symbsInRows){
        DF_parsed <- as.data.frame(t(DF_parsed))
}

outFileName <- sprintf("%s%s_%s_IDs.csv",
                       outDir,
                       name4Outs,
                       whatOutID)


print(outFileName)

write.csv(DF_parsed, file = outFileName)
