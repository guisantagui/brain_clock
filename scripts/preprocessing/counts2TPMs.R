################################################################################
# Brain clock: Compute TPMs of the merged counts dataset.                      #
################################################################################
if(!require("biomaRt", quietly = T)) BiocManager::install("biomaRt", update = F)
library(biomaRt)
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Plot surrogate variables from SVA")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--outDir"),
                       help = c("Input counts matrix, features in columns, samples rows. CSV.",
                                "Output directory where placing the results."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

countsFile <- parsed$input
outDir <- parsed$outDir
outName <- sprintf("%s%s", outDir, gsub("counts", "TPMs", basename(countsFile)))
geneInfoOutName <- gsub("TPMs", "geneInfo", outName)

# Functions
################################################################################

# Function for reading csv faster
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

# This function accepts a vector of ensembl IDs and returns info about 
# the updated version, symbols, start pos, end pos and size
getENSEMBL_ID_info <- function(ensembl_list,
                               mart,
                               removeNAs = T,
                               filt = "ensembl_gene_id_version"){
        #ensembl_list <- gsub("\\_.*", "", ensembl_list)
        print("Mapping the ensembl_gene_id_version...")
        if(filt == "ensembl_gene_id_version"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "ensembl_gene_id_version",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }else if(filt == "ensembl_gene_id"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }
        gene_coords <- getBM(attributes = atts,
                             filters = filt,
                             values = gsub("_PAR_Y", "", ensembl_list),
                             mart = mart)
        
        gene_coords$size <- gene_coords$end_position - gene_coords$start_position
        if(filt == "ensembl_gene_id_version"){
                mapped_matchVec <- match(gsub("_PAR_Y", "", ensembl_list),
                                         gene_coords$ensembl_gene_id_version)
                outDF <- data.frame(ensembl_gene_id_version_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    ensembl_gene_id_version_new = gene_coords$ensembl_gene_id_version[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                unMapped <- outDF$ensembl_gene_id_version_orig[is.na(outDF$ensembl_gene_id_version_new)]
                
                unMapped_toMap <- gsub("\\..*", "", unMapped)
                
                print("Mapping the ensembl_gene_id_version IDs that couldn't be mapped to any of those to the respective ensembl_gene_id...")
                gene_coords_unmapped <- getBM(attributes = atts,
                                              filters = "ensembl_gene_id",
                                              values = unMapped_toMap,
                                              mart = mart)
                
                gene_coords_unmapped$size <- gene_coords_unmapped$end_position - gene_coords_unmapped$start_position
                
                gene_coords_unmapped <- gene_coords_unmapped[match(unMapped_toMap, gene_coords_unmapped$ensembl_gene_id), ]
                
                
                isUnmap <- is.na(outDF$ensembl_gene_id_version_new)
                outDF$ensembl_gene_id[isUnmap] <- gene_coords_unmapped$ensembl_gene_id
                outDF$ensembl_gene_id_version_new[isUnmap] <- gene_coords_unmapped$ensembl_gene_id_version
                outDF$symbol[isUnmap] <- gene_coords_unmapped$hgnc_symbol
                outDF$chromosome[isUnmap] <- gene_coords_unmapped$chromosome_name
                outDF$start[isUnmap] <- gene_coords_unmapped$start_position
                outDF$end[isUnmap] <- gene_coords_unmapped$end_position
                outDF$size[isUnmap] <- gene_coords_unmapped$size
                naVec <- is.na(outDF$ensembl_gene_id_version_new)
        }else if(filt == "ensembl_gene_id"){
                mapped_matchVec <- match(gsub("_PAR_Y", "", ensembl_list),
                                         gene_coords$ensembl_gene_id)
                outDF <- data.frame(ensembl_gene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
        }
        if(removeNAs){
                nNAs <- sum(naVec)
                propNAs <- round(nNAs/nrow(outDF) * 100, digits = 2)
                outDF <- outDF[!naVec, ]
                print(sprintf("%s ensembl_gene_ids were removed due to inability to map to ENSEMBL database.", as.character(nNAs)))
                print(sprintf("This represents the %s%% of the total number of ensembl_gene_ids submitted (n = %s)",
                              as.character(propNAs),
                              as.character(length(ensembl_list))))
        }
        return(outDF)
}

# Given a counts matrix with genes in rows and samples in columns, and a 
# geneInfo dataframe containing gene sizes, returns a matrix of TPMs.
counts2TPMs <- function(countDF,
                        geneInfo,
                        IDsIn = "ensembl_gene_id_version_orig",
                        genesInCols = F){
        if(genesInCols){
                inDF <- t(countDF)
        }else{
                inDF <- countDF
        }
        # Filter out genes that are not in geneInfo DF (don't have size)
        inDF <- inDF[rownames(inDF) %in% geneInfo[, IDsIn], ]
        # Get a vector of gene sizes
        sizeVec <- geneInfo$size[match(rownames(inDF), geneInfo[, IDsIn])]
        # Obtain TPMs
        preTPMs <- inDF/sizeVec
        TPMs <- t(t(preTPMs) * 1e6 / colSums(preTPMs))
        if(genesInCols){
                TPMs <- t(TPMs)
        }
        return(TPMs)
}

# Load the data
################################################################################
counts <- readCsvFast(countsFile)

# Get gene info (including sizes) and transform to TPMs
################################################################################

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

geneInfo <- getENSEMBL_ID_info(colnames(counts),
                               mart = human,
                               removeNAs = T,
                               filt = "ensembl_gene_id")

# Transform to TPMs
TPMs <- counts2TPMs(counts,
                    geneInfo = geneInfo,
                    IDsIn = "ensembl_gene_id_orig",
                    genesInCols = T)


# Save TPMs dataset
writeCsvFst(TPMs, outName)
print(sprintf("%s saved at %s.",
              basename(outName),
              dirname(outName)))

# Save geneInfo Dataset
writeCsvFst(geneInfo, geneInfoOutName)
print(sprintf("%s saved at %s.",
              basename(geneInfoOutName),
              dirname(geneInfoOutName)))
