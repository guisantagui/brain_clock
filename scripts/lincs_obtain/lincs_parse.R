################################################################################
# Brain clock: obtain expression matrix and metadata file from LINCS gctx file.#
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("cmapR", quietly = T)) BiocManager::install("cmapR", update = F)
library(cmapR)
if(!require("biomaRt", quietly = T)) BiocManager::install("biomaRt", update = F)
library(biomaRt)
if(!require("org.Hs.eg.db", quietly = T)) BiocManager::install("org.Hs.eg.db",
                                                               update = F)
library(org.Hs.eg.db)

# Terminal argument parser
################################################################################
parser <- arg_parser("Read LINCS gctx filtered files and generate expression matrices and metadata files.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--outDir"),
                       help = c("Input GCTX file.",
                                "Output directory where placing the results."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

lincs_file <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/gctx/cp_predicted_RNAseq_profiles_NPC_2.gctx"
outDir <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/"

lincs_file <- parsed$input
outDir <- parsed$outDir

outExpName <- sprintf("%s%s_expMat.csv", outDir, gsub(".gctx", "", basename(lincs_file)))
outMetDatName <- sprintf("%s%s_metDat.csv", outDir, gsub(".gctx", "", basename(lincs_file)))

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

isCtl <- grepl("ctl", basename(lincs_file))

# Functions
################################################################################

# This function accepts a vector of ensembl IDs and returns info about 
# the updated version, symbols, start pos, end pos and size
getENSEMBL_ID_info <- function(ensembl_list,
                               mart,
                               removeNAs = T,
                               filt = "ensembl_gene_id_version",
                               remapUnmap = T){
        #ensembl_list <- gsub("X", "", colnames(countMat))
        #mart <- human
        #filt <- "entrezgene_id"
        
        
        print(sprintf("Mapping the %s...", filt))
        if(filt == "ensembl_gene_id_version"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "ensembl_gene_id_version",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }else if(filt == "ensembl_gene_id" | filt == "entrezgene_id" | filt == "hgnc_symbol"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }
        
        gene_coords <- getBM(attributes = atts,
                             filters = filt,
                             values = ensembl_list,
                             mart = mart)
        
        gene_coords$size <- gene_coords$end_position - gene_coords$start_position
        if(filt == "ensembl_gene_id_version"){
                mapped_matchVec <- match(ensembl_list,
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
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$ensembl_gene_id)
                outDF <- data.frame(ensembl_gene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
        }else if(filt == "entrezgene_id"){
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$entrezgene_id)
                outDF <- data.frame(entrezgene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    entrezgene_id = gene_coords$entrezgene_id[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
                unMapped <- outDF$entrezgene_id_orig[is.na(outDF$entrezgene_id)]
                # Try out with org.Hs.eg.db for the unmapped ones
                unMapped_ensembl <- mapIds(org.Hs.eg.db,
                                           keys = unMapped,
                                           keytype = "ENTREZID",
                                           column = c("ENSEMBL"),
                                           fuzzy = T)
                unMapped_ensembl_noNAs <- unMapped_ensembl[!is.na(unMapped_ensembl)]
                gene_coords_unMap <- getBM(attributes = atts,
                                           filters = "ensembl_gene_id",
                                           values = unMapped_ensembl_noNAs,
                                           mart = mart)
                gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                gene_coords_unMap$entrezgene_id_orig <- names(unMapped_ensembl_noNAs)[match(gene_coords_unMap$ensembl_gene_id, unMapped_ensembl_noNAs)]
                outDF[match(gene_coords_unMap$entrezgene_id_orig,
                            outDF$entrezgene_id_orig),
                      c("ensembl_gene_id",
                        "symbol",
                        "entrezgene_id",
                        "chromosome",
                        "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                          "start_position", "end_position", "size")]
                # Repeat mapping of unMapped but this time with symbols
                unMapped <- outDF$entrezgene_id_orig[is.na(outDF$entrezgene_id)]
                unMapped_symbols <- mapIds(org.Hs.eg.db,
                                           keys = unMapped,
                                           keytype = "ENTREZID",
                                           column = c("SYMBOL"),
                                           fuzzy = T)
                unMapped_symbols_noNAs <- unMapped_symbols[!is.na(unMapped_symbols)]
                gene_coords_unMap <- getBM(attributes = atts,
                                           filters = "hgnc_symbol",
                                           values = unMapped_symbols_noNAs,
                                           mart = mart)
                gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                gene_coords_unMap$entrezgene_id_orig <- names(unMapped_symbols_noNAs)[match(gene_coords_unMap$hgnc_symbol, unMapped_symbols_noNAs)]
                outDF[match(gene_coords_unMap$entrezgene_id_orig,
                            outDF$entrezgene_id_orig),
                      c("ensembl_gene_id",
                        "symbol",
                        "entrezgene_id",
                        "chromosome",
                        "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                          "start_position", "end_position", "size")]
                
                
                # pass the original entrez ids to the entrezgene_id column if they don't have id there
                boolCompEntrezIDs <- !is.na(outDF$entrezgene_id_orig) & is.na(outDF$entrezgene_id)
                outDF$entrezgene_id[boolCompEntrezIDs] <- outDF$entrezgene_id_orig[boolCompEntrezIDs]
                naVec <- is.na(outDF$ensembl_gene_id)
                
        }else if(filt == "hgnc_symbol"){
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$hgnc_symbol)
                outDF <- data.frame(symbol_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    entrezgene_id = gene_coords$entrezgene_id[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                # If selected, another round of mapping is done on the genes 
                # that are unmapped using org.Hs.eg.db, to maximize number of
                # maps.
                if(remapUnmap){
                        #naVec <- is.na(outDF$ensembl_gene_id)
                        unMapped <- outDF$symbol_orig[is.na(outDF$entrezgene_id)]
                        unMapped_ensembl <- mapIds(org.Hs.eg.db,
                                                   keys = unMapped,
                                                   keytype = "SYMBOL",
                                                   column = c("ENSEMBL"),
                                                   fuzzy = T)
                        unMapped_ensembl_noNAs <- unMapped_ensembl[!is.na(unMapped_ensembl)]
                        gene_coords_unMap <- getBM(attributes = atts,
                                                   filters = "ensembl_gene_id",
                                                   values = unMapped_ensembl_noNAs,
                                                   mart = mart)
                        gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                        gene_coords_unMap$symbol_orig <- names(unMapped_ensembl_noNAs)[match(gene_coords_unMap$ensembl_gene_id, unMapped_ensembl_noNAs)]
                        
                        outDF[match(gene_coords_unMap$symbol_orig,
                                    outDF$symbol_orig),
                              c("ensembl_gene_id",
                                "symbol",
                                "entrezgene_id",
                                "chromosome",
                                "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                                  "start_position", "end_position", "size")]
                }
                naVec <- is.na(outDF$size)
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

# This function adds up columns or rows that ended up with same name 
# after transformation to EMSEMBL.
sumDupsMat <- function(m, dups_in){
        if(dups_in == "cols"){
                m <- t(m)
        }
        featsOrd <- unique(rownames(m))
        dups <- unique(rownames(m)[duplicated(rownames(m))])
        
        for(i in seq_along(dups)){
                dup <- dups[i]
                dupMat <- m[rownames(m) == dup, ]
                sumsMat <- matrix(colSums(dupMat),
                                  nrow = 1,
                                  ncol = ncol(dupMat),
                                  dimnames = list(dup,
                                                  colnames(dupMat)))
                m <- m[rownames(m) != dup, ]
                m <- rbind(m, sumsMat)
        }
        m <- m[match(featsOrd, rownames(m)), ]
        if(dups_in == "cols"){
                m <- t(m)
        }
        return(m)
}

# Load data
################################################################################
lincs <- parse_gctx(lincs_file)


# Parse dataset into a table with proper column and row names 
################################################################################

nuColNames <- paste(make.names(lincs@cdesc$id),
                    paste(lincs@cdesc$cell,
                          lincs@cdesc$pertname,
                          gsub(" ", "", lincs@cdesc$timepoint),
                          gsub(" ", "", lincs@cdesc$dose),
                          sep = "_"),
                    sep = "..")
names(nuColNames) <- make.names(lincs@cdesc$id)

# Order nuColNames just in case
nuColNames <- nuColNames[match(make.names(colnames(lincs@mat)),
                               names(nuColNames))]

# Obtain gene info file, if it doesn't exist, to convert gene symbols into
# ensembl IDs.
geneInfoFile <- sprintf("%slincs_geneInfo.csv", outDir)
if(!file.exists(geneInfoFile)){
        # Convert gene symbols to ensembl IDs
        human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        
        symbList <- lincs@rdesc$id
        
        geneInfo <- getENSEMBL_ID_info(symbList, human, removeNAs = T,
                                       filt = "hgnc_symbol",
                                       remapUnmap = F)
        write.csv(geneInfo, file = geneInfoFile)
}else{
        geneInfo <- read.csv(geneInfoFile, row.names = 1)
}

lincs_mat <- lincs@mat
nuRowNames <- geneInfo$ensembl_gene_id[match(rownames(lincs_mat),
                                             geneInfo$symbol_orig)]

# Remove rows that have NAs in the names:
lincs_mat <- lincs_mat[!is.na(nuRowNames), ]
nuRowNames <- nuRowNames[!is.na(nuRowNames)]
rownames(lincs_mat) <- nuRowNames
colnames(lincs_mat) <- nuColNames


# Check if after rename there are duplicates
dups <- nuRowNames[(duplicated(nuRowNames))]
if(length(dups) == 0){
        print("After converting to ENSEMBL IDs there are no duplicates.")
}
if(length(dups) > 0){
        print(sprintf("After converting to ENSEMBL IDs there are %s duplicates. Adding rows with the same resulting ID together.",
                      as.character(length(dups))))
        lincs_mat <- sumDupsMat(lincs_mat, dups_in = "rows")
}

metDat <- lincs@cdesc
metDat$specimenID <- colnames(lincs_mat)
if(isCtl){
        metDat$group <- rep("ctrl", nrow(metDat))
}else{
        metDat$group <- rep("expr", nrow(metDat))
}

write.csv(lincs_mat, outExpName)
write.csv(metDat, outMetDatName)