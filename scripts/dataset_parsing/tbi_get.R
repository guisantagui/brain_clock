################################################################################
# Brain clock: Download TBI data and save it as a CSV file.                    #
################################################################################

if (!require("BiocManager", quietly = T)){
        install.packages("BiocManager")
}
if (!require("biomaRt", quietly = T)){
        BiocManager::install("biomaRt", update = F)
}
if (!require("org.Hs.eg.db", quietly = T)){
        BiocManager::install("org.Hs.eg.db", update = F)
}
library(biomaRt)
library(org.Hs.eg.db)


#if(!require("mygene", quietly = T)) BiocManager::install("mygene", update = F)
#library(mygene)
#if(!require("rentrez", quietly = T)) install.packages("rentrez")
#library(mygene)

# Directory stuff
################################################################################
baseURL <- "https://aging.brain-map.org"

setwd("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock_paper_repo/brain_clock/scripts/dataset_parsing/")
tbiURLFile <- "../../data/dataset_parsing/tbi_data_files.csv"
#outDir <- "../../results/dataset_parsi/tbi_raw/"
outDir <- sprintf("%s/tbi_raw/", dirname(tbiURLFile))
parsedDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/"
metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat.csv"

#gene2ensemblURL <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
#gene2ensemblFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/gene2ensembl/gene2ensembl"

#gene2ensembl <- read.table(gene2ensemblFile, sep = "\t", header = F)
#colnames(gene2ensembl) <- c()

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}


# Functions
################################################################################

# This function accepts a vector of ensembl IDs and returns info about 
# the updated version, symbols, start pos, end pos and size
getENSEMBL_ID_info <- function(ensembl_list,
                               mart,
                               removeNAs = T,
                               filt = "ensembl_gene_id_version"){
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
        }else if(filt == "ensembl_gene_id" | filt == "entrezgene_id"){
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
                # pass the original entrez ids to the entrezgene_id column if they don't have id there
                boolCompEntrezIDs <- !is.na(outDF$entrezgene_id_orig) & is.na(outDF$entrezgene_id)
                outDF$entrezgene_id[boolCompEntrezIDs] <- outDF$entrezgene_id_orig[boolCompEntrezIDs]
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

# Function to fetch ENSEMBL ID for a single Entrez ID
fetch_ensembl_id <- function(entrez_id) {
        links <- entrez_link(dbfrom = "gene", id = entrez_id, db = "ensembl")
        if (!is.null(links$links$gene_ensembl)) {
                return(links$links$gene_ensembl[1])
        } else {
                return(NA)
        }
}
#tst2 <- entrez_link(dbfrom = "gene", id = notMapped[1], db = "all")
#tst3 <- entrez_link(dbfrom = "gene", id = "1", db = "all")

#lapply(tst3$links, function(x) x)

#tst2$links$ense
#fetch_ensembl_id(notMapped[1])

#notMappedInfo <- getGenes(notMapped, fields = "all")


# Parse metadata
################################################################################
#metDat <- read.csv(metDatFile)

tbiURLs <- read.csv(tbiURLFile)
# Check duplicates in specimen name
dupSpecs <- tbiURLs$specimen_name[duplicated(tbiURLs$specimen_name)]
tbiURLs[tbiURLs$specimen_name %in% dupSpecs, ]
tbiURLs$specimen_name[tbiURLs$specimen_name %in% dupSpecs]
tbiURLs$structure_acronym[tbiURLs$specimen_name %in% dupSpecs]

# All the specimen_name duplicates have PCx, while in structure acronym they 
# don't have it: half of them are PCx (parietal neocortex) and half FWM (white
# matter of the forebrain). It seems looks like an error, so let's correct the
# strucure acronym in the specimen_name
tbiURLs$specimen_name <- sapply(1:nrow(tbiURLs),
                                function(x) gsub(strsplit(tbiURLs$specimen_name[x],
                                                          split = ".",
                                                          fixed = T)[[1]][4],
                                                 tbiURLs$structure_acronym[x],
                                                 tbiURLs$specimen_name[x]))

# Download the files
################################################################################
tabList <- list()
#countMat <- data.frame()
for(i in 1:nrow(tbiURLs)){
        specimen <- tbiURLs$specimen_name[i]
        #struc <- tbiURLs$structure_acronym[i]
        print(specimen)
        print(sprintf("%s %%", as.character(i/nrow(tbiURLs) * 100)))
        #fileName <- sprintf("%s%s_%s.txt", outDir, specimen, struc)
        fileName <- sprintf("%s%s.txt", outDir, specimen)
        urlEnd <- tbiURLs$gene_level_fpkm_file_link[i]
        if(!file.exists(fileName)){
                comm <- sprintf("curl -o %s %s%s", fileName, baseURL, urlEnd)
                system(comm)
                if(!file.exists(fileName)){
                        print(sprintf("%s could not be downloaded",
                                      basename(fileName)))
                        tab <- NA
                }else{
                        tab <- read.table(fileName, sep = "\t", header = T)
                }
        }else{
                print(sprintf("%s already exists.", basename(fileName)))
                tab <- read.table(fileName, sep = "\t", header = T)
        }
        tabList[[specimen]] <- tab
}

# Generate the counts matrix
################################################################################

# Keep the matrices that were downloaded successfully
tabList <- tabList[!unlist(lapply(tabList, function(x) all(is.na(x))))]

# Obtain unique genes from all the files 
uniqueGenes <- unique(unlist(lapply(tabList, function(x) x$gene_id)))

# Create the counts matrix
countMat <- matrix(nrow = length(tabList),
                   ncol = length(uniqueGenes),
                   dimnames = list(names(tabList),
                                   make.names(uniqueGenes)))
for(i in seq_along(tabList)){
        samp <- names(tabList)[i]
        df <- tabList[[samp]]
        countMat[samp,
                 match(make.names(df$gene_id),
                       colnames(countMat))] <- df$expected_count
}

countMat[1:10, 1:10]

# Check if there are NAs in some sample or features
sampsWNAs <- rownames(countMat)[apply(countMat, 1, function(x) any(is.na(x)))]
featsWNAs <- colnames(countMat)[apply(countMat, 2, function(x) any(is.na(x)))]
countMat[sampsWNAs, featsWNAs]

# H14.09.022.TCx.01_TCx has missing values in 22769, so let's just remove this
# sample
countMat <- countMat[rownames(countMat) != sampsWNAs, ]

# Get the gene info for the colnames of countMat

# Obtain gene info dataframe with bioMart
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl",
                 host = "https://www.ensembl.org")
geneInfo <- getENSEMBL_ID_info(gsub("X", "", colnames(countMat)),
                               mart = human,
                               removeNAs = T,
                               filt = "entrezgene_id")

notMapped <- gsub("X", "", colnames(countMat))[!gsub("X", "", colnames(countMat)) %in% geneInfo$entrezgene_id]


# Convert colnames of countMat to ensembl IDs
ensemblColnames <- geneInfo$ensembl_gene_id[match(colnames(countMat),
                                                  make.names(geneInfo$entrezgene_id_orig))]

# Remove columns whose names are NAs after conversion to ensembl.
countMat <- countMat[, !is.na(ensemblColnames)]
colnames(countMat) <- ensemblColnames[!is.na(ensemblColnames)]

# Check if after rename there are duplicates
ensemblDups <- unique(colnames(countMat)[duplicated(colnames(countMat))])
geneInfo_dups <- geneInfo[geneInfo$ensembl_gene_id %in% ensemblDups, ]
geneInfo_dups <- geneInfo_dups[order(geneInfo_dups$ensembl_gene_id), ]

countMat_dups <- countMat[, colnames(countMat) %in% ensemblDups]

countMat_dups[, order(colnames(countMat_dups))]


# Get the sum of the counts of the genes that after conversion have the same 
# ensembl ID.
for(i in seq_along(ensemblDups)){
        #i <- 1
        dup <- ensemblDups[i]
        dupMat <- countMat[, colnames(countMat) == dup]
        sumsMat <- matrix(rowSums(dupMat),
                          nrow = nrow(dupMat),
                          ncol = 1,
                          dimnames = list(rownames(dupMat),
                                          dup))
        countMat <- countMat[, colnames(countMat) != dup]
        countMat <- cbind(countMat, sumsMat)
}
unique(colnames(countMat)[duplicated(colnames(countMat))])

# Check the overlap with the integrated TPM matrix
dataFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts.csv"

dat <- data.frame(data.table::fread(dataFile))
rownames(dat) <- dat$V1
dat <- dat[, colnames(dat) != "V1"]

sum(colnames(countMat) %in% colnames(dat))/ncol(countMat)
sum(colnames(countMat) %in% colnames(dat))/ncol(dat)

dim(dat)
dim(countMat)
length(intersect(colnames(countMat), colnames(dat)))

# Load the model to see if the features significant in the training are among 
# the ones kept.
mod <- readRDS("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/sascha_procs/modFuncsAlpha0.5.rds")

modelGenes <- mod@model$coefficients_table$names[mod@model$coefficients_table$coefficients != 0 & mod@model$coefficients_table$names != "Intercept"]

sum(!modelGenes %in% geneInfo$symbol)
modelGenes[!modelGenes %in% geneInfo$symbol]
# Only two genes are left out, so let's try to merge the datasets and see how preprocessing and training
# with these samples works out.
commGenes <- intersect(colnames(countMat), colnames(dat))

counts_tbi_rest <- rbind.data.frame(data.frame(countMat[, commGenes]),
                                    dat[, commGenes])

dim(counts_tbi_rest)

# Create a metadata file of the TBI samples
metDat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat.csv", row.names = 1)
tbiMetDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/DonorInformation.csv"
tbiMetDat <- read.csv(tbiMetDatFile)

tbiMetDat_parsed <- data.frame(specimenID = rownames(countMat),
                               platform = rep(NA, nrow(countMat)),
                               RIN = tbiURLs$rna_integrity_number[match(rownames(countMat),
                                                                        tbiURLs$specimen_name)],
                               libraryPrep = rep(NA, nrow(countMat)),
                               libraryPreparationMethod = rep(NA, nrow(countMat)),
                               runType = rep(NA, nrow(countMat)),
                               readLength = rep(NA, nrow(countMat)),
                               individualID = sapply(rownames(countMat),
                                                     function(x) paste(strsplit(x,
                                                                                split = ".",
                                                                                fixed = T)[[1]][1:3],
                                                                       collapse = ".")),
                               organ = rep("brain", nrow(countMat)),
                               tissue = tbiURLs$structure_name[match(rownames(countMat),
                                                                     tbiURLs$specimen_name)],
                               assay = rep("rnaSeq", nrow(countMat)),
                               exclude = rep(NA, nrow(countMat)),
                               excludeReason = rep(NA, nrow(countMat)))

tbiMetDat_parsed$sex <- tbiMetDat$sex[match(tbiMetDat_parsed$individualID,
                                            tbiMetDat$name)]

tbiMetDat_parsed$sex <- gsub("M", "male", tbiMetDat_parsed$sex)
tbiMetDat_parsed$sex <- gsub("F", "female", tbiMetDat_parsed$sex)
tbiMetDat_parsed$race <- tbiMetDat$race[match(tbiMetDat_parsed$individualID,
                                              tbiMetDat$name)]
tbiMetDat_parsed$ageDeath <- tbiMetDat$age[match(tbiMetDat_parsed$individualID,
                                                  tbiMetDat$name)]

tbiMetDat_parsed$apoeGenotype <- tbiMetDat$apo_e4_allele[match(tbiMetDat_parsed$individualID,
                                                               tbiMetDat$name)]

tbiMetDat_parsed$pmi <- rep(NA, nrow(tbiMetDat_parsed))
tbiMetDat_parsed$Braak <- tbiMetDat$braak[match(tbiMetDat_parsed$individualID,
                                                tbiMetDat$name)]

# Get the diagnosis for the training. We consider controls the samples
# with a braak index of less than 4 and a cerad score of 0.

tbiDiagnConsensus <- ((tbiMetDat$dsm_iv_clinical_diagnosis[match(tbiMetDat_parsed$individualID, tbiMetDat$name)] == "No Dementia") & (tbiMetDat$nincds_arda_diagnosis[match(tbiMetDat_parsed$individualID, tbiMetDat$name)] == "No Dementia"))
tbiMetDat_parsed$diagn_4BrainClck <- (tbiMetDat_parsed$Braak <= 3 & tbiMetDat$cerad[match(tbiMetDat_parsed$individualID, tbiMetDat$name)] == 0) & tbiDiagnConsensus
tbiMetDat_parsed$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(tbiMetDat_parsed$diagn_4BrainClck) + 1]

table(tbiMetDat_parsed$Braak[tbiMetDat_parsed$diagn_4BrainClck == "Control"])

tbiMetDat_parsed$substudy <- rep("TBI", nrow(tbiMetDat_parsed))
tbiMetDat_parsed$batch_rna <- rep(NA, nrow(tbiMetDat_parsed))
tbiMetDat_parsed$batch_lib <- rep(NA, nrow(tbiMetDat_parsed))
tbiMetDat_parsed$batch_seq <- rep(NA, nrow(tbiMetDat_parsed))
tbiMetDat_parsed$mmse30_first_ad_dx <- rep(NA, nrow(tbiMetDat_parsed))
tbiMetDat_parsed$mmse30_lv <- rep(NA, nrow(tbiMetDat_parsed))

# Parse age column
# Remove the rows that have plus in ages
tbiMetDat_parsed <- tbiMetDat_parsed[!grepl("+",
                                            tbiMetDat_parsed$ageDeath,
                                            fixed = T), ]

tbiMetDat_parsed$ageDeath <- sapply(tbiMetDat_parsed$ageDeath,
                                     function(x){
                                             if(grepl("-", x, fixed = T)){
                                                     mean(as.numeric(strsplit(x,
                                                                              split = "-",
                                                                              fixed = T)[[1]]))
                                             }else{
                                                     as.numeric(x)
                                             }
                                     })

# Remove the samples that have traumatic brain injuries.
tbiMetDat_parsed <- tbiMetDat_parsed[tbiMetDat$age_at_first_tbi[match(tbiMetDat_parsed$individualID,
                                                                      tbiMetDat$name)] == 0, ]#$diagn_4BrainClck

# Merge metadata with the rest
metDat_merged <- rbind.data.frame(metDat, tbiMetDat_parsed)

# Remove from the dounts dataframe the samples that are not in the metadata
# dataframe
counts_tbi_rest <- counts_tbi_rest[rownames(counts_tbi_rest) %in% make.names(metDat_merged$specimenID), ]

dim(metDat_merged)
dim(counts_tbi_rest)

metDat_merged[!make.names(metDat_merged$specimenID) %in% rownames(counts_tbi_rest), ]

ncol(metDat)
ncol(tbiMetDat_parsed)

# There are some samples from Mayo study that are not in the counts DF. These
# samples were not present in the Mayo_gene_all_counts_matrix_clean.txt that
# we obtained from synapse.org. So no problem.
metDat_merged[!make.names(metDat_merged$specimenID) %in% rownames(counts_tbi_rest), ]

write.csv(tbiMetDat_parsed, file = paste0(parsedDir, "TBI_metadata_unified.csv"))
write.csv(metDat_merged, file = paste0(parsedDir, "combined_metDat_wTBI.csv"))
write.csv(counts_tbi_rest, file = paste0(parsedDir, "combined_counts_wTBI.csv"))


entrez_to_symb_ncbi <- function(entrez){
        entrez <- notMapped[1:5]
        entrezCollps <- paste(entrez, collapse = ",")
        print(entrezCollps)
        rootDir <- outDir
        tempDir <- paste0(rootDir, "temp/")
        if(!dir.exists(tempDir)){
                dir.create(tempDir)
        }
        
        dwComm <- sprintf("datasets download gene gene-id %s --filename %sgenes.zip",
                          entrezCollps,
                          tempDir)
        system(dwComm)
        
        commUnzip <- sprintf("unzip -o %sgenes.zip -d %s",
                             tempDir,
                             tempDir)
        
        system(commUnzip)
        
        jsonFile <- paste0(tempDir, "ncbi_dataset/data/data_report.jsonl")
        
        jsonDF <- as.data.frame(ndjson::stream_in(jsonFile, cls = "tbl"))
        
        synonyms <- jsonDF[, "symbol"]
        
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