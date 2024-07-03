################################################################################
# CureMILS project:                                                            #
# PROTEOMICS: NEAT Enrichment. The input must be a dataframe (.csv file) where #
# the rownames are ENSEMBL IDs, as well as the file of the FunCoup database.   #
#  Output is a CSV file with the results of the NEAT enrichment for each one   #
# of the pathways.                                                             #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser",
                         repos = 'http://cran.us.r-project.org')
}
library(argparser)
if(!require("org.Hs.eg.db", quietly = T)){
        BiocManager::install("org.Hs.eg.db", update = F)
}
if(!require("neat", quietly = T)){
        install.packages("neat",
                         repos = 'http://cran.us.r-project.org')
}
library(neat)
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager",
                         repos = 'http://cran.us.r-project.org')
}
if(!require(dplyr, quietly = T)){
        install.packages("dplyr",
                         repos = 'http://cran.us.r-project.org')
}
library(dplyr)
if(!require("GO.db", quietly = T)){
        BiocManager::install("GO.db", update = F)
}


library(GO.db)
# Terminal argument parser
################################################################################
parser <- arg_parser("NEAT functional erichment analyisis.")

parser <- add_argument(parser = parser,
                       arg = c("--DE", "--outName", "--FCFile"),
                       help = c("Differentially expressed genes file",
                                "Name given to the output file",
                                "Path to FunCoup file"),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

inFile <- parsed$DE
outName <- parsed$outName
funCoupFile <- parsed$FCFile

rootDir <- inFile %>% dirname() %>% dirname() %>% dirname() %>% paste0("/")

resuDir <- paste0(rootDir, "results")
resuEnrichDir <- paste0(resuDir, "enrichment/")
resuNEATDir <- paste0(resuEnrichDir, "NEAT/")

if(!dir.exists(resuDir)){
        dir.create(resuDir)
}

if(!dir.exists(resuEnrichDir)){
        dir.create(resuEnrichDir)
}

if(!dir.exists(resuNEATDir)){
        dir.create(resuNEATDir)
}

# Load and parse data
################################################################################

# DE genes file
DE_genes <- read.csv(inFile, row.names = 1)

# funCoup interaction file
funCoup <- read.table(file = funCoupFile,
                      sep = "\t", header = T, comment.char = "")

colnames(funCoup) <- sapply(colnames(funCoup),
                            function(x) strsplit(x,
                                                 ".",
                                                 fixed = T)[[1]][2])

colnames(funCoup)[1] <- "PFC"

# Filter the interactions to keep only high quality ones
funCoup_filt <- funCoup[funCoup$PFC >= 0.75, ]

# Get gene sets for each GO biological process in H. sapiens

# Create cache directory to store filtered GO gene sets
cacheDir <- paste0(rootDir, "cache/")
if(!dir.exists(cacheDir)){
        dir.create(cacheDir)
}
goGeneSets_file <- paste0(cacheDir, "goGeneSets.RData")
if(!file.exists(paste0(goGeneSets_file))){
        goGeneSets <- mapIds(org.Hs.eg.db,
                             keys(org.Hs.eg.db, "GO"),
                             "ENSEMBL",
                             "GO",
                             multiVals = "list")
        # Remove gene GO sets that don't have any gene included in the
        # FunCoup network
        filtGeneSets <- function(geneSets, network){
                pb <- txtProgressBar(min = 0,
                                     max = length(names(geneSets)),
                                     initial = 0,
                                     width = 80,
                                     style = 3)
                sets <- names(geneSets)
                for(i in seq_along(sets)){
                        n <- sets[i]
                        setTxtProgressBar(pb, i)
                        set <- geneSets[[n]]
                        if(sum(set %in% unique(c(network$Gene1,
                                                 network$Gene2))) == 0){
                                geneSets <- geneSets[names(geneSets) != n]
                        }
                }
                close(pb)
                return(geneSets)
        }
        
        print("Filtering gene GO sets to keep just the ones that are
              represented in FunCoup network...")
        goGeneSets_filtered <- filtGeneSets(goGeneSets, funCoup_filt)
        print(sprintf("Done, goGeneSets.RData saved at %s.", cacheDir))
        save(goGeneSets_filtered, file = paste0(cacheDir, "goGeneSets.RData"))
}else{
        print(sprintf("Filtered gene GO sets already exists. Loading
        goGeneSets.RData from %s", cacheDir))
        load(file = paste0(cacheDir, "goGeneSets.RData"))
}


# Create NEAT target genes list, converting symbols to ENSEMBL IDs and
# removing any gene that is not included in the FunCoup interaction file.
neat_input <- rownames(DE_genes)

#neat_input <- mapIds(org.Hs.eg.db,
#                     keys = neat_input,
#                     keytype = "SYMBOL",
#                     column = "ENSEMBL",
#                     fuzzy = T)

neat_input <- neat_input[!is.na(neat_input)]

neat_input <- neat_input[neat_input %in% c(funCoup_filt$Gene1,
                                           funCoup_filt$Gene2)]

neat_input <- list(DE = neat_input)


# Perform NEAT enrichment
################################################################################

print("Running NEAT...")
neat_result <- neat(alist = neat_input,
                    blist = goGeneSets_filtered,
                    network = as.matrix(funCoup_filt[, 3:4]),
                    nettype = "undirected",
                    nodes = sort(unique(c(funCoup_filt$Gene1,
                                          funCoup_filt$Gene2))))
neat_result <- as.data.frame(neat_result)

neat_result$GO_name <- sapply(neat_result$B, Term)
neat_result$ontology <- mapIds(org.Hs.eg.db,
                               neat_result$B,
                               column = "ONTOLOGY",
                               keytype = "GO")

write.csv(neat_result, paste0(resuNEATDir, outName))
print(sprintf("NEAT analysis done. %s saved at %s",
              outName, resuNEATDir))
