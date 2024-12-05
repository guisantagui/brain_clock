################################################################################
# Brain clock: NEAT Enrichment. The input must be a dataframe (.csv file)      #
# where the rownames are ENSEMBL IDs, as well as the file of the FunCoup       #
# database. Output is a CSV file with the results of the NEAT enrichment for   #
# each one of the pathways.                                                    #
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
                       arg = c("--DE",
                               "--outName",
                               "--extGeneSets",
                               "--FCFile",
                               "--FCVersion",
                               "--outDir"),
                       help = c("Differentially expressed genes file",
                                "Name given to the output file",
                                "File to external gene sets object (RDS file consisting in a list, with the names being the names of the gene sets, and the elements the ENSEMBL IDs).",
                                "Path to FunCoup file",
                                "Version of the FunCoup file (5 or 6)",
                                "Output directory to save the results."),
                       flag = c(F, F, F, F, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

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

# Directory stuff
################################################################################

inFile <- parsed$DE
outName <- parsed$outName
extGeneSets <- parsed$extGeneSets
funCoupFile <- parsed$FCFile
FCVersion <- as.numeric(parsed$FCVersion)
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

if (extGeneSets == "NULL"){
        extGeneSets <- NULL
}

outName <- paste0(outDir, outName)

# Load and parse data
################################################################################

# DE genes file
DE_genes <- read.csv(inFile, row.names = 1)

# funCoup interaction file
funCoup <- read.table(file = funCoupFile,
                      sep = "\t", header = T, comment.char = "")

colnames(funCoup) <- gsub("X.0", "X0", colnames(funCoup))

colnames(funCoup) <- sapply(colnames(funCoup),
                            function(x) strsplit(x,
                                                 ".",
                                                 fixed = T)[[1]][2])

# Filter the interactions to keep only high quality ones
if(FCVersion == 5){
        funCoup_filt <- funCoup[funCoup$PFC >= 0.75, ]
}else if (FCVersion == 6){
        funCoup_filt <- funCoup[funCoup$PPV >= 0.75, ]
}

# Get gene sets for each GO biological process in H. sapiens

# Create cache directory to store filtered GO gene sets
if(is.null(extGeneSets)){
        cacheDir <- paste0(outDir, "cache/")
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
}else{
        print(sprintf("Loading external gene sets from %s, at %s.",
                      basename(extGeneSets),
                      dirname(extGeneSets)))
        goGeneSets <- readRDS(extGeneSets)
        print("Filtering gene sets to keep just the ones that are
                      represented in FunCoup network...")
        goGeneSets_filtered <- filtGeneSets(goGeneSets, funCoup_filt)
}

# Create NEAT target genes list, converting symbols to ENSEMBL IDs and
# removing any gene that is not included in the FunCoup interaction file.
neat_input <- rownames(DE_genes)

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

if (is.null(extGeneSets)){
        neat_result$GO_name <- sapply(neat_result$B, Term)
        neat_result$ontology <- mapIds(org.Hs.eg.db,
                                       neat_result$B,
                                       column = "ONTOLOGY",
                                       keytype = "GO")
}
write.csv(neat_result, outName)
print(sprintf("NEAT analysis done. %s saved at %s",
              basename(outName), dirname(outName)))
