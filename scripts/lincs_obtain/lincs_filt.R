################################################################################
# Brain clock: read metadata of LINCS gctx files for selecting only a given    #
# cell type (in brain clock case, NPCs) and be able to extract only them from  #
# these huge files.                                                            #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("cmapR", quietly = T)) BiocManager::install("cmapR", update = F)
library(cmapR)

# Terminal argument parser
################################################################################
parser <- arg_parser("Filter GCTX LINCS file for desired cell types and write in new files.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--cellType",
                               "--batchSize",
                               "--outDir"),
                       help = c("Input GCTX file.",
                                "Desired cell type.",
                                "Batch size to divide result object in several files.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
#gctxFile <- "/mnt/lscratch/users/gsantamaria/test_large_files/cp_predicted_RNAseq_profiles.gctx"
#cellType <- "NPC"
#sizeBatch <- 4000
#outDir <- "/mnt/lscratch/users/gsantamaria/test_large_files/"

gctxFile <- parsed$input
cellType <- parsed$cellType
sizeBatch <- as.numeric(parsed$batchSize)
outDir <- parsed$outDir

outName <- gsub(".gctx", "", basename(gctxFile), fixed = T)
outName <- sprintf("%s%s_%s.gctx", outDir, outName, cellType)

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Read metadata from gctx
################################################################################
print(sprintf("Loading metadata from %s...", basename(gctxFile)))
metDat <- read_gctx_meta(gctxFile, dim = "column")

# Filter desired celltype and write new gctx files with only desired cell type
################################################################################

# Get indexes for the desired cell type.
cellIdxs <- which(metDat$cell == cellType)

# Load only desired samples based on indexes
nBatches <- ceiling(length(cellIdxs)/sizeBatch)

print(sprintf("Loading %s samples from %s...",
              cellType,
              basename(gctxFile)))
print(sprintf("As n = %s, samples will be split in %s files of %s samples each.",
              as.character(length(cellIdxs)),
              as.character(nBatches),
              as.character(sizeBatch)))
for(i in 1:nBatches){
        idx1 <- (i - 1) * sizeBatch + 1
        idx2 <- idx1 + sizeBatch - 1
        if(idx2 > length(cellIdxs)){
                idx2 <- length(cellIdxs)
        }
        gctx_cellType <- parse_gctx(gctxFile, cid = cellIdxs[idx1:idx2])
        #if(nBatches == 1){
        #        outNameFile <- outName
        #}else{
                outNameFile <- gsub(".gctx", "", outName, fixed = T)
                outNameFile <- sprintf("%s_%s.gctx", outNameFile,
                                       as.character(i))
        #}
        # Write the result in a new gctx file
        write_gctx(gctx_cellType, outNameFile, appenddim = F)
        print(sprintf("%s saved at %s.", basename(outNameFile), dirname(outNameFile)))
}