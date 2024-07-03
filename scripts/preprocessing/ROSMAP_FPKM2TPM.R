library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Transforms ROSMAP data (in FPKMs) to TPM")

parser <- add_argument(parser = parser,
                       arg = "input",
                       help = "ROSMAP dataset in FPKM",
                       flag = F)

parsed <- parse_args(parser)

# Directory stuff
################################################################################
ROSMAP_file <- parsed$input
outName <- gsub("FPKM", "TPM", ROSMAP_file)


# Load data
################################################################################

ROSMAP <- read.table(ROSMAP_file, sep = "\t", header = T)

FPKM2TPM <- function(FPKM){
        TPM <- sapply(FPKM, function(x) x/sum(FPKM) * 10^6)
        return(TPM)
}

print("Transforming ROSMAP data from FPKM to TPM...")
ROSMAP_TPM <- ROSMAP[, 1:2]
for(i in 3:ncol(ROSMAP)){
        samp <- colnames(ROSMAP)[i]
        FPKM <- ROSMAP[, i]
        TPM <- data.frame(matrix(FPKM2TPM(FPKM), ncol = 1))
        colnames(TPM) <- samp
        ROSMAP_TPM <- cbind.data.frame(ROSMAP_TPM, TPM)
}

write.table(ROSMAP_TPM, file = outName, row.names=FALSE, sep="\t")
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))