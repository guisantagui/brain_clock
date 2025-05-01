################################################################################
# Script for running bigPCA                                                    #
################################################################################

if(!require(bigstatsr, quietly = T)){
        install.packages("bigstatsr",
                         repos='http://cran.us.r-project.org')
}
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(bigstatsr)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Run PCA for big matrices.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--nPCs",
                               "--stand",
                               "--outDir"),
                       help = c("Input matrix, features in columns, samples rows. CSV.",
                                "Number of principal components.",
                                "If variables should be standarized before PCA.",
                                "Output directory where placing the results."),
                       flag = c(F, F, T, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
dat <- as.character(parsed$input)
nPCs <- as.numeric(parsed$nPCs)
stnd <- parsed$stand
outDir <- as.character(parsed$outDir)

outBaseName <- gsub(".csv", "", basename(dat))

# Load data
################################################################################
dat <- data.frame(data.table::fread(dat))
rownames(dat) <- dat$V1
dat <- dat[, colnames(dat) != "V1"]

# Function
################################################################################
big_PCA <- function(df, nComps = 20, stand = F){
        if(stand){
                print("Standarizing the input dataset...")
                # Remove genes with SD == 0
                df <- df[, apply(dat, 2, sd) != 0]
                scaled_DF <- apply(df, 2, function(x) (x - mean(x))/sd(x))
                scal = T
                centr = T
        }else{
                scaled_DF <- df
                scal = F
                centr = F
        }
        
        scaled_DF_fbm <- FBM(nrow(scaled_DF),
                             ncol(scaled_DF),
                             init = scaled_DF)
        
        print("Computing SVD...")
        scaled_DF_svd <- big_SVD(scaled_DF_fbm,
                                 fun.scaling = big_scale(center = F,
                                                         scale = F),
                                 k = nComps)
        
        # Multiply U and D to get the scores
        print("Generating PCA output...")
        print("Computing scores matrix...")
        scores <- scaled_DF_svd$u %*% diag(scaled_DF_svd$d)
        rownames(scores) <- rownames(scaled_DF)
        colnames(scores) <- paste0("PC", as.character(1:ncol(scores)))
        rotation <- scaled_DF_svd$v
        rownames(rotation) <- colnames(scaled_DF)
        colnames(rotation) <- paste0("PC", as.character(1:ncol(rotation)))
        print("Computing eigenvalues...")
        eigenVals <- scaled_DF_svd$d^2
        
        print("Computing total variance of the initial dataset...")
        totVar <- sum(diag(t(scaled_DF) %*% scaled_DF))
        varExp <- eigenVals/totVar
        
        cumulative_sum <- function(v, index = 1, current_sum = 0) {
                if (index > length(v)){
                        return(c())
                }
                
                current_sum <- current_sum + v[index]
                
                return(c(current_sum, cumulative_sum(v, index + 1, current_sum)))
        }
        
        print("Computing PC cumulative variance...")
        cumVar <- cumulative_sum(varExp)
        pcSD <- apply(scores, 2, sd)
        
        summMat <- rbind(pcSD, varExp, cumVar)
        rownames(summMat) <- c("Standard deviation",
                               "Proportion of Variance",
                               "Cumulative Proportion")
        outList <- list(sdev = pcSD,
                        rotation = rotation,
                        scale = scal,
                        center = centr,
                        x = scores,
                        summary = summMat)
        print("Done")
        return(outList)
}

# Run big PCA
################################################################################
bigPCAres <- big_PCA(df = dat,
                     nComps = nPCs,
                     stand = stnd)

outName <- sprintf("%s%s_pca.rds", outDir, outBaseName)
saveRDS(bigPCAres, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))