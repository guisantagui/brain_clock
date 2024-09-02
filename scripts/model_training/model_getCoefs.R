################################################################################
# Brain clock: given a model, obtain the non-zero coefficients table and       #  
# save it.                                                                     #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}

library(argparser)
library(h2o)

# Parser
################################################################################
parser <- arg_parser("Given a model, obtain the non-zero coefficients table and
save it.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--outDir"),
                       help = c("RDS file with the model.",
                                "Directory where the resulting CSV file will be stored."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

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

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Directory stuff
################################################################################
modRDS <- parsed$input
outDir <- parsed$outDir

outDir <- addSlashIfNot(outDir)
createIfNot(outDir)


outName <- sprintf("%s%s_coefs.csv", outDir, gsub(".rds", "", basename(modRDS)))

mod <- readRDS(modRDS)

isNonZero <- mod@model$coefficients_table$coefficients != 0

nonZeroCoefs <- mod@model$coefficients_table[isNonZero, ]
nonZeroCoefs <- as.data.frame(nonZeroCoefs)

nonZeroCoefs <- nonZeroCoefs[nonZeroCoefs$names != "Intercept", ]

colnames(nonZeroCoefs) <- gsub("names",
                               "ensembl_gene_id",
                               colnames(nonZeroCoefs))

rownames(nonZeroCoefs) <- nonZeroCoefs$ensembl_gene_id

writeCsvFst(nonZeroCoefs, outName)

print(sprintf("%s saved in %s.", basename(outName), dirname(outName)))