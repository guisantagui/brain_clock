################################################################################
# Brain clock: test the perturbations in the pert samples.                     #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}
library(h2o)

# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--modFile",
                               "--metDat",
                               "--ageTransPars",
                               "--sizeBatch",
                               "--mem",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features columns, samples rows. CSV.",
                                "File of the h2o model (generated with mod_train_and_test.R)",
                                "Metadata file (including ages).",
                                "RDS file with Gompertz-Makeham parameters for transformation.",
                                "Size of each batch for predictions.", 
                                "Memory to allocate to h2o instance.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F, F, F))

parsed <- parse_args(parser)


# Directory stuff
################################################################################

datFile <- parsed$input
modFile <- parsed$modFile
metDatFile <- parsed$metDat
ageTransParFile <- parsed$ageTransPars
sizeBatch <- as.numeric(parsed$sizeBatch)
mem <- parsed$mem
outDir <- parsed$outDir


outName <- sprintf("%spred_ages.csv", outDir)

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Functions
################################################################################

# Read csv faster
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

# Functions for doing age transformation
test <- function(x){
        exp(-vals[3]*x-(vals[1]/vals[2])*(exp(vals[2]*x)-1))
}

inverse = function (f, lower = 0, upper = 110) {
        function (y) uniroot((function (x) f(x) - y),
                             lower = lower,
                             upper = upper)[1]
}
agetrans_inverse = function(lower,upper) inverse(test, lower, upper)

# Clips transformed ages to be in a range that can be backtransformed to
# chronological ages, transforms back to chronological age and returns the 
# result.
back2Age <- function(transAge){
        transAge[transAge > 1] <- 1
        transAge[transAge < 0] <- .000012
        transAge2Age <- sapply(transAge,
                               function(x) agetrans_inverse(0, 110)(x))
        transAge2Age <- unlist(transAge2Age)
        names(transAge2Age) <- gsub(".root", "", names(transAge2Age))
        return(transAge2Age)
}

# Predict the age viven the expression matrix and the model
predictAge <- function(model,
                       predExpr
){
        predExpr <- predExpr[model@parameters$x,]
        predAges <- c()
        maxiter <- ceiling(ncol(predExpr)/1000)
        for(i in 1:maxiter){
                #i <- 1
                from <- ((i-1)*1000+1)
                to <- min(i*1000,ncol(predExpr))
                
                inp_h2o <- as.h2o(t(predExpr[,from:to]))
                
                pred <- h2o.predict(model,inp_h2o)
                pred <- as.data.frame(pred)$predict
                names(pred) <- colnames(predExpr)[from:to]
                
                predAges <- c(predAges,pred)
        }
        return(predAges)
}

# Load the data
################################################################################
dat <- readCsvFast(datFile)
metDat <- readCsvFast(metDatFile)

# Load Gompertz-Makeham parameters
vals <- readRDS(ageTransParFile)
vals <- round(vals,6) #Makes the parameters more readable with little loss of accuracy

# Keep only perturbated samples
pertSamps <- metDat$specimenID[metDat$diagn_4BrainClck == "perturbations"]
dat <- dat[make.names(pertSamps), ]

# Initialize h2o and load model
conn <- h2o.init(max_mem_size=mem)

mod <- h2o.loadModel(modFile)

# Do predictions
################################################################################

# Load only desired samples based on indexes
nBatches <- ceiling(nrow(dat)/sizeBatch)

predVec <- c()
for(i in 1:nBatches){
        idx1 <- (i - 1) * sizeBatch + 1
        idx2 <- idx1 + sizeBatch - 1
        if(idx2 > nrow(dat)){
                idx2 <- nrow(dat)
        }
        pred <- predictAge(mod, t(dat[idx1:idx2, ]))
        predVec <- c(predVec, pred)
}

predsDF <- data.frame(specimenID = names(predVec),
                      trans_age = predVec,
                      chron_age = back2Age(predVec))

writeCsvFst(predsDF, file = outName)
print(sprintf("%s saved in %s.", basename(outName), dirname(outName)))
h2o.shutdown(prompt = F)

################################################################################