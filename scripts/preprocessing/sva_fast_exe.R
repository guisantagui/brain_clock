################################################################################
#                                                                              #
# Script for running SVA and detemining the surrogate variables                #
#                                                                              #
################################################################################

if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
if(!require("SmartSVA", quietly = T)){
        install.packages("SmartSVA", repos='http://cran.us.r-project.org')
}
if(!require(sva, quietly = T)) BiocManager::install("sva", update = F)
library(argparser)
library(limma)
library(SmartSVA)
library(sva)



# Terminal argument parser
################################################################################
parser <- arg_parser("Run SVA.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--mod",
                               "--mod0",
                               "--nSV_method",
                               "--saveSVrem",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features in columns, samples rows. CSV.",
                                "Full model RDS object.",
                                "Null model RDS object.",
                                "Method used to compute the number of surrogate variables (leek or RMT)",
                                "If a csv file with the SVs regressed out should be created",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, T, F))

parsed <- parse_args(parser)

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

# Directory stuff
################################################################################
datFile <- parsed$input
mod <- parsed$mod
mod0 <- parsed$mod0
nSV_method <- parsed$nSV_method
saveSVrem <- parsed$saveSVrem
outDir <- parsed$outDir

outBaseName <- gsub(".csv", "", basename(datFile))
suffixMod <- strsplit(gsub(".rds", "", basename(mod)), split = "_")[[1]]
suffixMod <- suffixMod[length(suffixMod)]

# Load data
################################################################################

mod <- readRDS(mod)
mod0 <- readRDS(mod0)
dat <- readCsvFast(datFile)

# Run SVA
################################################################################
print("Running SVA...")
# Filter edata to keep only what's in mod (in case it was subsampled)
dat_rnames <- rownames(dat)
dat_notSubs <- dat[!make.names(rownames(dat)) %in% make.names(rownames(mod)), ]
dat <- dat[match(make.names(rownames(mod)),
                 make.names(rownames(dat))), ]
edata <- t(dat)

print("Computing the number of SVs...")
if(nSV_method == "RMT"){
        df <- data.frame(matrix(mod[, colnames(mod) != "(Intercept)"],
                        nrow = nrow(mod),
                        ncol = ncol(mod) - 1,
                        dimnames = list(rownames(mod),
                                        colnames(mod)[colnames(mod) != "(Intercept)"])))
        ## Determine the number of SVs
        df_vars <- paste(colnames(df), collapse = " + ")
        formula_str <- sprintf("t(edata) ~ %s", df_vars)

        Y.r <- t(resid(lm(as.formula(formula_str), data=df)))
        n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
}else if(nSV_method == "leek"){
        n.sv = num.sv(edata, mod, method="leek")
}


print(sprintf("Number of surrogate variables detected in %s: %s.",
              basename(datFile),
              as.character(n.sv)))

print("Computing the SVs...")
svobj = smartsva.cpp(edata,
                     mod,
                     mod0,
                     n.sv=n.sv)

rownames(svobj$sv) <- rownames(dat)

outName <- sprintf("%s%s_%s_svobj.rds",
                   outDir,
                   outBaseName,
                   suffixMod)

saveRDS(svobj, file = outName)
print(sprintf("%s saved at %s",
              basename(outName),
              dirname(outName)))

# Remove effect of the surrogate variables from the input matrix and save it
################################################################################

if (nrow(dat_notSubs) > 0){
        print("Fitting linear model to infer SVs in samples that were not subsampled...")
        sv_mat <- svobj$sv
        colnames(sv_mat) <- paste("SV", 1:ncol(sv_mat), sep = "_")
        surrogate_model <- lm(sv_mat ~ ., data = as.data.frame(dat))
        saveRDS(surrogate_model, file = sprintf("%ssv_lmod.rds", outDir))
        print(sprintf("sv_lmod.rds saved at %s", outDir))
        print("Inferring SVs in not-subsampled samples...")
        sv_rest <- predict(surrogate_model, dat_notSubs)
        if (is.vector(sv_rest)){
                sv_rest <- data.frame(SV_1 = sv_rest)
        }
        # Merge SV matrices and data matrices and reorder row names
        sv_all <- rbind.data.frame(sv_mat, sv_rest)
        sv_all <- sv_all[match(make.names(dat_rnames), make.names(rownames(sv_all))), , drop = FALSE]
        dat <- rbind.data.frame(dat, dat_notSubs)
        dat <- dat[match(make.names(dat_rnames), make.names(rownames(dat))), ]
}else{
        sv_all <- svobj$sv
}

if(saveSVrem){
        # Create new model matrix with only age, as it is the only covariate we want to
        # preserve.
        #mod4adj <- model.matrix(~age_death,
        #                        data = data.frame(age_death = mod[, "age_death"]))

        edata_adj <- removeBatchEffect(t(dat), covariates = sv_all)

        edata_adj <- t(edata_adj)

        outNameAdj <- sprintf("%s%s_%s_svaAdj.csv",
                              outDir,
                              outBaseName,
                              suffixMod)

        writeCsvFst(edata_adj, file = outNameAdj)
        print(sprintf("%s saved in %s",
                      basename(outNameAdj),
                      dirname(outNameAdj)))
}
