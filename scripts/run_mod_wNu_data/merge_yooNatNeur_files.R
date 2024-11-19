library(GEOquery)

# Functions 
################################################################################

# Read csv faster
readCsvFast <- function(f){
        df <- data.frame(data.table::fread(f))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

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

parse_soft <- function(sft){
        gsms <- c()
        samps <- c()
        for(g in names(sft@gsms)){
                gsms <- c(gsms, g)
                samps <- c(samps, sft@gsms[[g]]@header$title)
        }
        out_df <- data.frame(GSM = gsms, sample = samps)
        return(out_df)
}

# Directory stuff
################################################################################

GSE194241_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/GSE194241/GSE194241_raw_counts_GRCh38.p13_NCBI.tsv.gz"

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/"
# Load data
################################################################################
GSE194241 <- readCsvFast(GSE194241_file)

GSE194241_soft <- getGEO(GEO = "GSE194241", GSEMatrix = F)

GSE194241_soft <- parse_soft(GSE194241_soft)

GSE194241_soft <- GSE252932_soft[!grepl("day", GSE252932_soft$sample), ]
GSE194241_soft$pheno <- gsub("\\_.*",
                             "",
                             GSE194241_soft$sample)

dat_all <- GSE194241
soft_all <- GSE194241_soft
writeCsvFst(dat_all, file = sprintf("%snatNeur_all.csv", outDir))
writeCsvFst(soft_all, file = sprintf("%snatNeur_soft_all.csv", outDir))
