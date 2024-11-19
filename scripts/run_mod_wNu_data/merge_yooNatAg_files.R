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

GSE241430_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/GSE241430/GSE241430_raw_counts_GRCh38.p13_NCBI.tsv.gz"

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/"
# Load data
################################################################################
GSE241430 <- readCsvFast(GSE241430_file)

GSE241430_soft <- getGEO(GEO = "GSE241430", GSEMatrix = F)

GSE241430_soft <- parse_soft(GSE241430_soft)

GSE252932_soft <- GSE252932_soft[!grepl("day", GSE252932_soft$sample), ]
GSE252932_soft$pheno <- gsub("-",
                             "",
                             gsub("[0-9]",
                                  "",
                                  GSE252932_soft$sample))
GSE252932_soft$age <- as.numeric(gsub("[^0-9]",
                                      "",
                                      gsub("\\-.*",
                                           "",
                                           GSE252932_soft$sample)))

GSE253174_soft <- GSE253174_soft[grepl("H2O", GSE253174_soft$sample), ]
GSE253174_soft$sample <- gsub(" H2O", "", GSE253174_soft$sample)
GSE253174_soft$pheno <- gsub("-",
                             "",
                             gsub("[0-9]",
                                  "",
                                  GSE253174_soft$sample))
GSE253174_soft$age <- as.numeric(gsub("[^0-9]",
                                      "",
                                      gsub("\\-.*",
                                           "",
                                           GSE253174_soft$sample)))

GSE267613_soft$pheno <- gsub(" ",
                             "",
                             gsub("[0-9]",
                                  "",
                                  GSE267613_soft$sample))
GSE267613_soft$age <- as.numeric(gsub("[^0-9]",
                                      "",
                                      gsub("\\ .*",
                                           "",
                                           GSE267613_soft$sample)))

soft_all <- rbind.data.frame(GSE252932_soft,
                             GSE253174_soft,
                             GSE267613_soft)

dat_all <- cbind.data.frame(GSE252932,
                            GSE253174[, colnames(GSE253174) != "GeneID"],
                            GSE267613[, colnames(GSE267613) != "GeneID"])
dat_all <- dat_all[, colnames(dat_all) %in% c("GeneID", soft_all$GSM)]
writeCsvFst(dat_all, file = sprintf("%sscience_all.csv", outDir))
writeCsvFst(soft_all, file = sprintf("%sscience_soft_all.csv", outDir))
