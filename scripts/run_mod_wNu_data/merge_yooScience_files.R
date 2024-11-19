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
        access <- sft@header$geo_accession
        out_df <- data.frame(GSE = rep(access, length(gsms)),
                             GSM = gsms,
                             sample = samps)
        return(out_df)
}

# Directory stuff
################################################################################

GSE252932_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/GSE252932/GSE252932_raw_counts_GRCh38.p13_NCBI.tsv.gz"
GSE253174_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/GSE253174/GSE253174_raw_counts_GRCh38.p13_NCBI.tsv.gz"
GSE267613_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/GSE267613/GSE267613_raw_counts_GRCh38.p13_NCBI.tsv.gz"

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/"
# Load data
################################################################################
GSE252932 <- readCsvFast(GSE252932_file)
GSE253174 <- readCsvFast(GSE253174_file)
GSE267613 <- readCsvFast(GSE267613_file)

GSE252932_soft <- getGEO(GEO = "GSE252932", GSEMatrix = F)
GSE253174_soft <- getGEO(GEO = "GSE253174", GSEMatrix = F)
GSE267613_soft <- getGEO(GEO = "GSE267613", GSEMatrix = F)

GSE252932_soft <- parse_soft(GSE252932_soft)
GSE253174_soft <- parse_soft(GSE253174_soft)
GSE267613_soft <- parse_soft(GSE267613_soft)

GSE252932_soft <- GSE252932_soft[!grepl("day", GSE252932_soft$sample), ]
GSE252932_soft$treatment <- "H2O"
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

#GSE253174_soft <- GSE253174_soft[grepl("H2O", GSE253174_soft$sample), ]
GSE253174_soft$treatment <- gsub("\\-.*", "", gsub(".* ", "", GSE253174_soft$sample))
GSE253174_soft$sample <- gsub(" H2O| 3TC", "", GSE253174_soft$sample)
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

GSE267613_soft$treatment <- "H2O"
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
