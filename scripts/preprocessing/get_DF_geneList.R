################################################################################
# Brain clock: obtain file with genes in the integrated dataset.               #
################################################################################

if(!require("argparser")){
        install.packages("argparser",
                         repos = 'https://pbil.univ-lyon1.fr/CRAN/')
}

# Terminal argument parser
################################################################################
parser <- arg_parser("Create gene list file.")

parser <- add_argument(parser = parser,
		       arg = "input",
                       help = "Expression file with genes in columns and samples in rows.",
                       flag = F)

parsed <- parse_args(parser)

# Directory stuff
################################################################################
df_file <- parsed$input

outName <- gsub(".csv", "_genes.csv", df_file)

# Read file and save genes
################################################################################
readCsvFst <- function(path){
	df <- data.frame(data.table::fread(path))
	rownames(df) <- df$V1
	df <- df[, colnames(df) != "V1"]
	return(df)
}

df <- readCsvFst(df_file)

genes <- data.frame(ensembl_gene_id = colnames(df))

write.csv(genes, file = outName)
