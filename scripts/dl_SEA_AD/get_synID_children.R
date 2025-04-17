# Given a synapse ID, generates a CSV file with all the children files and their respective
# synIDs. Intended to use as way of passing it to download_synapse.py.
if (!require("synapser", quietly = T)){
        install.packages("synapser",
                         repos=c("http://ran.synapse.org",
                                 "https://cloud.r-project.org"))
}
if (!require("plotUtils", quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
library(synapser)
library(plotUtils)

parent_synID <- "syn51123517"
out_dir <- "/work/projects/age_sorter/data/synapse_data/"
create_dir_if_not(out_dir)

files <- synGetChildren(parent_synID)
files <- as.list(files)
files_df <- data.frame(file_name = unlist(lapply(files, function(x) x$name)),
                       synID = unlist(lapply(files, function(x) x$id)))

write_table_fast(files_df,
                 f = sprintf("%s%s_children.csv",
                             out_dir,
                             parent_synID))