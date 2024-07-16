library(DEXSeq)

# Set the file path
file_path <- "data/Epispliced/dxd.res.RDS"

# Read the file into the variable
normed_count_df <- readRDS(file_path)
head(normed_count_df)