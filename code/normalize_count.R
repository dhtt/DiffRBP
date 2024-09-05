library(DESeq2)

start_time <- Sys.time()
library("stringr", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("DEXSeq", quietly = TRUE)
library("optparse", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("ggpubr", quietly = TRUE)
library("pheatmap", quietly = TRUE)
library("BiocParallel", quietly = TRUE) 

################ PREPARATION ################ 
# Parse arguments for input/output settings
option_list <- list(
    make_option(c("-f", "--count_folder"),
                type = "character",
                help = "path to folder of mRNA-seq counts",
                metavar = "character", 
                default = "data/Epispliced/htseq_count"
    ),
    make_option(c("-a", "--epigenome1"),
                type = "character",
                help = "ID of first epigenome", 
                metavar = "character",
                default = "adipose"
    ),
    make_option(c("-b", "--epigenome2"),
                type = "character",
                help = "ID of second epigenome", 
                metavar = "character",
                default = "aorta"
    ),
    make_option(c("-n", "--num_cores"),
                type = "integer", 
                default = 1,
                help = "number of processing cores", 
                metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
epi1_epi2 <- paste(epi_id1, '_', epi_id2, sep='')
pair <- paste(paste("^", epi_id1, ".*.txt$", sep = ""), paste("^", epi_id2, ".*.txt$", sep = ""), sep = "|")
count_files <- list.files(opt$count_folder, pattern = pair, full.names = TRUE)
file_names <- as.data.table(str_split_fixed(basename(count_files), "\\_|\\.", 3))
result_dir <- strsplit(opt$count_folder, "/")[[1]]
result_dir <- paste(c(result_dir[1:(length(result_dir) - 1)], 'normalized_count'), collapse = "/")
print(result_dir)

# Write logs
log_name <- file(paste('code/logs', epi1_epi2, "log", sep='.'), sep='/'), open = "wt")
sink(log_name, type = c("output", "message"))

cat(paste("---> Working folder: ", opt$count_folder, sep=''), append = TRUE)
cat("\n---> Count files: ", append = TRUE)
cat(basename(count_files), append = TRUE)
cat(paste("\n---> Reference genome: ", gtf_files, sep=''), append = TRUE)


################ DESEQ2 ANALYSIS ################ 
counts <- lapply(count_files, function(x) data.frame(fread(x, sep='\t', header = TRUE))[[2]])
count_data <- do.call(cbind, counts)
colnames(count_data) = paste(file_names$V1, file_names$V2, sep="_")
rownames(count_data) = fread(count_files[1], sep='\t', header = TRUE)[[1]]
count_data = count_data[grep('_', rownames(count_data), invert = T), ]

# Create a metadata table with sample information
metadata <- data.frame(
    tissue = c(file_names$V1), 
    replicate = c(file_names$V2), 
    row.names = colnames(count_data)
    )

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(
    count_data = count_data,
    colData = metadata,
    design = ~ tissue
    )

# Perform normalization and differential expression analysis
dds_normed <- DESeq(dds)
dir.create(file.path(result_dir, 'DESeq2_RDS'), showWarnings = FALSE)
saveRDS(dds_normed, file = file.path(result_dir, 'DESeq2_RDS', epi1_epi2, "RDS", sep='.'))

# Save normalized counts
dds_result = results(dds_normed)
normed_data = DESeq2::counts(dds_normed, normalized = TRUE)
vsd = vst(dds)
dir.create(file.path(result_dir, epi1_epi2), showWarnings = FALSE)

for (i in 1:ncol(normed_data)){
    fwrite((data.frame(normed_data[, i])), 
           row.names = T, sep = '\t', col.names = F,
           file = file.path(result_dir, epi1_epi2, paste(colnames(normed_data)[i], ".txt")))
}
dds_result %>% 
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & !is.na(padj)) %>%
    dplyr::arrange(log2FoldChange)

dds_result[grep('ENST00000264705', rownames(dds)), ]

################ RESULTS ANALYSIS ################ 
# PCA plot for transformed counts 
png(file.path(result_dir, paste(epi_id1, '_', epi_id2, '_PCA.png', sep='')), width = 6, height = 4, unit = 'in', res = 200)
plotPCA(
    vsd, intgroup="tissue", ntop=500) +
    theme_bw() +
    geom_point(size = 4, alpha = 0.8) +
    ggtitle(label="Principal Component Analysis (PCA)", subtitle="Top 500 most variable genes")
dev.off()

#### MAplot for pair comparison with transformed counts ####
make_MAplot <- function(deseq2_result){
    plot = ggmaplot(deseq2_result, 
                    fdr = 0.05, fc = 1, size = 0.5, top = 0,
                    xlab = expression(bold("Log"["2"] ~ "mean expression")),
                    ylab = expression(bold("Log"["2"] ~ "fold change")),
                    palette = c("#64BF52", "#EA6B66", "lightgray"), alpha = 0.6,
                    legend = "top", font.legend = c("bold", 12), font.main = "bold", ggtheme = theme_light()) + 
        theme(
            panel.background = element_rect(colour = "grey", linewidth = 1),
            legend.position = c(0.75, 0.85),
            legend.background=element_blank(),
            legend.text = element_text(size=10),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            panel.grid = element_blank()
        ) + 
        guides(colour = guide_legend(override.aes = list(size=3, alpha=0.8)))
    return(plot)
}

png(file.path(result_dir, paste(epi_id1, '_', epi_id2, '_MA.png', sep='')), width = 6, height = 4, unit = 'in', res = 200)
make_MAplot(dds_result)
dev.off()

#### Heatmap for pair comparison with transformed counts ####
make_heatmap <- function(deseq2_result){
    comparison_type = paste(epi_id1, epi_id2, sep = '_')
    comparison_pair = strsplit(comparison_type, "_")[[1]]
    res_lfc <- subset(deseq2_result, padj < 0.05 & abs(log2FoldChange) > 1)
    
    genes <- rownames(res_lfc)[order(res_lfc$padj, decreasing=TRUE)] 
    
    vst_sig <- vsd[rownames(vsd) %in% genes, grep(paste(comparison_pair, collapse="|"), colnames(vsd))]
    heat <- t(scale(t(assay(vst_sig))))
    colnames(heat) = gsub("_", " ", colnames(heat))
    
    # png(paste("analyzed_result/heatmap_deseq2_toppadj_", comparison_type, ".png", sep=''), width = 7, height = 8, unit = 'in', res = 200)
    plot = pheatmap(
        heat, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE
        )
    return(plot)
}

png(file.path(result_dir, paste(epi_id1, '_', epi_id2, '_heatmap.png', sep='')), width = 6, height = 4, unit = 'in', res = 200)
make_heatmap(dds_result)
dev.off()


################ FINISH ################ 
cat("\n===> FINISHED!", append = TRUE)
end_time <- Sys.time()
cat(paste("\nTotal time:", end_time - start_time, sep = ' '), append = TRUE)

sink()