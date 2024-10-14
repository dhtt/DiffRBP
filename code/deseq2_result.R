library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

sig_transcripts <- readRDS("data/H1mes_CD48/normalized_count/deseq2/CD48_H1mes.sig_genes.RDS")
refgen <- import.gff("data/refgen/refgen.no_pseudogene.gtf")
sig_genes <- refgen$gene_name[match(sig_transcripts, refgen$transcript_id)]
table(!is.na(sig_genes))

mRNP_chrono_list <- c('all', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII')
RBP_list <- c('RBP2GO', 'INTACT')
network_dirs <- unlist(c(
    lapply(mRNP_chrono_list, function(x) file.path('data/H1mes_CD48/mRNP_chrono', x, 'PPICompare/ResultFiles')),
    lapply(RBP_list, function(x) file.path('data/H1mes_CD48', x, 'PPICompare/ResultFiles'))
    ))

for (network_dir in network_dirs){
    
    protein_attributes <- fread(file.path(network_dir, 'protein_attributes.txt'))
    
    uniprot = protein_attributes$UniProt_ACC[match(sig_genes, protein_attributes$Gene_name)]
    hgnc = protein_attributes$Gene_name
    id_table <- data.table('uniprot' = uniprot, 
                           'hgnc' = hgnc)
    
    diff_ppi <- fread(file.path(network_dir, 'differential_network.txt'))
    print("---------------")
    non_diff_ppi <- diff_ppi %>%
        dplyr::filter(`p-val_adj` <= 0.05) %>%
        dplyr::mutate(
            gene1 = protein_attributes$Gene_name[match(Protein1, protein_attributes$UniProt_ACC)],
            gene2 = protein_attributes$Gene_name[match(Protein2, protein_attributes$UniProt_ACC)]
        )
    print(paste(union(non_diff_ppi$gene1, non_diff_ppi$gene2), collapse = ','))
    
    diff_ppi$gene1 <- id_table$hgnc[match(diff_ppi$Protein1, id_table$uniprot)]
    diff_ppi$gene2 <- id_table$hgnc[match(diff_ppi$Protein2, id_table$uniprot)]
    
    diff_ppi <- diff_ppi %>%
        dplyr::filter((!is.na(gene1) | !is.na(gene2)) & `p-val_adj` <= 0.05)
    print(paste(union(diff_ppi$gene1, diff_ppi$gene2), collapse = ','))
    fwrite(diff_ppi, file.path(network_dir, 'differential_network_deseq2.txt'))
}
