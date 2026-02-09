#prepare table for beluga running

gene_table <- read.csv(file=paste0(PATH_scr10,"dat/asthma_gene_ranking_with_gene_ranges.csv"), header=TRUE)

gene_table$chr <- c(10,17,9,2,5,5,5,5,6,5,1)

#remove the row with gene_table$GWAS.SNP==""
gene_table <- gene_table[gene_table$GWAS.SNP != "", ]

#remove the row with gene_table$TWAS.eQTL.Evidence=="❌", ]
gene_table <- gene_table[gene_table$TWAS.eQTL.Evidence != "❌", ]

#keep only the columns: Full_Gene_Range_hg19 chr Gene.s.  GWAS.p.value
gene_table <- gene_table[, c("Full_Gene_Range_hg19", "chr", "Gene.s.", "GWAS.p.value", "Drug.Target", "GWAS.SNP")]


gene_table$StrongEvidence <-  "positive"


neg_gene_table <- read.csv(file=paste0(PATH_scr10,"dat/candidate_genes_summary.csv"), header=TRUE)

all_genes <- rbind(gene_table, neg_gene_table)


write.csv(all_genes, file=paste0(PATH_scr10,"dat/asthma_gene_table_test.csv"), row.names = FALSE)