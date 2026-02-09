#tomorrow: collect the csv files to prepare the final result table for report.

#ask GPT to prepare a ipynb report, add some descriptions of the GAM model, put in figures (excluding manhattan plot), and table. 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk,"results/10_sanity/")
PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")
PATH_scr10 <- paste0(PATH_wk, "scr/10_sanity/")



all_genes <- read.csv( file=paste0(PATH_scr10,"dat/asthma_gene_table_test.csv"))
all_genes <- read.csv( file=paste0(PATH_scr10,"dat/asthma_gene_ranking_with_gene_ranges.csv"))


gene_results <- gsub(".csv","", list.files(PATH_results, pattern = "*.csv", full.names = FALSE))

#intersect(gene_results, all_genes$Gene.s.)

#all_genes<- all_genes[-c(2,7),]

gene_list <- unique(all_genes$Gene.s.)

i_gene_list <- gene_list[1]

pval_list <- NULL

threshold <- 1e-5
for(i_gene_list in gene_list){

first_line <- readLines(paste0(PATH_results, i_gene_list, ".csv"), n = 1)

if("\"\""== first_line){
  #if the first line is empty, then skip this gene
  i_pval <- 1
}
else{

i_gene_dmr <- read.csv(paste0(PATH_results,i_gene_list,".csv"), header=TRUE)
  i_pval <- min(i_gene_dmr$pval_smooth_AA)

}

pval_list <- c(pval_list, i_pval)


}

all_genes$DMR_pval <- pval_list

all_genes$DMR_1e5 <- (pval_list < 1e-5)

write.csv(all_genes, file=paste0(PATH_scr10,"dat/5_col_asthma_gene_table.csv"), row.names = FALSE)


#all_genes <- read.csv(file=paste0(PATH_scr10,"dat/5_col_asthma_gene_table.csv"))