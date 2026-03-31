suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
paper_results <- file.path(base_dir, "scr", "14_paper", "results")
rerun_results <- file.path(base_dir, "scr", "8_rerun", "results")

dir.create(paper_results, recursive = TRUE, showWarnings = FALSE)

m1_info <- fread(file.path(rerun_results, "7_gsea_M1_genes.csv"))
m2_info <- fread(file.path(rerun_results, "7_gsea_M2_genes.csv"))
m3_info <- fread(file.path(rerun_results, "7_gsea_M3_genes.csv"))

load(file.path(rerun_results, "4_eval_M123.RData"))

m1_cpgs <- M1_manhattan[M1_manhattan$pval < 1e-5, ]
m3_cpgs <- M3_manhattan[M3_manhattan$pval < 1e-5, ]

table_m1 <- merge(
  m1_cpgs[, c("CpG", "pval", "coef")],
  m1_info[, c("CpG", "genes_entrez", "genes_symbol")],
  by = "CpG",
  all.x = TRUE
)
table_m2 <- merge(
  setNames(M2_cpgs[, c("CpG", "M2_pval", "M2_coef")], c("CpG", "pval", "coef")),
  m2_info[, c("CpG", "genes_entrez", "genes_symbol")],
  by = "CpG",
  all.x = TRUE
)
table_m3 <- merge(
  m3_cpgs[, c("CpG", "pval", "coef")],
  m3_info[, c("CpG", "genes_entrez", "genes_symbol")],
  by = "CpG",
  all.x = TRUE
)

table_m1$model <- "M1"
table_m2$model <- "M2"
table_m3$model <- "M3"
final_table <- rbind(table_m1, table_m2, table_m3, fill = TRUE)

final_table$overlap_gene <- NA_character_
for (i in seq_len(nrow(final_table))) {
  cpg <- final_table$CpG[i]
  gene <- final_table$genes_symbol[i]
  models_for_cpg <- unique(final_table$model[final_table$CpG == cpg])
  models_for_gene <- unique(final_table$model[final_table$genes_symbol == gene])
  overlap_models <- unique(c(models_for_cpg, models_for_gene))
  if (length(overlap_models) > 1) {
    final_table$overlap_gene[i] <- paste(overlap_models, collapse = ";")
  }
}
final_table$overlap_gene <- gsub(";NA", "", final_table$overlap_gene)

fwrite(final_table, file = file.path(paper_results, "9_table_significant_CpGs_M1_M2_M3.csv"))
