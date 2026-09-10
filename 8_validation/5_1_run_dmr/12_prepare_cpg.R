#collect tables and fgiures for the supplementary materials of the paper

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk, "/results/")
PATH_scr15 <- paste0(PATH_wk, "/scr/15_revision/")
PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")
PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"

library(data.table)
library(dplyr)
#collect the M1, M2, M3 results after rerun




# Load the data
M1_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M1_genes.csv"))
M2_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M2_genes.csv"))
M3_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M3_genes.csv"))

load(file = paste0(PATH_wk,"/scr/8_rerun/results/4_eval_M123.RData"),verbose=TRUE)

# Extract significant CpGs
M1_cpgs <- M1_manhattan[M1_manhattan$pval < 0.001, ]

M2_cpgs <- M2_infor[M2_infor$pval < 0.001, ]

#we first tried 0.05, then moved to 0.001
M3_cpgs <- M3_manhattan[M3_manhattan$pval < 0.05, ]

# Create table_M1 by merging M1_cpgs with M1_infor
table_M1 <- M1_cpgs[, c("CpG", "pval", "coef")]
table_M1 <- merge(table_M1, 
                  M1_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)

# Reorder columns if desired
table_M1 <- table_M1[, c("CpG", "pval", "coef", "genes_entrez", "genes_symbol")]

# Similarly for Model 2
table_M2 <- M2_cpgs[, c("CpG", "pval")]
#now rename columns
table_M2 <- merge(table_M2, 
                  M2_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)
table_M2 <- table_M2[, c("CpG", "pval",  "genes_entrez", "genes_symbol")]


table_M2$coef <- rep(1,nrow(table_M2))  # Add a 'coef' column with value 1 for all rows


# Similarly for Model 3
table_M3 <- M3_cpgs[, c("CpG", "pval", "coef")]
table_M3 <- merge(table_M3, 
                  M3_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)
table_M3 <- table_M3[, c("CpG", "pval", "coef", "genes_entrez", "genes_symbol")]


#now combine the three tables into a single table, add a column indicating the model
table_M1$model <- "M1"
table_M2$model <- "M2"
table_M3$model <- "M3"
final_table <- rbind(table_M1, table_M2, table_M3)

#write to a csv file
fwrite(final_table, file = paste0(PATH_results, "15_revision/5_cv/12_final_cpgs_table.csv"))




region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)


library(data.table)

setDT(final_table)
setDT(region_file)

# Preserve original row order
final_table[, row_id__ := .I]

# Extract chromosome and CpG start position
# Example: "1:10002719-10002720"
final_table[, cpg_chr := sub(":.*$", "", CpG)]

final_table[, cpg_start := as.integer(
  sub("^[^:]+:([0-9]+)-.*$", "\\1", CpG)
)]

# Make chromosome types consistent
final_table[, cpg_chr := as.character(cpg_chr)]
region_file[, chr := as.character(chr)]

# Attach data_chunk_id based on:
# chromosome match AND region_start <= cpg_start <= region_end
final_table[
  region_file,
  on = .(
    cpg_chr = chr,
    cpg_start >= region_start,
    cpg_start <= region_end
  ),
  data_chunk_id := i.data_chunk_id
]

# Restore original order
setorder(final_table, row_id__)
final_table[, row_id__ := NULL]


library(data.table)

setDT(final_table)
setDT(region_file)

# Make chromosome types consistent
final_table[, cpg_chr := as.character(cpg_chr)]
region_file[, chr := as.character(chr)]

# ------------------------------------------------------------------
# Previous region:
# closest region with region_end < cpg_start, on the same chromosome
# ------------------------------------------------------------------

prev_lookup <- region_file[
  ,
  .(
    chr,
    region_end,
    data_chunk_id_prev = data_chunk_id
  )
]

setkey(prev_lookup, chr, region_end)

final_table[, data_chunk_id_prev :=
  prev_lookup[
    final_table,
    on = .(
      chr = cpg_chr,
      region_end < cpg_start
    ),
    mult = "last",
    data_chunk_id_prev
  ]
]

# ------------------------------------------------------------------
# Next region:
# closest region with region_start > cpg_start, on the same chromosome
# ------------------------------------------------------------------

next_lookup <- region_file[
  ,
  .(
    chr,
    region_start,
    data_chunk_id_next = data_chunk_id
  )
]

setkey(next_lookup, chr, region_start)

final_table[, data_chunk_id_next :=
  next_lookup[
    final_table,
    on = .(
      chr = cpg_chr,
      region_start > cpg_start
    ),
    mult = "first",
    data_chunk_id_next
  ]
]

final_table[
  is.na(data_chunk_id) &
  !is.na(data_chunk_id_prev) &
  !is.na(data_chunk_id_next) &
  data_chunk_id_prev == data_chunk_id_next,
  data_chunk_id := data_chunk_id_prev
]

fwrite(final_table, file = paste0(PATH_results, "15_revision/5_cv/12_final_cpgs_table.csv"))


#now upload this file to server.


scp  /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/12_final_cpgs_table.csv  laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/15_revision/5_cv/


dir.create(
   paste0(PATH_results,
  "/15_revision/5_cv/6_single_cpg/" ),
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  final_table, paste0(PATH_results,
  "/15_revision/5_cv/6_single_cpg/final_table_for_cv.rds" )
)


scp  /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/6_single_cpg/final_table_for_cv.rds  laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/


scp  /mnt/c/Per/LaiJiang/Project/UQAC/meth/scr/15_revision/5_cv/5_1_run_dmr/12_run_single_cpg_training_associations.R  laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/15_revision/5_cv/5_1_run_dmr/12_run_single_cpg_training_associations.R  



scp  /mnt/c/Per/LaiJiang/Project/UQAC/meth/scr/15_revision/5_cv/5_1_run_dmr/12_run_single_cpg_training_associations.sh  laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/15_revision/5_cv/5_1_run_dmr/12_run_single_cpg_training_associations.sh  


library(data.table)

setDT(final_table)

smoke_M123 <- final_table[
  model %in% c("M1", "M2", "M3"),
  .SD[1L],
  by = model
]

smoke_M123


saveRDS(
  smoke_M123, 
  paste0(PATH_results,
  "/15_revision/5_cv/6_single_cpg/final_table_smoke_M123.rds")
)




scp  /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/6_single_cpg/final_table_smoke_M123.rds  laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/
