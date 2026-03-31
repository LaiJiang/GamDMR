suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(SOMNiBUS)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 6_test_chunk.R <task_id>")
}
task_id <- as.integer(args[1])
if (is.na(task_id)) {
  stop("task_id must be an integer.")
}

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
data_dir <- file.path(base_dir, "data")
script_dir <- file.path(base_dir, "9_regional")
results_dir <- file.path(base_dir, "results", "somnibus")

sel_cpg_file <- file.path(script_dir, "data", "4_sel_region.txt")
meth_file_loc <- file.path(data_dir, sprintf("chunk_%04d.csv", task_id))
pheno_rdata <- file.path(data_dir, "4_1_data.RData")
region_func_file <- file.path(script_dir, "0_func_region_6.R")
input_func_file <- file.path(script_dir, "0_func_somnibus_input.R")
results_loc <- file.path(results_dir, sprintf("AA_%04d_somnibus.RData", task_id))

if (!file.exists(sel_cpg_file)) stop("Missing selected CpG file: ", sel_cpg_file)
if (!file.exists(meth_file_loc)) stop("Missing methylation chunk file: ", meth_file_loc)
if (!file.exists(pheno_rdata)) stop("Missing phenotype file: ", pheno_rdata)
if (!file.exists(region_func_file)) stop("Missing region helper script: ", region_func_file)
if (!file.exists(input_func_file)) stop("Missing SOMNiBUS input helper script: ", input_func_file)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sel_covariate <- c("AA_only","AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5", "BMI")


sel_cpgs <- read.table(sel_cpg_file, stringsAsFactors = FALSE)
sel_cpgs_chr_start <- sub("-.*$", "", sel_cpgs[[1]])

meth_file <- read.csv(meth_file_loc, header = TRUE, sep = "\t", check.names = FALSE)
load(pheno_rdata)

chunk_cpg_ids <- paste0(meth_file$chr, ":", meth_file$start)
overlap_cpgs <- intersect(sel_cpgs_chr_start, chunk_cpg_ids)

if (length(overlap_cpgs) == 0) {
  message("No selected CpGs found in this chunk.")
  quit(save = "no", status = 0)
}

source(region_func_file, local = TRUE)
source(input_func_file, local = TRUE)

somnibus_input_meth <- convert_to_somnibus_format(meth_file)
rm(meth_file)

pheno_subset <- pheno_file %>%
  select(ID, FID, all_of(sel_covariate))

common_ids <- intersect(somnibus_input_meth$ID, pheno_subset$ID)

somnibus_input_filtered <- somnibus_input_meth %>%
  filter(ID %in% common_ids) %>%
  left_join(pheno_subset, by = "ID") %>%
  select(-FID) %>%
  mutate(
    Meth_Counts = round(Meth_Counts * Total_Counts),
    Meth_Counts = as.numeric(Meth_Counts),
    Total_Counts = as.numeric(Total_Counts)
  )

rm(somnibus_input_meth)

somnibus_input <- somnibus_input_filtered %>%
  filter(Total_Counts != 0) %>%
  na.omit()

rm(somnibus_input_filtered)

if (nrow(somnibus_input) == 0) {
  message("No valid rows remain after filtering.")
  quit(save = "no", status = 0)
}

n_k_dim <- max(5L, as.integer(min(10, length(unique(somnibus_input$Position)) / 20)))

outs <- runSOMNiBUS(
  dat = somnibus_input,
  split = list(approach = "region", gap = 250),
  n.k = rep(n_k_dim, length(sel_covariate) + 1),
  p0 = 0.003,
  p1 = 0.9,
  min.cpgs = 50,
  max.cpgs = 2000,
  verbose = TRUE
)

save(outs, file = results_loc)
message("Saved SOMNiBUS results to: ", results_loc)
