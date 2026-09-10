#!/usr/bin/env Rscript

# =============================================================================
# 1_run_sva_threshold_sensitivity.R
#
# Run one SVA threshold:
#   task 1 -> vfilter = 10,000  (primary setting)
#   task 2 -> vfilter = 25,000
#   task 3 -> vfilter = 50,000
#
# The model specification is intentionally kept identical to the user's
# original SVA code. In particular, AA_only is NOT added here, because this
# sensitivity analysis is intended to isolate the effect of vfilter.
#
# Original model:
#   mod  ~ AgeCalc + Sex + Non-smoker + EOSINOpc + LYMPHOpc +
#          MONOpc + NEUTROpc + BMI
#   mod0 ~ 1
#
# SVA calls:
#   n.sv <- num.sv(beta_sub, mod, vfilter=..., method="be")
#   sva(beta_sub, mod, mod0, n.sv=n.sv, vfilter=...)
# =============================================================================

suppressPackageStartupMessages({
  library(sva)
  library(data.table)
})

options(mc.cores = 1L)
data.table::setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript 1_run_sva_threshold_sensitivity.R <task_id_or_vfilter>")
}

x <- as.integer(args[1])
if (is.na(x)) stop("Argument must be an integer.")

vfilters <- c(10000L, 25000L, 50000L)

if (x %in% seq_along(vfilters)) {
  vfilter <- vfilters[x]
} else if (x %in% vfilters) {
  vfilter <- x
} else {
  stop("Argument must be 1, 2, 3, 10000, 25000, or 50000.")
}

PATH_wk <- path.expand(
  Sys.getenv("METH_BASE_DIR", unset = "~/scratch/UQAC/meth")
)

out_dir <- path.expand(
  Sys.getenv(
    "SVA_SENS_OUT",
    unset = file.path(
      PATH_wk, "results", "15_revision", "7_sva_threshold_sensitivity"
    )
  )
)

input_rds <- file.path(out_dir, "sva_sensitivity_input.rds")
if (!file.exists(input_rds)) {
  stop("Missing prepared input: ", input_rds)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading:", input_rds, "\n")
obj <- readRDS(input_rds)

beta_sub <- obj$beta_sub
pheno_sub <- obj$pheno_sub

cat("beta_sub dimensions:", nrow(beta_sub), "CpGs x", ncol(beta_sub), "subjects\n")
cat("Running vfilter =", vfilter, "\n")

if (nrow(beta_sub) < vfilter) {
  stop("beta_sub has fewer rows than requested vfilter.")
}

if (ncol(beta_sub) != nrow(pheno_sub)) {
  stop("beta_sub columns and pheno_sub rows do not match.")
}

if (!identical(colnames(beta_sub), as.character(pheno_sub$ID))) {
  stop("Subject order mismatch between beta_sub and pheno_sub.")
}

# ------------------------- Original SVA model -------------------------------

mod <- model.matrix(
  ~ AgeCalc + Sex + Non_smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + BMI,
  data = pheno_sub
)

mod0 <- model.matrix(~ 1, data = pheno_sub)

if (qr(mod)$rank < ncol(mod)) {
  stop("Full SVA model matrix is rank deficient.")
}

# Reproducibility for the Buja-Eyuboglu permutation procedure.
# Use the SAME seed at every threshold so threshold is the intended difference.
SVA_SEED <- as.integer(Sys.getenv("SVA_SEED", unset = "20260816"))
set.seed(SVA_SEED)

cat("Estimating n.sv with num.sv(method='be', vfilter=", vfilter, ")...\n", sep = "")

n_sv <- num.sv(
  beta_sub,
  mod,
  vfilter = vfilter,
  method = "be"
)

cat("Estimated n.sv:", n_sv, "\n")

set.seed(SVA_SEED)

cat("Running sva()...\n")

svobj <- sva(
  beta_sub,
  mod,
  mod0,
  n.sv = n_sv,
  vfilter = vfilter
)

sv <- svobj$sv

if (is.null(dim(sv))) {
  sv <- matrix(sv, ncol = 1L)
}

colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
rownames(sv) <- pheno_sub$ID

# ------------------------------ Save ----------------------------------------

result <- list(
  vfilter = vfilter,
  n_sv = n_sv,
  sv = sv,
  svobj = svobj,
  mod = mod,
  mod0 = mod0,
  subject_ids = pheno_sub$ID,
  seed = SVA_SEED,
  sva_version = as.character(packageVersion("sva")),
  R_version = R.version.string,
  input_rds = input_rds
)

rds_out <- file.path(
  out_dir,
  sprintf("sva_vfilter_%05d.rds", vfilter)
)

saveRDS(result, rds_out, compress = FALSE)

sv_dt <- as.data.table(sv)
sv_dt[, ID := pheno_sub$ID]
setcolorder(sv_dt, c("ID", setdiff(names(sv_dt), "ID")))

csv_out <- file.path(
  out_dir,
  sprintf("SV_vfilter_%05d.csv", vfilter)
)

fwrite(sv_dt, csv_out)

summary_out <- data.table(
  vfilter = vfilter,
  n_cpg_available = nrow(beta_sub),
  n_subjects = ncol(beta_sub),
  n_sv = n_sv,
  seed = SVA_SEED,
  sva_version = as.character(packageVersion("sva"))
)

fwrite(
  summary_out,
  file.path(out_dir, sprintf("summary_vfilter_%05d.csv", vfilter))
)

cat("\nCompleted vfilter =", vfilter, "\n")
cat("n.sv =", n_sv, "\n")
cat("Saved:", rds_out, "\n")
cat("Saved:", csv_out, "\n")
