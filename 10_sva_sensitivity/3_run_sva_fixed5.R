#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sva)
  library(data.table)
})
options(mc.cores = 1L)
data.table::setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript 3_run_sva_fixed5.R <1|2|3|10000|25000|50000>")

x <- as.integer(args[1])
vfilters <- c(10000L, 25000L, 50000L)
if (x %in% 1:3) {
  vfilter <- vfilters[x]
} else if (x %in% vfilters) {
  vfilter <- x
} else {
  stop("Argument must be 1,2,3,10000,25000,50000")
}

FIXED_N_SV <- 5L
PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", unset="~/scratch/UQAC/meth"))
out_dir <- path.expand(Sys.getenv(
  "SVA_SENS_OUT",
  unset=file.path(PATH_wk,"results","15_revision","7_sva_threshold_sensitivity")
))
input_rds <- path.expand(Sys.getenv(
  "SVA_PREP_RDS",
  unset=file.path(out_dir,"sva_sensitivity_input.rds")
))
fixed_dir <- file.path(out_dir, "fixed_nsv5")
dir.create(fixed_dir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(input_rds)) stop("Missing prepared input: ", input_rds)

obj <- readRDS(input_rds)
beta_sub <- as.matrix(obj$beta_sub)
pheno_sub <- as.data.frame(obj$pheno_sub)
storage.mode(beta_sub) <- "double"

if (nrow(beta_sub) < vfilter) stop("Not enough CpGs for vfilter=", vfilter)
if (ncol(beta_sub) != nrow(pheno_sub)) stop("Subject count mismatch")
if (!all(is.finite(beta_sub))) stop("beta_sub contains NA/Inf; use prepared complete matrix")

if (!"Non_smoker" %in% names(pheno_sub)) {
  cand <- c("Non-smoker","Non.smoker","Non_smoker")
  smoke_col <- cand[cand %in% names(pheno_sub)][1]
  if (is.na(smoke_col)) stop("Smoking variable not found")
  pheno_sub$Non_smoker <- pheno_sub[[smoke_col]]
}

req <- c("AgeCalc","Sex","Non_smoker","EOSINOpc","LYMPHOpc","MONOpc","NEUTROpc","BMI")
if (length(setdiff(req, names(pheno_sub))) > 0L) stop("Missing SVA covariates")
if (any(!complete.cases(pheno_sub[,req,drop=FALSE]))) stop("Missing covariate values")

if ("ID" %in% names(pheno_sub) && !is.null(colnames(beta_sub))) {
  if (!identical(as.character(pheno_sub$ID), colnames(beta_sub))) stop("Subject order mismatch")
}

mod <- model.matrix(
  ~ AgeCalc + Sex + Non_smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + BMI,
  data=pheno_sub
)
mod0 <- model.matrix(~1, data=pheno_sub)

SEED <- as.integer(Sys.getenv("SVA_SEED", unset="20260816"))
set.seed(SEED)

cat("Running fixed n.sv=5 with vfilter=", vfilter, "\n", sep="")
svobj <- sva(beta_sub, mod, mod0, n.sv=FIXED_N_SV, vfilter=vfilter)
sv <- as.matrix(svobj$sv)
if (ncol(sv) != 5L) stop("sva() did not return 5 SVs")
colnames(sv) <- paste0("SV",1:5)

ids <- if ("ID" %in% names(pheno_sub)) as.character(pheno_sub$ID) else as.character(seq_len(nrow(pheno_sub)))
rownames(sv) <- ids

res <- list(
  analysis="fixed_nsv5",
  vfilter=vfilter,
  fixed_n_sv=5L,
  sv=sv,
  svobj=svobj,
  mod=mod,
  mod0=mod0,
  subject_ids=ids,
  seed=SEED,
  sva_version=as.character(packageVersion("sva")),
  R_version=R.version.string,
  input_rds=input_rds
)

saveRDS(res, file.path(fixed_dir, sprintf("sva_fixed5_vfilter_%05d.rds", vfilter)), compress=FALSE)

sv_dt <- as.data.table(sv)
sv_dt[, ID := ids]
setcolorder(sv_dt, c("ID",paste0("SV",1:5)))
fwrite(sv_dt, file.path(fixed_dir, sprintf("SV_fixed5_vfilter_%05d.csv", vfilter)))

fwrite(
  data.table(
    vfilter=vfilter,
    fixed_n_sv=5L,
    n_cpg_available=nrow(beta_sub),
    n_subjects=ncol(beta_sub),
    seed=SEED,
    sva_version=as.character(packageVersion("sva"))
  ),
  file.path(fixed_dir, sprintf("run_summary_fixed5_vfilter_%05d.csv", vfilter))
)

cat("Completed vfilter=", vfilter, "\n", sep="")
