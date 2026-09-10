#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

PATH_WK <- path.expand(Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/"))
assoc_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_ASSOC_OUTPUT",
  file.path(PATH_WK, "results/15_revision/5_cv/6_single_cpg/1_training_associations")
))
collect_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_COLLECT_ROOT",
  file.path(PATH_WK, "results/15_revision/5_cv/6_single_cpg/2_collected_significant")
))
n_assoc_jobs <- as.integer(Sys.getenv("N_ASSOC_JOBS", "990"))
n_collect_jobs <- as.integer(Sys.getenv("N_COLLECT_JOBS", "99"))
candidate_p_max <- as.numeric(Sys.getenv("CANDIDATE_P_MAX", "0.001"))
allow_incomplete <- Sys.getenv("ALLOW_INCOMPLETE", "0") %in% c("1","TRUE","true","T","yes","YES")
overwrite <- Sys.getenv("OVERWRITE_COLLECT", "1") %in% c("1","TRUE","true","T","yes","YES")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1L) stop("Pass collection array task ID.")
collect_id <- as.integer(args[1])
if (is.na(collect_id) || collect_id<1L || collect_id>n_collect_jobs) stop("Invalid collect_id.")
if (!is.finite(candidate_p_max) || candidate_p_max<=0 || candidate_p_max>1) stop("Invalid CANDIDATE_P_MAX.")

result_dir <- file.path(assoc_root,"job_results")
upstream_qc_dir <- file.path(assoc_root,"job_qc")
partial_dir <- file.path(collect_root,"partial")
dir.create(partial_dir, recursive=TRUE, showWarnings=FALSE)

assoc_start <- floor((collect_id-1L)*n_assoc_jobs/n_collect_jobs)+1L
assoc_end <- floor(collect_id*n_assoc_jobs/n_collect_jobs)
assoc_job_ids <- if (assoc_end>=assoc_start) seq.int(assoc_start,assoc_end) else integer()

summary_file <- file.path(partial_dir,sprintf("association_summary_collect_%03d.tsv",collect_id))
candidate_file <- file.path(partial_dir,sprintf("association_candidates_collect_%03d.tsv",collect_id))
qc_file <- file.path(partial_dir,sprintf("association_collect_%03d_qc.tsv",collect_id))
if (overwrite) unlink(c(summary_file,candidate_file,qc_file), force=TRUE)

required_cols <- c(
  "job_id","final_table_row_id","splitID","CpG","requested_model","cpg_chr","cpg_start",
  "data_chunk_id","resolved_chunk_id","chunk_source","cpg_row_in_chunk",
  "association_term","coefficient","std_error","statistic","df","p_value","FID_sd",
  "optimal_lambda","stage1_coefficient","stage1_p_value","selected_stage1",
  "final_fit_stage","n_train_total","n_train_complete","n_train_FIDs","n_cases",
  "n_controls","n_meth_observed","meth_sd","fit_status","error_message"
)

thresholds <- c(0.05,0.01,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
threshold_names <- c("n_p_le_0.05","n_p_le_0.01","n_p_le_1e-3","n_p_le_1e-4",
                     "n_p_le_1e-5","n_p_le_1e-6","n_p_le_1e-7","n_p_le_1e-8")

append_tsv <- function(dt,path){
  if (is.null(dt) || nrow(dt)==0L) return(invisible(NULL))
  ex <- file.exists(path)
  fwrite(dt,path,sep="\t",quote=FALSE,na="NA",append=ex,col.names=!ex)
}

run_start <- Sys.time()
n_found <- 0L; n_completed <- 0L; n_raw <- 0; n_fit_ok <- 0; n_finite <- 0; n_candidates <- 0
missing_jobs <- integer(); bad_qc_jobs <- integer()

cat("Collector",collect_id,"upstream jobs",assoc_start,"to",assoc_end,"\n")

for (k in seq_along(assoc_job_ids)) {
  jid <- assoc_job_ids[k]
  raw_file <- file.path(result_dir,sprintf("single_cpg_training_assoc_job_%04d.tsv",jid))
  raw_qc <- file.path(upstream_qc_dir,sprintf("single_cpg_training_assoc_job_%04d_qc.tsv",jid))
  cat(sprintf("[%03d] upstream %04d (%d/%d)\n",collect_id,jid,k,length(assoc_job_ids)))

  if (!file.exists(raw_file)) {
    missing_jobs <- c(missing_jobs,jid)
    if (!allow_incomplete) stop("Missing upstream result: ",raw_file)
    next
  }
  n_found <- n_found+1L

  completed <- FALSE
  if (file.exists(raw_qc)) {
    q <- tryCatch(fread(raw_qc), error=function(e) NULL)
    completed <- !is.null(q) && "job_status"%in%names(q) && nrow(q)>=1L && q$job_status[1]=="COMPLETED"
  }
  if (completed) n_completed <- n_completed+1L else bad_qc_jobs <- c(bad_qc_jobs,jid)
  if (!completed && !allow_incomplete) stop("Upstream QC not COMPLETED for job ",jid)

  hdr <- fread(raw_file,nrows=0L)
  miss <- setdiff(required_cols,names(hdr))
  if (length(miss)>0L) stop("Missing columns in ",raw_file,": ",paste(miss,collapse=", "))

  dt <- fread(raw_file,select=required_cols,showProgress=FALSE)
  dt[, source_assoc_job_id:=as.integer(jid)]
  dt[, requested_model:=toupper(as.character(requested_model))]
  dt[, splitID:=as.integer(splitID)]
  dt[, p_value:=suppressWarnings(as.numeric(p_value))]
  dt[, coefficient:=suppressWarnings(as.numeric(coefficient))]

  n_raw <- n_raw+nrow(dt)
  n_fit_ok <- n_fit_ok+sum(dt$fit_status=="FIT_OK",na.rm=TRUE)
  n_finite <- n_finite+sum(dt$fit_status=="FIT_OK" & is.finite(dt$p_value),na.rm=TRUE)

  sum_dt <- dt[, {
    ok <- fit_status=="FIT_OK"
    fp <- ok & is.finite(p_value)
    vals <- lapply(thresholds,function(th) sum(fp & p_value<=th,na.rm=TRUE))
    names(vals) <- threshold_names
    c(list(
      source_assoc_job_id=as.integer(jid),
      n_rows=.N,
      n_fit_ok=sum(ok,na.rm=TRUE),
      n_finite_p=sum(fp,na.rm=TRUE),
      n_model_error=sum(fit_status=="MODEL_ERROR",na.rm=TRUE),
      n_skipped=sum(!is.na(fit_status) & fit_status!="FIT_OK" & fit_status!="MODEL_ERROR"),
      min_p_value=if(any(fp)) min(p_value[fp],na.rm=TRUE) else NA_real_
    ),vals)
  },by=.(splitID,requested_model)]
  append_tsv(sum_dt,summary_file)

  cand <- dt[fit_status=="FIT_OK" & is.finite(p_value) & p_value<=candidate_p_max]
  if (nrow(cand)>0L) {
    setcolorder(cand,c("source_assoc_job_id",setdiff(names(cand),"source_assoc_job_id")))
    append_tsv(cand,candidate_file)
    n_candidates <- n_candidates+nrow(cand)
  }
  rm(dt,sum_dt,cand); invisible(gc(FALSE))
}

qc <- data.table(
  collect_id=collect_id,n_collect_jobs=n_collect_jobs,n_assoc_jobs=n_assoc_jobs,
  assoc_job_start=assoc_start,assoc_job_end=assoc_end,
  n_upstream_expected=length(assoc_job_ids),n_upstream_found=n_found,
  n_upstream_completed=n_completed,n_raw_rows_read=n_raw,n_fit_ok=n_fit_ok,
  n_finite_p=n_finite,candidate_p_max=candidate_p_max,n_candidate_rows=n_candidates,
  missing_upstream_jobs=if(length(missing_jobs)) paste(missing_jobs,collapse=",") else NA_character_,
  bad_upstream_qc_jobs=if(length(bad_qc_jobs)) paste(bad_qc_jobs,collapse=",") else NA_character_,
  summary_file=summary_file,candidate_file=candidate_file,
  elapsed_minutes=as.numeric(difftime(Sys.time(),run_start,units="mins")),
  completed_at=format(Sys.time(),tz="America/Toronto",usetz=TRUE),
  job_status=if(n_found==length(assoc_job_ids) && n_completed==length(assoc_job_ids)) "COMPLETED" else "COMPLETED_INCOMPLETE_UPSTREAM"
)
fwrite(qc,qc_file,sep="\t",quote=FALSE,na="NA")
print(qc)
