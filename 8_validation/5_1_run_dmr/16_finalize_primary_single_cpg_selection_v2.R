#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

PATH_WK <- path.expand(Sys.getenv("PATH_WK","~/scratch/UQAC/meth/"))
collect_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_COLLECT_ROOT",
  file.path(PATH_WK,"results/15_revision/5_cv/6_single_cpg/2_collected_significant")
))
n_collect_jobs <- as.integer(Sys.getenv("N_COLLECT_JOBS","99"))
significance_p <- as.numeric(Sys.getenv("SIGNIFICANCE_P","1e-5"))
candidate_p_max <- as.numeric(Sys.getenv("CANDIDATE_P_MAX","0.001"))
bonf_alpha <- as.numeric(Sys.getenv("BONF_ALPHA","0.05"))
expected_splits <- as.integer(Sys.getenv("EXPECTED_SPLITS","100"))
allow_incomplete <- Sys.getenv("ALLOW_INCOMPLETE","0") %in% c("1","TRUE","true","T","yes","YES")

if (!is.finite(significance_p) || significance_p<=0 || significance_p>1) stop("Invalid SIGNIFICANCE_P.")
if (!is.finite(candidate_p_max) || significance_p>candidate_p_max)
  stop("SIGNIFICANCE_P must be <= CANDIDATE_P_MAX.")

partial_dir <- file.path(collect_root,"partial")
table_dir <- file.path(collect_root,"tables")
selected_dir <- file.path(collect_root,"selected_by_split")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(selected_dir,recursive=TRUE,showWarnings=FALSE)

old <- list.files(selected_dir,pattern="^single_cpg_selected_split_[0-9]{3}\\.tsv$",full.names=TRUE)
if(length(old)) unlink(old,force=TRUE)

run_start <- Sys.time()

# ----------------------------- validate Stage 1 ------------------------------
qc_files <- file.path(partial_dir,sprintf("association_collect_%03d_qc.tsv",seq_len(n_collect_jobs)))
missing_qc <- qc_files[!file.exists(qc_files)]
if(length(missing_qc) && !allow_incomplete) stop("Missing Stage-1 QC files: ",length(missing_qc))
qc_files <- qc_files[file.exists(qc_files)]
qc_all <- rbindlist(lapply(qc_files,fread),use.names=TRUE,fill=TRUE)
if(!allow_incomplete && any(qc_all$job_status!="COMPLETED"))
  stop("Some Stage-1 collector jobs are not COMPLETED.")

summary_files <- file.path(partial_dir,sprintf("association_summary_collect_%03d.tsv",seq_len(n_collect_jobs)))
missing_summary <- summary_files[!file.exists(summary_files)]
if(length(missing_summary) && !allow_incomplete) stop("Missing Stage-1 summary files: ",length(missing_summary))
summary_files <- summary_files[file.exists(summary_files)]

candidate_files <- file.path(partial_dir,sprintf("association_candidates_collect_%03d.tsv",seq_len(n_collect_jobs)))
candidate_files <- candidate_files[file.exists(candidate_files)]

cat("Stage-1 QC:",length(qc_files),"\n")
cat("Stage-1 summaries:",length(summary_files),"\n")
cat("Candidate files:",length(candidate_files),"\n")
cat("Primary p threshold:",significance_p,"\n")

# ---------------------- exact split/model significance counts ----------------
sum_all <- rbindlist(lapply(summary_files,function(f) fread(f,showProgress=FALSE)),
                     use.names=TRUE,fill=TRUE)

count_cols <- c("n_rows","n_fit_ok","n_finite_p","n_model_error","n_skipped",
                "n_p_le_0.05","n_p_le_0.01","n_p_le_1e-3","n_p_le_1e-4",
                "n_p_le_1e-5","n_p_le_1e-6","n_p_le_1e-7","n_p_le_1e-8")

for(cc in intersect(count_cols,names(sum_all)))
  set(sum_all,j=cc,value=as.numeric(sum_all[[cc]]))

summary_split_model <- sum_all[, {
  z <- lapply(.SD,sum,na.rm=TRUE)
  z$min_p_value <- if(all(is.na(min_p_value))) NA_real_ else min(min_p_value,na.rm=TRUE)
  z
},by=.(splitID,requested_model),.SDcols=count_cols]

summary_split_model[, bonferroni_threshold :=
  fifelse(n_finite_p>0,bonf_alpha/n_finite_p,NA_real_)]
setorder(summary_split_model,splitID,requested_model)

fwrite(summary_split_model,
       file.path(table_dir,"association_significance_summary_by_split_model.tsv"),
       sep="\t",quote=FALSE,na="NA")

summary_overall <- summary_split_model[, .(
  n_splits=uniqueN(splitID),
  mean_n_tested=mean(n_finite_p,na.rm=TRUE),
  min_n_tested=min(n_finite_p,na.rm=TRUE),
  max_n_tested=max(n_finite_p,na.rm=TRUE),
  mean_n_p_le_0_05=mean(get("n_p_le_0.05"),na.rm=TRUE),
  mean_n_p_le_0_01=mean(get("n_p_le_0.01"),na.rm=TRUE),
  mean_n_p_le_1e_3=mean(get("n_p_le_1e-3"),na.rm=TRUE),
  mean_n_p_le_1e_4=mean(get("n_p_le_1e-4"),na.rm=TRUE),
  mean_n_p_le_1e_5=mean(get("n_p_le_1e-5"),na.rm=TRUE),
  median_n_p_le_1e_5=median(get("n_p_le_1e-5"),na.rm=TRUE),
  min_n_p_le_1e_5=min(get("n_p_le_1e-5"),na.rm=TRUE),
  max_n_p_le_1e_5=max(get("n_p_le_1e-5"),na.rm=TRUE),
  mean_n_p_le_1e_6=mean(get("n_p_le_1e-6"),na.rm=TRUE),
  mean_n_p_le_1e_7=mean(get("n_p_le_1e-7"),na.rm=TRUE),
  mean_n_p_le_1e_8=mean(get("n_p_le_1e-8"),na.rm=TRUE),
  median_min_p=median(min_p_value,na.rm=TRUE)
),by=requested_model]
setorder(summary_overall,requested_model)

fwrite(summary_overall,
       file.path(table_dir,"association_significance_summary_overall_model.tsv"),
       sep="\t",quote=FALSE,na="NA")

# -------------------------- compact candidate tail ---------------------------
if(length(candidate_files)){
  cand <- rbindlist(lapply(candidate_files,function(f) fread(f,showProgress=FALSE)),
                    use.names=TRUE,fill=TRUE)
  cand[, splitID:=as.integer(splitID)]
  cand[, requested_model:=toupper(as.character(requested_model))]
  cand[, p_value:=as.numeric(p_value)]
  cand[, coefficient:=as.numeric(coefficient)]
  cand <- cand[fit_status=="FIT_OK" & is.finite(p_value)]
}else{
  cand <- data.table()
}

# Bonferroni is exact because all such hits must be far below CANDIDATE_P_MAX.
if (nrow(cand) > 0L) {
  bonf_map <- summary_split_model[
    ,
    .(splitID, requested_model, bonferroni_threshold)
  ]
  bonf <- merge(
    cand,
    bonf_map,
    by = c("splitID", "requested_model")
  )
  bonf <- bonf[
    is.finite(bonferroni_threshold) &
      p_value <= bonferroni_threshold
  ]
} else {
  bonf <- data.table()
}

fwrite(bonf,file.path(table_dir,"bonferroni_significant_associations.tsv"),
       sep="\t",quote=FALSE,na="NA")

# -------------------------- primary selection rule ---------------------------
sig <- if (nrow(cand) > 0L) {
  cand[p_value <= significance_p]
} else {
  data.table()
}
fwrite(sig,file.path(table_dir,"primary_significant_associations.tsv"),
       sep="\t",quote=FALSE,na="NA")

models <- c("M1","M2","M3")
split_ids <- sort(unique(summary_split_model$splitID))

sig_counts <- if(nrow(sig)){
  sig[,.(n_significant_associations=.N,n_unique_cpgs=uniqueN(CpG)),
      by=.(splitID,requested_model)]
} else {
  data.table(
    splitID = integer(),
    requested_model = character(),
    n_significant_associations = integer(),
    n_unique_cpgs = integer()
  )
}

grid <- CJ(splitID=split_ids,requested_model=models,unique=TRUE)
sig_counts <- merge(grid,sig_counts,by=c("splitID","requested_model"),all.x=TRUE)
sig_counts[is.na(n_significant_associations),
           `:=`(n_significant_associations=0L,n_unique_cpgs=0L)]
fwrite(sig_counts,file.path(table_dir,"primary_significant_counts_by_split_model.tsv"),
       sep="\t",quote=FALSE,na="NA")

# -------------------- deduplicated CpG union within split --------------------
template <- data.table(
  splitID=integer(),CpG=character(),cpg_chr=character(),cpg_start=integer(),
  data_chunk_id=integer(),resolved_chunk_id=integer(),cpg_row_in_chunk=integer(),
  min_p_value=numeric(),best_model=character(),best_coefficient=numeric(),
  n_models_significant=integer(),significant_models=character(),
  p_M1=numeric(),p_M2=numeric(),p_M3=numeric(),
  coef_M1=numeric(),coef_M2=numeric(),coef_M3=numeric()
)

if(nrow(sig)){
  setorder(sig,splitID,CpG,requested_model,p_value)
  su <- sig[,.SD[1L],by=.(splitID,CpG,requested_model)]

  best <- su[order(splitID,CpG,p_value),.SD[1L],by=.(splitID,CpG)][,
    .(splitID,CpG,cpg_chr,cpg_start,data_chunk_id,resolved_chunk_id,cpg_row_in_chunk,
      min_p_value=p_value,best_model=requested_model,best_coefficient=coefficient)]

  modagg <- su[,.(n_models_significant=uniqueN(requested_model),
                  significant_models=paste(sort(unique(requested_model)),collapse=";")),
               by=.(splitID,CpG)]

  pw <- dcast(su,splitID+CpG~requested_model,value.var="p_value",
              fun.aggregate=min,fill=NA_real_)
  for(m in models) if(!m%in%names(pw)) pw[,(m):=NA_real_]
  setnames(pw,models,paste0("p_",models))

  cw <- dcast(su,splitID+CpG~requested_model,value.var="coefficient",
              fun.aggregate=function(x)x[1L],fill=NA_real_)
  for(m in models) if(!m%in%names(cw)) cw[,(m):=NA_real_]
  setnames(cw,models,paste0("coef_",models))

  selected <- Reduce(function(x,y) merge(x,y,by=c("splitID","CpG"),all=TRUE),
                     list(best,modagg,pw,cw))
  selected <- selected[,names(template),with=FALSE]
  setorder(selected,splitID,min_p_value,CpG)
} else {
  selected <- copy(template)
}

fwrite(selected,file.path(table_dir,"primary_selected_cpg_union_all_splits.tsv"),
       sep="\t",quote=FALSE,na="NA")

# Per-split files are the ONLY files to feed the outer train/test prediction.
for(sid in split_ids){
  z <- selected[splitID==sid]
  if(nrow(z)==0L) z <- copy(template)
  fwrite(z,file.path(selected_dir,sprintf("single_cpg_selected_split_%03d.tsv",sid)),
         sep="\t",quote=FALSE,na="NA")
}

# ----------------------------- split summaries -------------------------------
union_counts <- selected[,.(n_unique_selected_cpgs=.N),by=splitID]
split_summary <- merge(data.table(splitID=split_ids),union_counts,by="splitID",all.x=TRUE)
split_summary[is.na(n_unique_selected_cpgs),n_unique_selected_cpgs:=0L]

for(m in models){
  x <- sig_counts[requested_model==m,.(splitID,n_unique_cpgs)]
  setnames(x,"n_unique_cpgs",paste0("n_selected_",m))
  split_summary <- merge(split_summary,x,by="splitID",all.x=TRUE)
}
for(cc in paste0("n_selected_",models))
  split_summary[is.na(get(cc)),(cc):=0L]
setorder(split_summary,splitID)

fwrite(split_summary,file.path(table_dir,"primary_selected_cpg_counts_by_split.tsv"),
       sep="\t",quote=FALSE,na="NA")

patterns <- if (nrow(selected) > 0L) {
  selected[
    ,
    .(n_cpgs = .N),
    by = .(splitID, significant_models)
  ][order(splitID, significant_models)]
} else {
  data.table(
    splitID = integer(),
    significant_models = character(),
    n_cpgs = integer()
  )
}

fwrite(patterns,file.path(table_dir,"primary_selected_cpg_model_overlap_by_split.tsv"),
       sep="\t",quote=FALSE,na="NA")

# DESCRIPTIVE ONLY. Never use across-split frequency to choose predictors.
if(nrow(selected)){
  stability <- selected[,.(n_splits_selected=uniqueN(splitID),
                           selection_frequency=uniqueN(splitID)/length(split_ids),
                           min_p_across_selected_splits=min(min_p_value,na.rm=TRUE),
                           median_p_across_selected_splits=median(min_p_value,na.rm=TRUE),
                           models_ever_significant=paste(sort(unique(unlist(strsplit(significant_models,";",fixed=TRUE)))),collapse=";"),
                           cpg_chr=cpg_chr[1L],cpg_start=cpg_start[1L],
                           data_chunk_id=data_chunk_id[1L],
                           cpg_row_in_chunk=cpg_row_in_chunk[1L]),by=CpG]
  setorder(stability,-selection_frequency,min_p_across_selected_splits,CpG)
} else {
  stability <- data.table()
}

fwrite(stability,
       file.path(table_dir,"primary_selected_cpg_stability_DESCRIPTIVE_ONLY.tsv"),
       sep="\t",quote=FALSE,na="NA")

# -------------------------------- final QC -----------------------------------
final_qc <- data.table(
  n_collect_jobs_expected=n_collect_jobs,
  n_collect_qc_files=length(qc_files),
  n_collect_summary_files=length(summary_files),
  n_candidate_files=length(candidate_files),
  n_observed_splits=length(split_ids),
  expected_splits=expected_splits,
  observed_models=paste(sort(unique(summary_split_model$requested_model)),collapse=","),
  candidate_p_max=candidate_p_max,
  primary_significance_p=significance_p,
  bonferroni_alpha=bonf_alpha,
  n_candidate_rows_loaded=nrow(cand),
  n_primary_significant_association_rows=nrow(sig),
  n_primary_unique_split_cpg_predictors=nrow(selected),
  n_unique_cpgs_selected_in_any_split=if(nrow(selected)) uniqueN(selected$CpG) else 0L,
  min_selected_cpgs_per_split=min(split_summary$n_unique_selected_cpgs),
  median_selected_cpgs_per_split=median(split_summary$n_unique_selected_cpgs),
  max_selected_cpgs_per_split=max(split_summary$n_unique_selected_cpgs),
  selected_by_split_dir=selected_dir,
  elapsed_minutes=as.numeric(difftime(Sys.time(),run_start,units="mins")),
  completed_at=format(Sys.time(),tz="America/Toronto",usetz=TRUE),
  job_status = if (
    length(qc_files) == n_collect_jobs &&
      length(split_ids) == expected_splits
  ) {
    "COMPLETED"
  } else {
    "COMPLETED_WITH_QC_WARNING"
  }
)
fwrite(final_qc,file.path(collect_root,"primary_single_cpg_selection_qc.tsv"),
       sep="\t",quote=FALSE,na="NA")

cat("\nPrimary single-CpG selection finalized.\n")
print(final_qc)
cat("\nSelected CpG counts by split:\n")
print(split_summary)
cat("\nThe stability table is DESCRIPTIVE ONLY; do not use it to define outer-split predictors.\n")
