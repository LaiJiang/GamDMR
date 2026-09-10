#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

PATH_wk<-path.expand(Sys.getenv("PATH_WK","~/scratch/UQAC/meth/"))
root<-path.expand(Sys.getenv(
  "SINGLE_CPG_PREDICTION_OUTPUT",
  file.path(PATH_wk,"results/15_revision/5_cv/5_prediction/6_single_cpg_primary")
))
outdir<-file.path(root,"collected");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
expected<-as.integer(Sys.getenv("EXPECTED_SPLITS","100"))
allow<-Sys.getenv("ALLOW_INCOMPLETE","0")%in%c("1","TRUE","true","T","yes","YES")

read_all<-function(dir,fmt,label){
  files<-file.path(dir,sprintf(fmt,seq_len(expected)));missing<-files[!file.exists(files)]
  if(length(missing)&&!allow)stop(label,": missing ",length(missing)," files. First: ",missing[1L])
  files<-files[file.exists(files)]
  if(!length(files))return(data.table())
  rbindlist(lapply(files,function(f)fread(f,showProgress=FALSE)),use.names=TRUE,fill=TRUE)
}

perf<-read_all(file.path(root,"performance"),
               "single_cpg_primary_performance_split_%03d.tsv","performance")
pred<-read_all(file.path(root,"predictions"),
               "single_cpg_primary_predictions_split_%03d.tsv","predictions")
coefdt<-read_all(file.path(root,"coefficients"),
                 "single_cpg_primary_coefficients_split_%03d.tsv","coefficients")
fqc<-read_all(file.path(root,"feature_qc"),
              "single_cpg_primary_feature_qc_split_%03d.tsv","feature_qc")
sqc<-read_all(file.path(root,"split_qc"),
              "single_cpg_primary_split_qc_%03d.tsv","split_qc")

fwrite(perf,file.path(outdir,"single_cpg_primary_performance_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(pred,file.path(outdir,"single_cpg_primary_predictions_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(coefdt,file.path(outdir,"single_cpg_primary_coefficients_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(fqc,file.path(outdir,"single_cpg_primary_feature_qc_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(sqc,file.path(outdir,"single_cpg_primary_split_qc_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")

metrics<-intersect(c("AUROC_test","PR_AUC_test","PR_AUC_lift_test","Brier_test","log_loss_test",
                     "calibration_intercept_test","calibration_slope_test","sensitivity_test",
                     "specificity_test","balanced_accuracy_test"),names(perf))
rows<-list()
if(nrow(perf)){
  gr<-unique(perf[,.(DMR_method,model)])
  for(i in seq_len(nrow(gr))){
    g<-gr[i];d<-perf[DMR_method==g$DMR_method&model==g$model]
    r<-data.table(DMR_method=g$DMR_method,model=g$model,n_splits=uniqueN(d$splitID),
                  n_fallback_splits=sum(d$model_fallback%in%TRUE,na.rm=TRUE),
                  fallback_fraction=mean(d$model_fallback%in%TRUE,na.rm=TRUE),
                  mean_CpGs_selected=mean(d$n_CpG_selected_input,na.rm=TRUE),
                  median_CpGs_selected=median(d$n_CpG_selected_input,na.rm=TRUE),
                  mean_CpGs_usable=mean(d$n_CpG_usable,na.rm=TRUE),
                  median_CpGs_usable=median(d$n_CpG_usable,na.rm=TRUE),
                  mean_nonzero_CpGs=mean(d$n_methylation_features_nonzero,na.rm=TRUE))
    for(v in metrics){
      r[,(paste0(v,"_mean")):=mean(d[[v]],na.rm=TRUE)]
      r[,(paste0(v,"_sd")):=sd(d[[v]],na.rm=TRUE)]
    }
    rows[[length(rows)+1L]]<-r
  }
}
summ<-if(length(rows))rbindlist(rows,fill=TRUE)else data.table()
fwrite(summ,file.path(outdir,"single_cpg_primary_performance_summary.tsv"),sep="\t",quote=FALSE,na="NA")

qc<-data.table(
  expected_splits=expected,
  performance_splits=if(nrow(perf))uniqueN(perf$splitID)else 0L,
  prediction_splits=if(nrow(pred))uniqueN(pred$splitID)else 0L,
  coefficient_splits=if(nrow(coefdt))uniqueN(coefdt$splitID)else 0L,
  feature_qc_splits=if(nrow(fqc))uniqueN(fqc$splitID)else 0L,
  split_qc_splits=if(nrow(sqc))uniqueN(sqc$splitID)else 0L,
  observed_models=if(nrow(perf))paste(sort(unique(perf$model)),collapse=",")else NA_character_,
  total_performance_rows=nrow(perf),
  completed_at=format(Sys.time(),tz="America/Toronto",usetz=TRUE),
  job_status=if(nrow(perf)&&uniqueN(perf$splitID)==expected)"COMPLETED"else"COMPLETED_WITH_MISSING_SPLITS"
)
fwrite(qc,file.path(outdir,"single_cpg_primary_collection_qc.tsv"),sep="\t",quote=FALSE,na="NA")
print(summ);print(qc)
