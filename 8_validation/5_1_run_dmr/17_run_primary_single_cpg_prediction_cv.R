#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
  library(PRROC)
})
data.table::setDTthreads(1L)
options(mc.cores=1L)

# One array task = one outer family split.
# M1/M2/M3 CpG sets are predicted separately. A zero-CpG set is retained via
# the covariate-only fallback so all 100 paired outer splits remain evaluable.

PATH_wk <- path.expand(Sys.getenv("PATH_WK","~/scratch/UQAC/meth/"))
selected_cpg_dir <- path.expand(Sys.getenv(
  "SINGLE_CPG_SELECTED_DIR",
  file.path(PATH_wk,"results/15_revision/5_cv/6_single_cpg/2_collected_significant/selected_by_split")
))
train_split_file <- path.expand(Sys.getenv(
  "TRAIN_SPLIT_FILE",file.path(PATH_wk,"scr/11_mgcv/dat/100_family_training_splits.csv")
))
pheno_file_path <- path.expand(Sys.getenv(
  "PHENO_FILE",file.path(PATH_wk,"scr/11_mgcv/dat/18_pheno_BMI.RData")
))
chunk_dir <- path.expand(Sys.getenv(
  "METH_CHUNK_DIR",file.path(PATH_wk,"data/meth_split")
))
output_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_PREDICTION_OUTPUT",
  file.path(PATH_wk,"results/15_revision/5_cv/5_prediction/6_single_cpg_primary")
))
shared_fold_dir <- file.path(PATH_wk,"results/15_revision/5_cv/5_prediction/shared_inner_folds")
performance_dir <- file.path(output_root,"performance")
prediction_dir <- file.path(output_root,"predictions")
coefficient_dir <- file.path(output_root,"coefficients")
feature_qc_dir <- file.path(output_root,"feature_qc")
split_qc_dir <- file.path(output_root,"split_qc")
for(d in c(performance_dir,prediction_dir,coefficient_dir,feature_qc_dir,split_qc_dir,shared_fold_dir))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)

inner_folds_requested <- as.integer(Sys.getenv("INNER_FOLDS","5"))
seed_base <- as.integer(Sys.getenv("SEED_BASE","20260803"))
min_obs <- as.numeric(Sys.getenv("MIN_TRAIN_CPG_OBSERVED_FRACTION","0.50"))
variance_epsilon <- as.numeric(Sys.getenv("VARIANCE_EPSILON","1e-8"))
max_cpgs_per_model <- as.integer(Sys.getenv("MAX_CPGS_PER_MODEL","0"))
overwrite <- Sys.getenv("OVERWRITE","0") %in% c("1","TRUE","true","T","yes","YES")
if(is.na(inner_folds_requested)||inner_folds_requested<3L) stop("INNER_FOLDS must be >=3.")
if(!is.finite(min_obs)||min_obs<=0||min_obs>1) stop("MIN_TRAIN_CPG_OBSERVED_FRACTION must be in (0,1].")
if(!is.finite(variance_epsilon)||variance_epsilon<0) stop("Invalid VARIANCE_EPSILON.")
if(is.na(max_cpgs_per_model)||max_cpgs_per_model<0L) stop("MAX_CPGS_PER_MODEL must be >=0.")

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1L) stop("Pass splitID.")
split_id <- as.integer(args[1L])
if(is.na(split_id)||split_id<1L||split_id>100L) stop("splitID must be 1...100.")
run_start <- Sys.time()
cat("Primary single-CpG prediction, outer split",split_id,"\n")

performance_file <- file.path(performance_dir,sprintf("single_cpg_primary_performance_split_%03d.tsv",split_id))
prediction_file <- file.path(prediction_dir,sprintf("single_cpg_primary_predictions_split_%03d.tsv",split_id))
coefficient_file <- file.path(coefficient_dir,sprintf("single_cpg_primary_coefficients_split_%03d.tsv",split_id))
feature_qc_file <- file.path(feature_qc_dir,sprintf("single_cpg_primary_feature_qc_split_%03d.tsv",split_id))
split_qc_file <- file.path(split_qc_dir,sprintf("single_cpg_primary_split_qc_%03d.tsv",split_id))
fold_file <- file.path(shared_fold_dir,sprintf("inner_family_folds_split_%03d.tsv",split_id))
outs <- c(performance_file,prediction_file,coefficient_file,feature_qc_file,split_qc_file)
if(!overwrite && all(file.exists(outs))) quit(save="no",status=0L)
if(overwrite) unlink(outs,force=TRUE)

normalize_id <- function(x){x<-as.character(x);x<-sub("^X","",x);gsub("\\.","-",x)}
normalize_chr <- function(x){x<-as.character(x);x<-sub("^chr","",x,ignore.case=TRUE);sub("\\.0$","",x)}
safe_feature_name <- function(x) paste0("CPG_",gsub("[^A-Za-z0-9]+","_",as.character(x)))
clean_probability <- function(p,eps=1e-8) pmin(pmax(as.numeric(p),eps),1-eps)
safe_divide <- function(a,b) if(!is.finite(b)||b==0) NA_real_ else as.numeric(a/b)

safe_roc_auc <- function(y,p){
  y<-as.integer(y);p<-as.numeric(p);k<-is.finite(y)&is.finite(p);y<-y[k];p<-p[k]
  if(length(unique(y))<2L) return(NA_real_)
  z<-tryCatch(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE),error=function(e)NULL)
  if(is.null(z)) NA_real_ else as.numeric(pROC::auc(z))
}
safe_pr_auc <- function(y,p){
  y<-as.integer(y);p<-as.numeric(p);k<-is.finite(y)&is.finite(p);y<-y[k];p<-p[k]
  if(length(unique(y))<2L) return(NA_real_)
  z<-tryCatch(PRROC::pr.curve(scores.class0=p[y==1L],scores.class1=p[y==0L],curve=FALSE),
              error=function(e)NULL)
  if(is.null(z)) NA_real_ else as.numeric(z$auc.integral)
}
choose_youden_threshold <- function(y,p){
  y<-as.integer(y);p<-as.numeric(p);k<-is.finite(y)&is.finite(p);y<-y[k];p<-p[k]
  if(length(unique(y))<2L||length(unique(p))<2L) return(0.5)
  z<-tryCatch(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE),error=function(e)NULL)
  if(is.null(z)) return(0.5)
  q<-tryCatch(pROC::coords(z,x="best",best.method="youden",ret="threshold",transpose=FALSE),
              error=function(e)NULL)
  if(is.null(q)) return(0.5)
  q<-suppressWarnings(as.numeric(unlist(q,use.names=FALSE)));q<-q[is.finite(q)]
  if(!length(q)) 0.5 else q[1L]
}
calculate_metrics <- function(y,p,threshold){
  y<-as.integer(y);p<-as.numeric(p);k<-is.finite(y)&is.finite(p);y<-y[k];p<-p[k]
  if(!length(y)){
    nn<-c("n","n_case","n_control","prevalence","AUROC","PR_AUC","PR_baseline","PR_AUC_lift",
          "Brier","log_loss","calibration_intercept","calibration_slope","sensitivity",
          "specificity","balanced_accuracy","PPV","NPV")
    x<-as.list(rep(NA_real_,length(nn)));names(x)<-nn;return(x)
  }
  pc<-as.integer(p>=threshold);tp<-sum(pc==1&y==1);tn<-sum(pc==0&y==0)
  fp<-sum(pc==1&y==0);fn<-sum(pc==0&y==1)
  se<-safe_divide(tp,tp+fn);sp<-safe_divide(tn,tn+fp);ppv<-safe_divide(tp,tp+fp);npv<-safe_divide(tn,tn+fn)
  ba<-mean(c(se,sp),na.rm=TRUE);if(!is.finite(ba))ba<-NA_real_
  pp<-clean_probability(p);lp<-qlogis(pp);ci<-NA_real_;cs<-NA_real_
  if(length(unique(y))>=2L&&is.finite(sd(lp))&&sd(lp)>0){
    fi<-tryCatch(suppressWarnings(glm(y~1+offset(lp),family=binomial())),error=function(e)NULL)
    fs<-tryCatch(suppressWarnings(glm(y~lp,family=binomial())),error=function(e)NULL)
    if(!is.null(fi))ci<-unname(coef(fi)[1L])
    if(!is.null(fs)&&length(coef(fs))>=2L)cs<-unname(coef(fs)[2L])
  }
  pr<-safe_pr_auc(y,p)
  list(n=length(y),n_case=sum(y==1L),n_control=sum(y==0L),prevalence=mean(y),
       AUROC=safe_roc_auc(y,p),PR_AUC=pr,PR_baseline=mean(y),PR_AUC_lift=pr-mean(y),
       Brier=mean((p-y)^2),log_loss=-mean(y*log(pp)+(1-y)*log(1-pp)),
       calibration_intercept=ci,calibration_slope=cs,sensitivity=se,specificity=sp,
       balanced_accuracy=ba,PPV=ppv,NPV=npv)
}
metrics_to_columns <- function(x,suffix){
  n<-names(x);names(x)<-ifelse(n=="n",paste0("n_evaluable_",suffix),paste0(n,"_",suffix));x
}
most_frequent <- function(x){
  x<-x[!is.na(x)&nzchar(as.character(x))];if(!length(x))return(NA_character_)
  names(sort(table(as.character(x)),decreasing=TRUE))[1L]
}
prep_num <- function(tr,te,name){
  tr<-suppressWarnings(as.numeric(as.character(tr)));te<-suppressWarnings(as.numeric(as.character(te)))
  med<-median(tr[is.finite(tr)],na.rm=TRUE);if(!is.finite(med))med<-0
  tr[!is.finite(tr)]<-med;te[!is.finite(te)]<-med
  list(train=matrix(tr,ncol=1,dimnames=list(NULL,name)),test=matrix(te,ncol=1,dimnames=list(NULL,name)))
}
prep_sex <- function(tr,te){
  tc<-as.character(tr);ec<-as.character(te);tn<-suppressWarnings(as.numeric(tc));en<-suppressWarnings(as.numeric(ec))
  ok<-!is.na(tc)&nzchar(tc)
  if(all(!ok|is.finite(tn)))return(prep_num(tn,en,"COV_Sex"))
  mode<-most_frequent(tc);if(is.na(mode))mode<-"Unknown"
  tc[is.na(tc)|!nzchar(tc)]<-mode;ec[is.na(ec)|!nzchar(ec)]<-mode
  lev<-sort(unique(tc));ec[!ec%in%lev]<-mode
  if(length(lev)<=1L)return(list(train=matrix(numeric(0),nrow=length(tc),ncol=0),
                                test=matrix(numeric(0),nrow=length(ec),ncol=0)))
  dl<-lev[-1L];a<-sapply(dl,function(z)as.numeric(tc==z));b<-sapply(dl,function(z)as.numeric(ec==z))
  if(is.null(dim(a))){a<-matrix(a,ncol=1);b<-matrix(b,ncol=1)}
  colnames(a)<-paste0("COV_Sex_",make.names(dl));colnames(b)<-colnames(a);list(train=a,test=b)
}
make_covariates <- function(tr,te){
  a<-prep_num(tr$AgeCalc,te$AgeCalc,"COV_AgeCalc")
  s<-prep_sex(tr$Sex,te$Sex)
  n<-prep_num(tr$Non.smoker,te$Non.smoker,"COV_NonSmoker")
  b<-prep_num(tr$BMI,te$BMI,"COV_BMI")
  x<-do.call(cbind,list(a$train,s$train,n$train,b$train))
  z<-do.call(cbind,list(a$test,s$test,n$test,b$test))
  rownames(x)<-tr$ID;rownames(z)<-te$ID;list(train=x,test=z)
}

make_family_inner_folds <- function(train_dt,requested_k,seed,fold_path){
  if(file.exists(fold_path)){
    x<-fread(fold_path);req<-c("splitID","FID","inner_fold")
    if(!all(req%in%names(x)))stop("Incompatible existing inner-fold file.")
    x[,FID:=as.character(FID)]
    if(!setequal(x$FID,sort(unique(as.character(train_dt$FID)))))
      stop("Existing inner-fold FIDs do not match training FIDs.")
    return(list(foldid=as.integer(x$inner_fold[match(as.character(train_dt$FID),x$FID)]),
                n_folds=max(x$inner_fold),family_folds=x))
  }
  fam<-train_dt[,.(fam_AA=max(AA_only),n_subjects=.N,n_cases=sum(AA_only==1L),
                    n_controls=sum(AA_only==0L)),by=FID]
  fam[,FID:=as.character(FID)]
  mk<-min(requested_k,fam[fam_AA==1L,.N],fam[fam_AA==0L,.N])
  if(!is.finite(mk)||mk<3L)stop("Too few case/control families for inner CV.")
  chosen<-NULL
  for(k in seq.int(mk,3L,by=-1L)){
    for(attempt in seq_len(200L)){
      set.seed(seed+1000L*k+attempt);z<-copy(fam);z[,inner_fold:=NA_integer_]
      for(cls in c(0L,1L)){
        idx<-which(z$fam_AA==cls);idx<-sample(idx,length(idx))
        z$inner_fold[idx]<-rep(seq_len(k),length.out=length(idx))
      }
      sf<-z$inner_fold[match(as.character(train_dt$FID),z$FID)]
      good<-all(vapply(seq_len(k),function(j)length(unique(train_dt$AA_only[sf==j]))==2L,logical(1)))
      if(good){chosen<-z;break}
    }
    if(!is.null(chosen))break
  }
  if(is.null(chosen))stop("Could not construct family-grouped inner folds.")
  out<-chosen[,.(splitID=split_id,FID,fam_AA,n_subjects,n_cases,n_controls,inner_fold)]
  setorder(out,inner_fold,FID);fwrite(out,fold_path,sep="\t",quote=FALSE,na="NA")
  list(foldid=as.integer(out$inner_fold[match(as.character(train_dt$FID),out$FID)]),
       n_folds=max(out$inner_fold),family_folds=out)
}

# ------------------------------- inputs --------------------------------------
selected_cpg_file <- file.path(selected_cpg_dir,sprintf("single_cpg_selected_split_%03d.tsv",split_id))
if(!file.exists(selected_cpg_file))stop("Missing selected CpG file: ",selected_cpg_file)
sel<-fread(selected_cpg_file)
req<-c("splitID","CpG","cpg_chr","cpg_start","data_chunk_id","cpg_row_in_chunk","p_M1","p_M2","p_M3")
miss<-setdiff(req,names(sel));if(length(miss))stop("Selected CpG file missing: ",paste(miss,collapse=", "))
sel<-sel[splitID==split_id]
sel[,cpg_chr:=normalize_chr(cpg_chr)]
sel[,cpg_start:=as.integer(cpg_start)]
sel[,data_chunk_id:=as.integer(data_chunk_id)]
sel[,cpg_row_in_chunk:=as.integer(cpg_row_in_chunk)]
for(cc in c("p_M1","p_M2","p_M3"))set(sel,j=cc,value=suppressWarnings(as.numeric(sel[[cc]])))
sel<-unique(sel,by="CpG")

model_sets<-list(M1=sel[is.finite(p_M1)],M2=sel[is.finite(p_M2)],M3=sel[is.finite(p_M3)])
if(max_cpgs_per_model>0L){
  for(m in names(model_sets)){
    z<-copy(model_sets[[m]]);pc<-paste0("p_",m)
    if(nrow(z)>max_cpgs_per_model){setorderv(z,c(pc,"CpG"),c(1L,1L),na.last=TRUE);z<-head(z,max_cpgs_per_model)}
    model_sets[[m]]<-z
  }
}
for(m in names(model_sets))cat("Selected",m,":",nrow(model_sets[[m]]),"\n")
union_cpgs<-rbindlist(model_sets,use.names=TRUE,fill=TRUE)
if(nrow(union_cpgs)){union_cpgs<-unique(union_cpgs,by="CpG");setorder(union_cpgs,data_chunk_id,cpg_chr,cpg_start,CpG)}

split_tab<-fread(train_split_file)
train_cols<-grep("^train_FID_[0-9]+$",names(split_tab),value=TRUE)
sr<-split_tab[splitID==split_id];if(nrow(sr)!=1L)stop("Expected one split row.")
train_fids<-unique(trimws(as.character(unlist(sr[,..train_cols],use.names=FALSE))))
train_fids<-train_fids[!is.na(train_fids)&nzchar(train_fids)&train_fids!="NA"]

e<-new.env(parent=emptyenv());load(pheno_file_path,envir=e)
if(!exists("pheno_file",envir=e,inherits=FALSE))stop("pheno_file not found.")
ph<-as.data.table(copy(get("pheno_file",envir=e)));rm(e)
rp<-c("ID","FID","AA_only","AgeCalc","Sex","Non.smoker","BMI")
if(length(setdiff(rp,names(ph))))stop("Phenotype columns missing.")
ph<-ph[,..rp];ph[,ID:=normalize_id(ID)];ph[,FID:=as.character(FID)]
ph[,AA_only:=suppressWarnings(as.integer(as.character(AA_only)))]
ph<-ph[!is.na(ID)&nzchar(ID)&!is.na(FID)&nzchar(FID)&AA_only%in%c(0L,1L)]
if(anyDuplicated(ph$ID))stop("Duplicated phenotype IDs.")
train_dt<-ph[FID%chin%train_fids];test_dt<-ph[!FID%chin%train_fids]
setorder(train_dt,ID);setorder(test_dt,ID)
if(!nrow(train_dt)||!nrow(test_dt))stop("Empty outer train/test set.")
if(length(unique(train_dt$AA_only))<2L||length(unique(test_dt$AA_only))<2L)stop("Both classes required.")
fam_overlap<-length(intersect(train_dt$FID,test_dt$FID));if(fam_overlap)stop("Family leakage detected.")
y_train<-as.integer(train_dt$AA_only);y_test<-as.integer(test_dt$AA_only)
covars<-make_covariates(train_dt,test_dt)
fold_obj<-make_family_inner_folds(train_dt,inner_folds_requested,seed_base+split_id,fold_file)
foldid<-fold_obj$foldid
cat("Train:",nrow(train_dt),"subjects,",uniqueN(train_dt$FID),"families. Test:",nrow(test_dt),"subjects.\n")

# -------------------------- CpG feature construction -------------------------
empty_features <- function(){
  q <- data.table(
    splitID=integer(),CpG=character(),feature=character(),cpg_chr=character(),
    cpg_start=integer(),data_chunk_id=integer(),cpg_row_in_chunk=integer(),
    feature_status=character(),feature_message=character(),train_observed_n=integer(),
    train_observed_fraction=numeric(),test_observed_n=integer(),
    training_median_asr=numeric(),training_sd_after_imputation=numeric()
  )
  list(train=matrix(numeric(0),nrow=nrow(train_dt),ncol=0,dimnames=list(train_dt$ID,NULL)),
       test=matrix(numeric(0),nrow=nrow(test_dt),ncol=0,dimnames=list(test_dt$ID,NULL)),
       qc=q)
}
build_cpg_features <- function(cpgs){
  if(!nrow(cpgs))return(empty_features())
  trlist<-list();telist<-list();qclist<-list();all_ids<-unique(c(train_dt$ID,test_dt$ID))
  addqc<-function(info,status,msg=NA_character_,ntr=NA_integer_,frac=NA_real_,nte=NA_integer_,med=NA_real_,sdv=NA_real_){
    qclist[[length(qclist)+1L]]<<-data.table(
      splitID=split_id,CpG=as.character(info$CpG),feature=safe_feature_name(info$CpG),
      cpg_chr=as.character(info$cpg_chr),cpg_start=as.integer(info$cpg_start),
      data_chunk_id=as.integer(info$data_chunk_id),cpg_row_in_chunk=as.integer(info$cpg_row_in_chunk),
      feature_status=status,feature_message=msg,train_observed_n=ntr,
      train_observed_fraction=frac,test_observed_n=nte,training_median_asr=med,
      training_sd_after_imputation=sdv)
  }
  chunks<-split(cpgs,cpgs$data_chunk_id)
  for(cn in names(chunks)){
    z<-as.data.table(chunks[[cn]]);cid<-as.integer(cn)
    ff<-file.path(chunk_dir,sprintf("chunk_%04d.csv",cid))
    if(!file.exists(ff)){for(i in seq_len(nrow(z)))addqc(z[i],"MISSING_CHUNK_FILE",ff);next}
    h<-tryCatch(fread(ff,nrows=0L,showProgress=FALSE),error=function(e)e)
    if(inherits(h,"error")){for(i in seq_len(nrow(z)))addqc(z[i],"CHUNK_HEADER_ERROR",conditionMessage(h));next}
    if(!all(c("chr","start")%in%names(h))){for(i in seq_len(nrow(z)))addqc(z[i],"CHUNK_COORDINATE_ERROR");next}
    mc<-grep("_meth$",names(h),value=TRUE)
    mmap<-data.table(source_column=mc,ID=normalize_id(sub("_meth$","",mc)))[ID%chin%all_ids]
    if(anyDuplicated(mmap$ID))stop("Duplicated methylation IDs in chunk ",cid)
    if(!nrow(mmap)){for(i in seq_len(nrow(z)))addqc(z[i],"NO_MATCHING_METHYLATION_SAMPLES");next}
    dt<-tryCatch(fread(ff,select=c("chr","start",mmap$source_column),showProgress=FALSE),error=function(e)e)
    if(inherits(dt,"error")){for(i in seq_len(nrow(z)))addqc(z[i],"CHUNK_READ_ERROR",conditionMessage(dt));next}
    dt[,chr:=normalize_chr(chr)];dt[,start:=as.integer(start)];src<-mmap$source_column
    for(i in seq_len(nrow(z))){
      info<-z[i];rr<-as.integer(info$cpg_row_in_chunk)
      ok<-is.finite(rr)&&rr>=1L&&rr<=nrow(dt)&&as.character(dt$chr[rr])==as.character(info$cpg_chr)&&
          !is.na(dt$start[rr])&&dt$start[rr]==as.integer(info$cpg_start)
      if(!ok){
        hit<-which(dt$chr==as.character(info$cpg_chr)&dt$start==as.integer(info$cpg_start))
        if(!length(hit)){addqc(info,"CPG_NOT_FOUND");next};rr<-hit[1L]
      }
      vm<-as.matrix(dt[rr,..src]);storage.mode(vm)<-"numeric";v<-as.numeric(vm[1L,])
      v[!is.finite(v)|v<0|v>1]<-NA_real_;names(v)<-mmap$ID
      a<-v[match(all_ids,names(v))];names(a)<-all_ids;a<-asin(sqrt(a))
      tr<-a[train_dt$ID];te<-a[test_dt$ID];ntr<-sum(is.finite(tr));nte<-sum(is.finite(te))
      frac<-mean(is.finite(tr))
      if(!is.finite(frac)||frac<min_obs){addqc(info,"LOW_TRAIN_OBSERVED_FRACTION",ntr=ntr,frac=frac,nte=nte);next}
      med<-median(tr[is.finite(tr)],na.rm=TRUE)
      if(!is.finite(med)){addqc(info,"NO_FINITE_TRAINING_MEDIAN",ntr=ntr,frac=frac,nte=nte);next}
      tr[!is.finite(tr)]<-med;te[!is.finite(te)]<-med;sdv<-sd(tr)
      if(!is.finite(sdv)||sdv<=variance_epsilon){addqc(info,"ZERO_TRAIN_VARIANCE",ntr=ntr,frac=frac,nte=nte,med=med,sdv=sdv);next}
      fn<-safe_feature_name(info$CpG);trlist[[fn]]<-tr;telist[[fn]]<-te
      addqc(info,"FEATURE_OK",ntr=ntr,frac=frac,nte=nte,med=med,sdv=sdv)
    }
    rm(dt,h);invisible(gc(FALSE))
  }
  qc<-if(length(qclist))rbindlist(qclist,use.names=TRUE,fill=TRUE) else data.table()
  if(!length(trlist)){x<-empty_features();x$qc<-qc;return(x)}
  a<-do.call(cbind,trlist);b<-do.call(cbind,telist)
  if(is.null(dim(a))){a<-matrix(a,ncol=1);b<-matrix(b,ncol=1)}
  rownames(a)<-train_dt$ID;rownames(b)<-test_dt$ID
  if(anyDuplicated(colnames(a)))stop("Duplicate sanitized CpG feature names.")
  list(train=a,test=b,qc=qc)
}
features<-build_cpg_features(union_cpgs)
if(!identical(rownames(features$train),train_dt$ID)||!identical(rownames(features$test),test_dt$ID))
  stop("Feature sample order mismatch.")

fq<-copy(features$qc)
if(nrow(fq)){
  ann<-sel[,.(CpG,p_M1,p_M2,p_M3,selected_M1=as.integer(is.finite(p_M1)),
             selected_M2=as.integer(is.finite(p_M2)),selected_M3=as.integer(is.finite(p_M3)))]
  fq<-merge(fq,ann,by="CpG",all.x=TRUE,sort=FALSE)
}
fwrite(fq,feature_qc_file,sep="\t",quote=FALSE,na="NA")

# ----------------------------- prediction models -----------------------------
fit_covariate_reference <- function(x_train,y_train,x_test,foldid){
  tr<-as.data.frame(x_train);tr$y<-as.integer(y_train);te<-as.data.frame(x_test)
  oof<-rep(NA_real_,length(y_train))
  for(fd in sort(unique(foldid))){
    fi<-foldid!=fd;va<-foldid==fd
    fit<-tryCatch(suppressWarnings(glm(y~.,data=tr[fi,,drop=FALSE],family=binomial())),
                  error=function(e)e)
    if(inherits(fit,"error"))stop("Covariate inner fold failed: ",conditionMessage(fit))
    oof[va]<-as.numeric(predict(fit,newdata=tr[va,,drop=FALSE],type="response"))
  }
  fit<-tryCatch(suppressWarnings(glm(y~.,data=tr,family=binomial())),error=function(e)e)
  if(inherits(fit,"error"))stop("Final covariate model failed: ",conditionMessage(fit))
  testp<-as.numeric(predict(fit,newdata=te,type="response"))
  th<-choose_youden_threshold(y_train,oof)
  co<-data.table(feature=names(coef(fit)),coefficient=as.numeric(coef(fit)))
  co[,feature_type:=fifelse(feature=="(Intercept)","intercept","covariate")]
  list(oof_prob=oof,test_prob=testp,threshold=th,coefficients=co,lambda=NA_real_,
       n_features_input=ncol(x_train),
       n_features_nonzero=sum(co$feature!="(Intercept)"&is.finite(co$coefficient)&co$coefficient!=0),
       n_methylation_nonzero=0L,fit_status="FIT_OK")
}

fit_cpg_lasso <- function(xc_tr,xc_te,xv_tr,xv_te,y,foldid,seed){
  if(!ncol(xc_tr))return(NULL)
  xtr<-cbind(xc_tr,xv_tr);xte<-cbind(xc_te,xv_te)
  s<-apply(xtr,2,sd);keep<-is.finite(s)&s>variance_epsilon
  xtr<-xtr[,keep,drop=FALSE];xte<-xte[,keep,drop=FALSE]
  cpg<-grepl("^CPG_",colnames(xtr));cov<-grepl("^COV_",colnames(xtr))
  if(!any(cpg))return(NULL)
  if(!all(cpg|cov))warning("Some columns are neither CPG_ nor COV_.")
  pf<-ifelse(cpg,1,0);set.seed(seed)
  cv<-tryCatch(glmnet::cv.glmnet(
    x=xtr,y=as.integer(y),family="binomial",alpha=1,foldid=as.integer(foldid),
    type.measure="auc",standardize=TRUE,intercept=TRUE,penalty.factor=pf,
    keep=TRUE,grouped=TRUE,maxit=1000000,parallel=FALSE
  ),error=function(e)e)
  if(inherits(cv,"error"))stop("Single-CpG LASSO failed: ",conditionMessage(cv))
  lam<-as.numeric(cv$lambda.1se);li<-which.min(abs(cv$lambda-lam))
  oof<-as.numeric(cv$fit.preval[,li])
  if(any(oof<0|oof>1,na.rm=TRUE))oof<-plogis(oof)
  testp<-as.numeric(predict(cv,newx=xte,s="lambda.1se",type="response"))
  th<-choose_youden_threshold(y,oof)
  cm<-coef(cv,s="lambda.1se")
  co<-data.table(feature=rownames(cm),coefficient=as.numeric(cm))
  co[,feature_type:=fcase(feature=="(Intercept)","intercept",
                          grepl("^CPG_",feature),"methylation",
                          grepl("^COV_",feature),"covariate",default="other")]
  nz<-co[feature!="(Intercept)"&is.finite(coefficient)&coefficient!=0]
  list(oof_prob=oof,test_prob=testp,threshold=th,coefficients=co,lambda=lam,
       n_features_input=ncol(xtr),n_features_nonzero=nrow(nz),
       n_methylation_nonzero=nz[feature_type=="methylation",.N],fit_status="FIT_OK")
}

baseline<-fit_covariate_reference(covars$train,y_train,covars$test,foldid)
model_info<-list()
for(m in c("M1","M2","M3")){
  meth<-paste0("SingleCpG_",m);mn<-paste0(meth,"_LASSO");z<-model_sets[[m]]
  reqf<-safe_feature_name(z$CpG);usable<-intersect(reqf,colnames(features$train))
  xt<-features$train[,usable,drop=FALSE];xe<-features$test[,usable,drop=FALSE]
  reason<-NA_character_
  if(!nrow(z))reason<-"NO_SELECTED_CPGS" else if(!length(usable))reason<-"NO_USABLE_CPGS_AFTER_TRAIN_QC"
  if(!length(usable)){
    fit<-baseline;fallback<-TRUE
  }else{
    fit<-fit_cpg_lasso(xt,xe,covars$train,covars$test,y_train,foldid,
                       seed_base+10000L*match(m,c("M1","M2","M3"))+split_id)
    if(is.null(fit)){fit<-baseline;fallback<-TRUE;reason<-"NO_CPGS_AFTER_MODEL_MATRIX_QC"}else fallback<-FALSE
  }
  model_info[[m]]<-list(method=meth,model=mn,fit=fit,fallback=fallback,reason=reason,
                        nsel=nrow(z),nuse=length(usable))
  cat(meth,": selected=",nrow(z),", usable=",length(usable),", fallback=",fallback,
      if(!is.na(reason))paste0(", reason=",reason)else"","\n",sep="")
}

make_outputs <- function(fit,method,model,model_type,fallback,reason,nsel,nuse){
  mi<-calculate_metrics(y_train,fit$oof_prob,fit$threshold)
  mt<-calculate_metrics(y_test,fit$test_prob,fit$threshold)
  perf<-data.table(
    splitID=split_id,DMR_method=method,model=model,model_type=model_type,
    alpha=if(model_type=="LASSO")1 else NA_real_,lambda=fit$lambda,
    lambda_rule=if(model_type=="LASSO")"lambda.1se" else NA_character_,
    threshold=fit$threshold,threshold_source="family_grouped_inner_CV_Youden",
    inner_folds=fold_obj$n_folds,fit_status=fit$fit_status,model_fallback=fallback,
    fallback_reason=reason,n_train=nrow(train_dt),n_test=nrow(test_dt),
    n_train_FIDs=uniqueN(train_dt$FID),n_test_FIDs=uniqueN(test_dt$FID),
    n_DMR_selected_input=NA_integer_,n_DMR_usable=NA_integer_,
    n_CpG_selected_input=as.integer(nsel),n_CpG_usable=as.integer(nuse),
    n_methylation_features_input=as.integer(nuse),n_features_model_input=fit$n_features_input,
    n_features_nonzero=fit$n_features_nonzero,
    n_methylation_features_nonzero=fit$n_methylation_nonzero
  )
  perf<-cbind(perf,as.data.table(metrics_to_columns(mi,"inner_oof")),
             as.data.table(metrics_to_columns(mt,"test")))
  ptr<-data.table(splitID=split_id,DMR_method=method,model=model,ID=train_dt$ID,FID=train_dt$FID,
                  AA_only=y_train,prediction_set="train_inner_oof",inner_fold=foldid,
                  pred_prob=fit$oof_prob,threshold=fit$threshold,
                  pred_class=as.integer(fit$oof_prob>=fit$threshold),
                  model_fallback=fallback,fallback_reason=reason)
  pte<-data.table(splitID=split_id,DMR_method=method,model=model,ID=test_dt$ID,FID=test_dt$FID,
                  AA_only=y_test,prediction_set="test",inner_fold=NA_integer_,
                  pred_prob=fit$test_prob,threshold=fit$threshold,
                  pred_class=as.integer(fit$test_prob>=fit$threshold),
                  model_fallback=fallback,fallback_reason=reason)
  co<-copy(fit$coefficients)
  co[,`:=`(splitID=split_id,DMR_method=method,model=model,
           selected_nonzero=as.integer(is.finite(coefficient)&coefficient!=0),
           model_fallback=fallback,fallback_reason=reason)]
  list(performance=perf,predictions=rbindlist(list(ptr,pte)),coefficients=co)
}

oo<-lapply(c("M1","M2","M3"),function(m){
  z<-model_info[[m]]
  make_outputs(z$fit,z$method,z$model,"LASSO",z$fallback,z$reason,z$nsel,z$nuse)
})
baseout<-make_outputs(baseline,"Covariates_only","Covariate_only","unpenalized_logistic",
                      FALSE,NA_character_,0L,0L)

performance_out<-rbindlist(c(lapply(oo,`[[`,"performance"),list(baseout$performance)),fill=TRUE)
prediction_out<-rbindlist(c(lapply(oo,`[[`,"predictions"),list(baseout$predictions)),fill=TRUE)
coefficient_out<-rbindlist(c(lapply(oo,`[[`,"coefficients"),list(baseout$coefficients)),fill=TRUE)
fwrite(performance_out,performance_file,sep="\t",quote=FALSE,na="NA")
fwrite(prediction_out,prediction_file,sep="\t",quote=FALSE,na="NA")
fwrite(coefficient_out,coefficient_file,sep="\t",quote=FALSE,na="NA")

split_qc<-rbindlist(lapply(c("M1","M2","M3"),function(m){
  z<-model_info[[m]]
  data.table(
    splitID=split_id,DMR_method=z$method,selected_cpg_file=selected_cpg_file,
    n_selected_CpGs=z$nsel,n_usable_CpGs=z$nuse,n_feature_failures=z$nsel-z$nuse,
    n_train=nrow(train_dt),n_test=nrow(test_dt),n_train_FIDs=uniqueN(train_dt$FID),
    n_test_FIDs=uniqueN(test_dt$FID),n_train_cases=sum(y_train==1L),
    n_train_controls=sum(y_train==0L),n_test_cases=sum(y_test==1L),
    n_test_controls=sum(y_test==0L),inner_folds=fold_obj$n_folds,
    family_overlap_count=fam_overlap,model_fallback=z$fallback,fallback_reason=z$reason,
    min_train_cpg_observed_fraction=min_obs,variance_epsilon=variance_epsilon,
    elapsed_minutes=as.numeric(difftime(Sys.time(),run_start,units="mins")),
    completed_at=format(Sys.time(),tz="America/Toronto",usetz=TRUE),job_status="COMPLETED")
}),fill=TRUE)
fwrite(split_qc,split_qc_file,sep="\t",quote=FALSE,na="NA")

cat("\nCompleted split",split_id,"\n")
print(performance_out[,.(splitID,DMR_method,model,model_fallback,fallback_reason,
                         n_CpG_selected_input,n_CpG_usable,AUROC_test,PR_AUC_test,Brier_test)])
