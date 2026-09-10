#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
  library(ggplot2)
})
data.table::setDTthreads(1L)

PATH_wk <- path.expand(Sys.getenv("PATH_WK","~/scratch/UQAC/meth/"))
pred_root <- file.path(PATH_wk,"results/15_revision/5_cv/5_prediction")
out_root <- path.expand(Sys.getenv(
  "ALL_METHOD_COMPARISON_OUTPUT",
  file.path(pred_root,"7_all_method_comparison")
))
tab_dir <- file.path(out_root,"tables")
fig_dir <- file.path(out_root,"figures")
dat_dir <- file.path(out_root,"data")
for(d in c(out_root,tab_dir,fig_dir,dat_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

expected_splits <- as.integer(Sys.getenv("EXPECTED_SPLITS","100"))
overwrite <- Sys.getenv("OVERWRITE","0") %in% c("1","TRUE","true","T","yes","YES")

roots <- list(
  GAM_DMR=file.path(pred_root,"1_mgcv_primary"),
  SOMNiBUS=file.path(pred_root,"2_somnibus_primary"),
  DMRcate=file.path(pred_root,"3_dmrcate_primary"),
  BSmooth=file.path(pred_root,"4_bsmooth_primary"),
  SingleCpG=file.path(pred_root,"6_single_cpg_primary")
)
roots$GAM_DMR <- path.expand(Sys.getenv("MGCV_PREDICTION_OUTPUT",roots$GAM_DMR))
roots$SOMNiBUS <- path.expand(Sys.getenv("SOMNIBUS_PREDICTION_OUTPUT",roots$SOMNiBUS))
roots$DMRcate <- path.expand(Sys.getenv("DMRCATE_PREDICTION_OUTPUT",roots$DMRcate))
roots$BSmooth <- path.expand(Sys.getenv("BSMOOTH_PREDICTION_OUTPUT",roots$BSmooth))
roots$SingleCpG <- path.expand(Sys.getenv("SINGLE_CPG_PREDICTION_OUTPUT",roots$SingleCpG))

method_order <- c(
  "GAM-DMR","SOMNiBUS","DMRcate","BSmooth",
  "Univariate method M1","Univariate method M2","Univariate method M3","Covariates only"
)
method_palette <- c(
  "GAM-DMR"="#D55E00","SOMNiBUS"="#0072B2","DMRcate"="#009E73",
  "BSmooth"="#CC79A7","Univariate method M1"="#E69F00",
  "Univariate method M2"="#56B4E9","Univariate method M3"="#000000",
  "Covariates only"="#777777"
)
method_linetype <- c(
  "GAM-DMR"="solid","SOMNiBUS"="solid","DMRcate"="solid","BSmooth"="solid",
  "Univariate method M1"="longdash","Univariate method M2"="dashed",
  "Univariate method M3"="dotdash","Covariates only"="dotted"
)
method_linewidth <- c(
  "GAM-DMR"=1.45,"SOMNiBUS"=1,"DMRcate"=1,"BSmooth"=1,
  "Univariate method M1"=1,"Univariate method M2"=1,"Univariate method M3"=1,"Covariates only"=.95
)

find_collected <- function(root,type,preferred=NULL){
  d<-file.path(root,"collected")
  if(!dir.exists(d)) stop("Missing collected directory: ",d)
  if(!is.null(preferred)){
    p<-file.path(d,preferred)
    if(file.exists(p)) return(p)
  }
  pat<-switch(type,
    performance="primary_performance_all_splits\\.tsv$",
    predictions="primary_predictions_all_splits\\.tsv$",
    stop("bad type"))
  h<-list.files(d,pattern=pat,full.names=TRUE)
  if(length(h)!=1L) stop("Expected one ",type," file under ",d,"; found ",length(h))
  h
}
safe_read <- function(f,label){
  if(!file.exists(f)) stop(label," missing: ",f)
  z<-fread(f)
  if(!nrow(z)) stop(label," empty: ",f)
  z
}
safe_roc_auc <- function(y,p){
  k<-is.finite(y)&is.finite(p);y<-as.integer(y[k]);p<-as.numeric(p[k])
  if(length(unique(y))<2) return(NA_real_)
  z<-tryCatch(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE),error=function(e)NULL)
  if(is.null(z)) NA_real_ else as.numeric(pROC::auc(z))
}
safe_pr_auc <- function(y,p){
  k<-is.finite(y)&is.finite(p);y<-as.integer(y[k]);p<-as.numeric(p[k])
  if(length(unique(y))<2) return(NA_real_)
  z<-tryCatch(PRROC::pr.curve(scores.class0=p[y==1L],scores.class1=p[y==0L],curve=FALSE),
              error=function(e)NULL)
  if(is.null(z)) NA_real_ else as.numeric(z$auc.integral)
}
aggregate_prob <- function(y,p){
  k<-is.finite(y)&is.finite(p);y<-as.integer(y[k]);p<-as.numeric(p[k])
  pc<-pmin(pmax(p,1e-8),1-1e-8);lp<-qlogis(pc);ci<-cs<-NA_real_
  if(length(unique(y))>=2 && is.finite(sd(lp)) && sd(lp)>0){
    fi<-tryCatch(suppressWarnings(glm(y~1+offset(lp),family=binomial())),error=function(e)NULL)
    fs<-tryCatch(suppressWarnings(glm(y~lp,family=binomial())),error=function(e)NULL)
    if(!is.null(fi))ci<-unname(coef(fi)[1])
    if(!is.null(fs)&&length(coef(fs))>=2)cs<-unname(coef(fs)[2])
  }
  pr<-safe_pr_auc(y,p)
  data.table(n_subjects=length(y),n_cases=sum(y==1),n_controls=sum(y==0),
             prevalence=mean(y),AUROC=safe_roc_auc(y,p),PR_AUC=pr,
             PR_baseline=mean(y),PR_AUC_lift=pr-mean(y),
             Brier=mean((p-y)^2),
             log_loss=-mean(y*log(pc)+(1-y)*log(1-pc)),
             calibration_intercept=ci,calibration_slope=cs)
}
corrected_cv <- function(d,ratio){
  k<-is.finite(d)&is.finite(ratio);d<-d[k];ratio<-ratio[k];n<-length(d)
  if(n<3)return(data.table(n_splits=n,mean_effect=if(n)mean(d)else NA_real_,
                           median_effect=if(n)median(d)else NA_real_,
                           corrected_SE=NA_real_,corrected_CI_low=NA_real_,
                           corrected_CI_high=NA_real_,corrected_t=NA_real_,
                           corrected_df=if(n)n-1L else NA_integer_,
                           corrected_pvalue=NA_real_,
                           mean_test_train_ratio=if(n)mean(ratio)else NA_real_,
                           win_rate_first_method=if(n)mean(d>0)else NA_real_))
  se<-sqrt((1/n+mean(ratio))*var(d));df<-n-1L;md<-mean(d)
  if(!is.finite(se)||se==0){
    tt<-pv<-NA_real_;lo<-hi<-md
  }else{
    tt<-md/se;pv<-2*pt(-abs(tt),df=df);cr<-qt(.975,df=df);lo<-md-cr*se;hi<-md+cr*se
  }
  data.table(n_splits=n,mean_effect=md,median_effect=median(d),corrected_SE=se,
             corrected_CI_low=lo,corrected_CI_high=hi,corrected_t=tt,
             corrected_df=df,corrected_pvalue=pv,mean_test_train_ratio=mean(ratio),
             win_rate_first_method=mean(d>0))
}
summ <- function(x){
  x<-as.numeric(x);x<-x[is.finite(x)]
  if(!length(x))return(c(mean=NA,sd=NA,median=NA,q025=NA,q975=NA))
  c(mean=mean(x),sd=sd(x),median=median(x),q025=quantile(x,.025,names=FALSE),
    q975=quantile(x,.975,names=FALSE))
}

files <- list(
  GAM_perf=find_collected(roots$GAM_DMR,"performance","mgcv_primary_performance_all_splits.tsv"),
  GAM_pred=find_collected(roots$GAM_DMR,"predictions","mgcv_primary_predictions_all_splits.tsv"),
  SOM_perf=find_collected(roots$SOMNiBUS,"performance","somnibus_primary_performance_all_splits.tsv"),
  SOM_pred=find_collected(roots$SOMNiBUS,"predictions","somnibus_primary_predictions_all_splits.tsv"),
  DMR_perf=find_collected(roots$DMRcate,"performance","dmrcate_primary_performance_all_splits.tsv"),
  DMR_pred=find_collected(roots$DMRcate,"predictions","dmrcate_primary_predictions_all_splits.tsv"),
  BS_perf=find_collected(roots$BSmooth,"performance","bsmooth_primary_performance_all_splits.tsv"),
  BS_pred=find_collected(roots$BSmooth,"predictions","bsmooth_primary_predictions_all_splits.tsv"),
  CPG_perf=find_collected(roots$SingleCpG,"performance","single_cpg_primary_performance_all_splits.tsv"),
  CPG_pred=find_collected(roots$SingleCpG,"predictions","single_cpg_primary_predictions_all_splits.tsv")
)
print(files)

gp<-safe_read(files$GAM_perf,"GAM perf"); gpr<-safe_read(files$GAM_pred,"GAM pred")
sp<-safe_read(files$SOM_perf,"SOM perf"); spr<-safe_read(files$SOM_pred,"SOM pred")
dp<-safe_read(files$DMR_perf,"DMRcate perf"); dpr<-safe_read(files$DMR_pred,"DMRcate pred")
bp<-safe_read(files$BS_perf,"BSmooth perf"); bpr<-safe_read(files$BS_pred,"BSmooth pred")
cp<-safe_read(files$CPG_perf,"CpG perf"); cpr<-safe_read(files$CPG_pred,"CpG pred")

pick <- function(perf,pred,target_model,display,source){
  a<-perf[as.character(model)==target_model]
  b<-pred[as.character(model)==target_model]
  if(!nrow(a)||!nrow(b))stop(source,": missing model ",target_model)
  a[,`:=`(benchmark_method=display,benchmark_source=source)]
  b[,`:=`(benchmark_method=display,benchmark_source=source)]
  list(performance=a,predictions=b)
}
parts<-list(
  pick(gp,gpr,"MGCV_PC12_LASSO","GAM-DMR","MGCV/GAM-DMR"),
  pick(sp,spr,"SOMNiBUS_PC12_LASSO","SOMNiBUS","SOMNiBUS"),
  pick(dp,dpr,"DMRcate_PC12_LASSO","DMRcate","DMRcate"),
  pick(bp,bpr,"BSmooth_PC12_LASSO","BSmooth","BSmooth"),
  pick(cp,cpr,"SingleCpG_M1_LASSO","Univariate method M1","Single-CpG"),
  pick(cp,cpr,"SingleCpG_M2_LASSO","Univariate method M2","Single-CpG"),
  pick(cp,cpr,"SingleCpG_M3_LASSO","Univariate method M3","Single-CpG"),
  pick(gp,gpr,"Covariate_only","Covariates only","GAM baseline")
)
perf<-rbindlist(lapply(parts,`[[`,"performance"),use.names=TRUE,fill=TRUE)
pred<-rbindlist(lapply(parts,`[[`,"predictions"),use.names=TRUE,fill=TRUE)
perf[,benchmark_method:=factor(benchmark_method,levels=method_order)]
pred[,benchmark_method:=factor(benchmark_method,levels=method_order)]

q<-perf[,.(n_splits=uniqueN(splitID),n_rows=.N),by=benchmark_method]
if(any(q$n_splits!=expected_splits)||any(q$n_rows!=expected_splits)){print(q);stop("Performance split QC failed.")}
q2<-pred[,.(n_splits=uniqueN(splitID)),by=benchmark_method]
if(any(q2$n_splits!=expected_splits)){print(q2);stop("Prediction split QC failed.")}
sizes<-perf[,.(ntr=uniqueN(n_train),nte=uniqueN(n_test)),by=splitID]
if(any(sizes$ntr!=1|sizes$nte!=1))stop("Outer split sizes differ across methods.")

fwrite(perf,file.path(dat_dir,"all_methods_primary_performance_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(pred,file.path(dat_dir,"all_methods_primary_predictions_all_splits.tsv"),sep="\t",quote=FALSE,na="NA")

metrics<-intersect(c("AUROC_test","PR_AUC_test","PR_AUC_lift_test","Brier_test","log_loss_test",
                     "calibration_intercept_test","calibration_slope_test","sensitivity_test",
                     "specificity_test","balanced_accuracy_test","PPV_test","NPV_test"),names(perf))
rows<-list()
for(m in method_order){
  d<-perf[benchmark_method==m]
  r<-data.table(benchmark_method=m,n_splits=uniqueN(d$splitID),
                n_fallback_splits=if("model_fallback"%in%names(d))sum(d$model_fallback%in%TRUE,na.rm=TRUE)else NA_integer_)
  for(v in metrics){
    s<-summ(d[[v]])
    for(nm in names(s))r[,(paste0(v,"_",nm)):=as.numeric(s[nm])]
  }
  rows[[length(rows)+1L]]<-r
}
ps<-rbindlist(rows,fill=TRUE)
ps[,benchmark_method:=factor(benchmark_method,levels=method_order)]
setorder(ps,benchmark_method)
fwrite(ps,file.path(tab_dir,"Table1_all_methods_primary_performance.tsv"),sep="\t",quote=FALSE,na="NA")
fmt<-ps[,.(Method=as.character(benchmark_method),
           `Test AUROC, mean (SD)`=sprintf("%.3f (%.3f)",AUROC_test_mean,AUROC_test_sd),
           `Test PR-AUC, mean (SD)`=sprintf("%.3f (%.3f)",PR_AUC_test_mean,PR_AUC_test_sd),
           `Test Brier, mean (SD)`=sprintf("%.3f (%.3f)",Brier_test_mean,Brier_test_sd),
           `Fallback splits`=n_fallback_splits)]
fwrite(fmt,file.path(tab_dir,"Table1_all_methods_primary_performance_formatted.tsv"),
       sep="\t",quote=FALSE,na="NA")

# Corrected repeated-CV pairwise comparisons. Effects are oriented so positive favors method1.
pw<-list();kk<-0L
for(metric in c("AUROC_test","PR_AUC_test","Brier_test")){
  for(i in 1:(length(method_order)-1)){
    for(j in (i+1):length(method_order)){
      m1<-method_order[i];m2<-method_order[j]
      a<-perf[benchmark_method==m1,.(splitID,v1=get(metric),n_train,n_test)]
      b<-perf[benchmark_method==m2,.(splitID,v2=get(metric))]
      z<-merge(a,b,by="splitID")
      z[,effect:=if(metric=="Brier_test")v2-v1 else v1-v2]
      ans<-corrected_cv(z$effect,z$n_test/z$n_train)
      ans[,`:=`(metric=metric,method1=m1,method2=m2,
                comparison=paste0(m1," vs ",m2),
                effect_orientation="positive_favors_method1")]
      kk<-kk+1L;pw[[kk]]<-ans
    }
  }
}
pw<-rbindlist(pw,fill=TRUE)
pw[,corrected_pvalue_Holm:=p.adjust(corrected_pvalue,method="holm"),by=metric]
setcolorder(pw,c("metric","method1","method2","comparison","effect_orientation",
                 "n_splits","mean_effect","median_effect","corrected_SE",
                 "corrected_CI_low","corrected_CI_high","corrected_t","corrected_df",
                 "corrected_pvalue","corrected_pvalue_Holm","mean_test_train_ratio",
                 "win_rate_first_method"))
fwrite(pw,file.path(tab_dir,"Table_S1_all_methods_pairwise_corrected_repeated_CV.tsv"),
       sep="\t",quote=FALSE,na="NA")
fwrite(pw[method1=="GAM-DMR"],
       file.path(tab_dir,"Table2_GAM_DMR_vs_all_alternatives_corrected_CV.tsv"),
       sep="\t",quote=FALSE,na="NA")

# Subject-average held-out predictions, same approach as previous regional comparison.
pt<-pred[prediction_set=="test"]
pt[,benchmark_method:=as.character(benchmark_method)]
oc<-pt[,.(n_outcomes=uniqueN(AA_only)),by=.(benchmark_method,ID,FID)]
if(any(oc$n_outcomes!=1L))stop("Inconsistent outcomes across held-out appearances.")
sa<-pt[,.(AA_only=unique(AA_only)[1L],
          mean_pred_prob=mean(pred_prob,na.rm=TRUE),
          median_pred_prob=median(pred_prob,na.rm=TRUE),
          sd_pred_prob=sd(pred_prob,na.rm=TRUE),
          n_test_appearances=.N,
          mean_split_threshold=mean(threshold,na.rm=TRUE)),
       by=.(benchmark_method,ID,FID)]
sa[,benchmark_method:=factor(benchmark_method,levels=method_order)]
setorder(sa,benchmark_method,ID)
fwrite(sa,file.path(dat_dir,"all_methods_subject_average_test_predictions.tsv"),
       sep="\t",quote=FALSE,na="NA")
sap<-sa[,aggregate_prob(AA_only,mean_pred_prob),by=benchmark_method]
sap[,benchmark_method:=factor(benchmark_method,levels=method_order)]
setorder(sap,benchmark_method)
fwrite(sap,file.path(tab_dir,"Table_S2_all_methods_subject_average_performance.tsv"),
       sep="\t",quote=FALSE,na="NA")

# ROC and PR curve coordinates.
rl<-list();pl<-list()
for(m in method_order){
  d<-sa[benchmark_method==m]
  ro<-pROC::roc(d$AA_only,d$mean_pred_prob,levels=c(0,1),direction="<",quiet=TRUE)
  au<-as.numeric(pROC::auc(ro))
  rl[[m]]<-data.table(benchmark_method=m,FPR=1-ro$specificities,TPR=ro$sensitivities,AUROC=au)
  po<-PRROC::pr.curve(scores.class0=d$mean_pred_prob[d$AA_only==1L],
                      scores.class1=d$mean_pred_prob[d$AA_only==0L],curve=TRUE)
  pl[[m]]<-data.table(benchmark_method=m,Recall=po$curve[,1],Precision=po$curve[,2],
                      AUPRC=as.numeric(po$auc.integral))
}
roc<-rbindlist(rl);pr<-rbindlist(pl)
fwrite(roc,file.path(dat_dir,"all_methods_subject_average_ROC_coordinates.tsv"),sep="\t",na="NA")
fwrite(pr,file.path(dat_dir,"all_methods_subject_average_PR_coordinates.tsv"),sep="\t",na="NA")

# EXACTLY TWO DIGITS in figure legends.
rlabels<-unique(roc[,.(benchmark_method,AUROC)])
rlabels[,curve_label:=paste0(benchmark_method," (AUROC=",sprintf("%.2f",AUROC),")")]
plabels<-unique(pr[,.(benchmark_method,AUPRC)])
plabels[,curve_label:=paste0(benchmark_method," (AUPRC=",sprintf("%.2f",AUPRC),")")]
rmap<-setNames(rlabels$curve_label,rlabels$benchmark_method)
pmap<-setNames(plabels$curve_label,plabels$benchmark_method)
roc[,curve_label:=rmap[benchmark_method]]
pr[,curve_label:=pmap[benchmark_method]]
rorder<-unname(rmap[method_order]);porder<-unname(pmap[method_order])
roc[,curve_label:=factor(curve_label,levels=rorder)]
pr[,curve_label:=factor(curve_label,levels=porder)]
rpal<-setNames(unname(method_palette[method_order]),rorder)
ppal<-setNames(unname(method_palette[method_order]),porder)
rlt<-setNames(unname(method_linetype[method_order]),rorder)
plt<-setNames(unname(method_linetype[method_order]),porder)
rlw<-setNames(unname(method_linewidth[method_order]),rorder)
plw<-setNames(unname(method_linewidth[method_order]),porder)

p_roc<-ggplot(roc,aes(FPR,TPR,color=curve_label,linetype=curve_label,linewidth=curve_label))+
  geom_line()+geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=.6,color="grey55")+
  scale_color_manual(values=rpal,drop=FALSE)+scale_linetype_manual(values=rlt,drop=FALSE)+
  scale_linewidth_manual(values=rlw,drop=FALSE)+
  coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_bw(base_size=13)+
  theme(legend.position="right",legend.title=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Subject-averaged out-of-fold ROC curves",x="False positive rate",y="True positive rate")

uo<-unique(sa[,.(ID,AA_only)]);prevalence<-mean(uo$AA_only)
p_pr<-ggplot(pr,aes(Recall,Precision,color=curve_label,linetype=curve_label,linewidth=curve_label))+
  geom_line()+geom_hline(yintercept=prevalence,linetype="dashed",linewidth=.6,color="grey55")+
  scale_color_manual(values=ppal,drop=FALSE)+scale_linetype_manual(values=plt,drop=FALSE)+
  scale_linewidth_manual(values=plw,drop=FALSE)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_bw(base_size=13)+
  theme(legend.position="right",legend.title=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Subject-averaged out-of-fold precision-recall curves",x="Recall",y="Precision")

ggsave(file.path(fig_dir,"Figure4A_all_methods_subject_average_ROC.pdf"),p_roc,width=9.5,height=6.5)
ggsave(file.path(fig_dir,"Figure4A_all_methods_subject_average_ROC.png"),p_roc,width=9.5,height=6.5,dpi=300)
ggsave(file.path(fig_dir,"Figure4B_all_methods_subject_average_PR.pdf"),p_pr,width=9.5,height=6.5)
ggsave(file.path(fig_dir,"Figure4B_all_methods_subject_average_PR.png"),p_pr,width=9.5,height=6.5,dpi=300)

# Split-performance distributions across all methods.
pv<-rbindlist(list(
  perf[,.(splitID,benchmark_method,metric="AUROC",value=AUROC_test)],
  perf[,.(splitID,benchmark_method,metric="PR-AUC",value=PR_AUC_test)],
  perf[,.(splitID,benchmark_method,metric="Brier",value=Brier_test)]
))
pv[,benchmark_method:=factor(benchmark_method,levels=method_order)]
pv[,metric:=factor(metric,levels=c("AUROC","PR-AUC","Brier"))]
means<-pv[,.(value=mean(value,na.rm=TRUE)),by=.(benchmark_method,metric)]
pdist<-ggplot(pv,aes(benchmark_method,value,fill=benchmark_method))+
  geom_violin(trim=FALSE,alpha=.45,linewidth=.4)+
  geom_boxplot(width=.14,outlier.shape=NA,linewidth=.4)+
  geom_point(data=means,shape=23,size=2.4,fill="white",color="black")+
  facet_wrap(~metric,scales="free_y",nrow=1)+scale_fill_manual(values=method_palette,drop=FALSE)+
  theme_bw(base_size=12)+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1),
                               panel.grid.minor=element_blank())+
  labs(title="Predictive performance across repeated family-aware test splits",x=NULL,y=NULL)
ggsave(file.path(fig_dir,"Figure1_all_methods_split_performance_distributions.pdf"),pdist,width=13,height=5.5)
ggsave(file.path(fig_dir,"Figure1_all_methods_split_performance_distributions.png"),pdist,width=13,height=5.5,dpi=300)

qc<-data.table(expected_splits=expected_splits,n_methods=length(method_order),
               methods=paste(method_order,collapse=";"),
               performance_rows=nrow(perf),prediction_rows=nrow(pred),
               subject_average_rows=nrow(sa),subject_average_unique_subjects=uniqueN(sa$ID),
               PR_baseline=prevalence,
               completed_at=format(Sys.time(),tz="America/Toronto",usetz=TRUE),
               job_status="COMPLETED")
fwrite(qc,file.path(out_root,"all_method_comparison_qc.tsv"),sep="\t",quote=FALSE,na="NA")

cat("\nAll-method comparison completed:",out_root,"\n")
print(ps[,.(benchmark_method,n_splits,n_fallback_splits,
            mean_AUROC=AUROC_test_mean,mean_PR_AUC=PR_AUC_test_mean,mean_Brier=Brier_test_mean)])
cat("\nSubject-average performance:\n")
print(sap[,.(benchmark_method,AUROC,PR_AUC,Brier)])
cat("\nROC legend values (2 decimals):\n");print(rlabels)
cat("\nPR legend values (2 decimals):\n");print(plabels)
