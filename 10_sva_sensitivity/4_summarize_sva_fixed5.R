#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", unset="~/scratch/UQAC/meth"))
out_dir <- path.expand(Sys.getenv(
  "SVA_SENS_OUT",
  unset=file.path(PATH_wk,"results","15_revision","7_sva_threshold_sensitivity")
))
fixed_dir <- file.path(out_dir,"fixed_nsv5")
fig_dir <- file.path(fixed_dir,"figures")
tab_dir <- file.path(fixed_dir,"tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

labels <- c("10k","25k","50k")
vfilters <- c(10000L,25000L,50000L)
files <- setNames(file.path(fixed_dir, sprintf("sva_fixed5_vfilter_%05d.rds",vfilters)), labels)
if (any(!file.exists(files))) stop("Missing one or more fixed-n.sv result files")

res <- lapply(files, readRDS)
ids <- lapply(res, function(x) as.character(x$subject_ids))
if (!all(vapply(ids[-1], identical, logical(1), ids[[1]]))) stop("Subject order mismatch")
svs <- lapply(res, function(x) as.matrix(x$sv))
if (!all(vapply(svs, ncol, integer(1)) == 5L)) stop("All solutions must contain 5 SVs")

orthobasis <- function(X) {
  X <- scale(as.matrix(X), center=TRUE, scale=FALSE)
  q <- qr(X)
  qr.Q(q)[,seq_len(q$rank),drop=FALSE]
}

all_perms <- function(x) {
  if (length(x)==1L) return(matrix(x,nrow=1L))
  do.call(rbind,lapply(seq_along(x), function(i) cbind(x[i], all_perms(x[-i]))))
}
P5 <- all_perms(1:5)

optimal_match <- function(C) {
  scores <- apply(P5,1,function(p) sum(C[cbind(1:5,p)]))
  p <- P5[which.max(scores),]
  data.table(
    SV_A=rownames(C),
    SV_B=colnames(C)[p],
    abs_correlation=C[cbind(1:5,p)]
  )
}

ARI <- function(x,y) {
  tab <- table(x,y)
  n <- sum(tab)
  c2 <- function(z) z*(z-1)/2
  a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); c <- sum(c2(colSums(tab))); t <- c2(n)
  if (t==0) return(NA_real_)
  e <- b*c/t
  m <- 0.5*(b+c)
  if (abs(m-e)<.Machine$double.eps) return(1)
  (a-e)/(m-e)
}

pair_metrics <- function(A,B,la,lb) {
  C <- abs(cor(A,B,use="pairwise.complete.obs"))
  QA <- orthobasis(A); QB <- orthobasis(B)

  pc <- svd(crossprod(QA,QB),nu=0,nv=0)$d
  pc <- pmin(pmax(pc,0),1)

  PA <- tcrossprod(QA); PB <- tcrossprod(QB)
  proj_sim <- sum(PA*PB)/sqrt(sum(PA*PA)*sum(PB*PB))

  dA <- as.vector(dist(QA)); dB <- as.vector(dist(QB))
  dist_rho <- cor(dA,dB,method="spearman")

  hcA <- hclust(dist(QA),method="ward.D2")
  hcB <- hclust(dist(QB),method="ward.D2")
  coph_rho <- cor(as.vector(cophenetic(hcA)),as.vector(cophenetic(hcB)),method="spearman")

  mm <- optimal_match(C)

  list(
    C=C, QA=QA, QB=QB, hcA=hcA, hcB=hcB, match=mm,
    summary=data.table(
      comparison=paste(la,"vs",lb),
      mean_matched_abs_SV_correlation=mean(mm$abs_correlation),
      min_matched_abs_SV_correlation=min(mm$abs_correlation),
      mean_principal_angle_cosine=mean(pc),
      min_principal_angle_cosine=min(pc),
      projection_matrix_similarity=proj_sim,
      subject_distance_spearman=dist_rho,
      clustering_cophenetic_spearman=coph_rho
    )
  )
}

overall <- rbindlist(lapply(labels, function(nm) {
  x <- res[[nm]]
  data.table(
    threshold=nm,
    vfilter=x$vfilter,
    fixed_n_sv=x$fixed_n_sv,
    n_subjects=nrow(x$sv),
    sva_version=x$sva_version,
    seed=x$seed
  )
}))
fwrite(overall,file.path(tab_dir,"Table_fixed5_overall_summary.csv"))

pairs <- list(c("10k","25k"),c("10k","50k"),c("25k","50k"))
pair_results <- list()
sum_list <- list()
match_list <- list()
cluster_list <- list()

for (i in seq_along(pairs)) {
  a <- pairs[[i]][1]; b <- pairs[[i]][2]
  cmp <- pair_metrics(svs[[a]],svs[[b]],a,b)
  key <- paste(a,b,sep="_")
  pair_results[[key]] <- cmp
  sum_list[[i]] <- cmp$summary

  mm <- copy(cmp$match)
  mm[,comparison:=paste(a,"vs",b)]
  mm[,SV_A:=paste0(a,"_",SV_A)]
  mm[,SV_B:=paste0(b,"_",SV_B)]
  setcolorder(mm,c("comparison","SV_A","SV_B","abs_correlation"))
  match_list[[i]] <- mm

  cdt <- as.data.table(cmp$C,keep.rownames="SV_A")
  fwrite(cdt,file.path(tab_dir,sprintf("Cor_fixed5_%s_vs_%s.csv",a,b)))

  cluster_list[[i]] <- rbindlist(lapply(2:6,function(k) {
    data.table(
      comparison=paste(a,"vs",b),
      k_clusters=k,
      adjusted_rand_index=ARI(cutree(cmp$hcA,k),cutree(cmp$hcB,k))
    )
  }))
}

pair_tab <- rbindlist(sum_list)
match_tab <- rbindlist(match_list)
cluster_tab <- rbindlist(cluster_list)

fwrite(pair_tab,file.path(tab_dir,"Table_fixed5_pairwise_stability.csv"))
fwrite(match_tab,file.path(tab_dir,"Table_fixed5_matched_SV_correlations.csv"))
fwrite(cluster_tab,file.path(tab_dir,"Table_fixed5_clustering_stability.csv"))

# Figure 1: pairwise absolute SV-correlation heatmaps
png(file.path(fig_dir,"Figure_fixed5_SV_correlations.png"),width=3000,height=1000,res=220)
par(mfrow=c(1,3),mar=c(5,5,4,2))
for (pp in pairs) {
  a <- pp[1]; b <- pp[2]
  C <- pair_results[[paste(a,b,sep="_")]]$C
  image(1:5,1:5,t(C[5:1,,drop=FALSE]),axes=FALSE,xlab=paste0(b," SVs"),
        ylab=paste0(a," SVs"),main=paste(a,"vs",b),zlim=c(0,1))
  axis(1,at=1:5,labels=paste0("SV",1:5))
  axis(2,at=1:5,labels=paste0("SV",5:1),las=1)
  for (r in 1:5) for (cc in 1:5)
    text(cc,6-r,sprintf("%.2f",C[r,cc]),cex=.8)
}
dev.off()

# Figure 2: pairwise subject-distance concordance
png(file.path(fig_dir,"Figure_fixed5_subject_distance.png"),width=3000,height=1000,res=220)
par(mfrow=c(1,3),mar=c(5,5,4,2))
for (pp in pairs) {
  a <- pp[1]; b <- pp[2]
  cmp <- pair_results[[paste(a,b,sep="_")]]
  dA <- as.vector(dist(cmp$QA)); dB <- as.vector(dist(cmp$QB))
  rho <- cor(dA,dB,method="spearman")
  plot(dA,dB,pch=16,cex=.25,xlab=paste0(a," subject distance"),
       ylab=paste0(b," subject distance"),
       main=paste0(a," vs ",b,"\nSpearman rho=",sprintf("%.3f",rho)))
  abline(lm(dB~dA),lwd=1.5)
}
dev.off()

# Figure 3: summary of rotation-invariant stability
metrics <- c("mean_matched_abs_SV_correlation","mean_principal_angle_cosine",
             "projection_matrix_similarity","subject_distance_spearman",
             "clustering_cophenetic_spearman")
M <- as.matrix(pair_tab[,..metrics])
rownames(M) <- pair_tab$comparison

png(file.path(fig_dir,"Figure_fixed5_subspace_similarity.png"),width=2400,height=1500,res=220)
par(mar=c(8,5,4,2))
barplot(t(M),beside=TRUE,ylim=c(0,1.05),las=2,ylab="Similarity",
        main="Stability of fixed five-dimensional SVA solutions")
legend("bottomleft",
       legend=c("Matched SV |r|","Mean principal-angle cosine","Projection similarity",
                "Subject-distance rho","Cophenetic rho"),
       bty="n",cex=.8)
abline(h=c(.8,.9,1),lty=c(3,3,2))
dev.off()

# Figure 4: descriptive dendrograms
png(file.path(fig_dir,"Figure_fixed5_dendrograms.png"),width=3000,height=1100,res=220)
par(mfrow=c(1,3),mar=c(3,4,4,1))
for (nm in labels) {
  hc <- hclust(dist(orthobasis(svs[[nm]])),method="ward.D2")
  plot(hc,labels=FALSE,hang=-1,main=paste0(nm," CpG vfilter"),xlab="Subjects",sub="")
}
dev.off()

# Human-readable summary
sink(file.path(fixed_dir,"SVA_fixed5_results_summary.txt"))
cat("Fixed n.sv=5 SVA threshold sensitivity analysis\n")
cat("=============================================\n\n")
cat("Overall design\n"); print(overall); cat("\n")
cat("Pairwise stability\n"); print(pair_tab); cat("\n")
cat("Matched SV correlations\n"); print(match_tab); cat("\n")
cat("Clustering stability (ARI)\n"); print(cluster_tab); cat("\n")
cat("Interpretation: individual SV signs/order/rotation are not unique; ",
    "the subspace, projection, subject-distance and cophenetic metrics are therefore ",
    "stronger evidence of latent-space stability than SV1-to-SV1 comparison.\n",sep="")
sink()

print(pair_tab)
cat("Tables:",tab_dir,"\n")
cat("Figures:",fig_dir,"\n")
