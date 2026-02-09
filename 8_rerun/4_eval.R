

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(qqman)

library(VennDiagram)

#before rerun, the first stage results
AA_simple_results <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/rerun_simple.txt"), header=FALSE)
colnames(AA_simple_results) <- c("Chunk_index","cpg_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda")                                
                                  
#remove dupliated rows
AA_simple_results <- AA_simple_results[!duplicated(AA_simple_results), ]

#replace the NA in M2_coef as 0 
AA_simple_results$M2_coef[is.na(AA_simple_results$M2_coef)] <- 0
AA_simple_results$M2_pval[is.na(AA_simple_results$M2_pval)] <- 1

min(AA_simple_results$M1_pval)
min(AA_simple_results$M3_pval)

#extract the AA_simple_results$CpG that need to be rerun: AA_simple_results$M1_coef !=0 OR AA_simple_results$M3_coef !=0
id_rerun <- AA_simple_results[which(AA_simple_results$M1_coef !=0 | AA_simple_results$M3_coef !=0),]



id_rerun <- id_rerun %>%
  mutate(
    Chunk_index = Chunk_index %>%
      str_remove("^AA_") %>%             # drop the prefix
      str_remove("_simple\\.txt$") %>%    # drop the suffix
      as.numeric()                        # convert to numeric
  )



#################################################################################################
rerun_simple_results <- fread(paste0(PATH_wk,"/scr/8_rerun/results/rerun_results_simple.txt"), header=FALSE)

rerun_simple_results[, chunk_id := sub(".*AA_(.*?)_simple\\.txt.*", "\\1", V1)]

#remove the first column
rerun_simple_results <- rerun_simple_results[, -1, with = FALSE]

rerun_simple_results <- rerun_simple_results[!duplicated(rerun_simple_results), ]
colnames(rerun_simple_results) <- c("cpg_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda", "Chunk_index")  


rerun_simple_results$M1_pval[rerun_simple_results$CpG=="1:214170993-214170994"]
id_rerun$M1_pval[id_rerun$CpG=="1:214170993-214170994"]


sum(rerun_simple_results$CpG %in% id_rerun$CpG)
sum(id_rerun$CpG %in% rerun_simple_results$CpG)

#id_rerun[id_rerun$CpG=="1:880109-880110",]


sel_cpgs_ori <- id_rerun$CpG[id_rerun$M1_pval<1e-5]

#the new sig cpgs are these pvalue < 1e-5 and also it has a non-zero coef in the first run (i.e. selected for re-run)
sel_cpgs_new <- intersect( rerun_simple_results$CpG[rerun_simple_results$M1_pval<1e-5], id_rerun$CpG[id_rerun$M1_coef!=0] )

#select these rows cpg from sel_cpgs_new, and Chunk_index from 0934 0693 0690 1180 1167 1116
#and M1_pval < 1e-5
sel_cpgs <- rerun_simple_results %>%
  filter(CpG %in% sel_cpgs_new & Chunk_index %in% c("0934", "0693", "0690", "1180", "1167", "1116") )




head(rerun_simple_results[rerun_simple_results$CpG%in% sel_cpgs_new,])
#these cpgs are significant in both the rerun and the original run 


sel_cpgs <- intersect(sel_cpgs_ori, sel_cpgs_new)

View(rerun_simple_results[rerun_simple_results$CpG%in% sel_cpgs,])


#merge the rerun_simple_results with id_rerun by CpG 
#and Chunk_index
id_rerun_merge <- id_rerun %>%
  left_join(rerun_simple_results, by = c("CpG"))



id_rerun_merge$M1_pval.y[id_rerun_merge$CpG=="1:3094655-3094656"]


#the subset where model 1 need to be rerun
id_rerun_merge_M1 <- id_rerun_merge %>%
  filter(M1_coef.x != 0 ) %>%
  select(Chunk_index.x, cpg_index.x, CpG, M1_coef.x, M1_pval.x, M1_lambda.x, M1_coef.y, M1_pval.y, M1_lambda.y ) %>%
  mutate(model = "M1")

id_rerun_merge_M1[id_rerun_merge_M1$CpG=="1:3094655-3094656",]


id_rerun_merge_M3 <- id_rerun_merge %>%
  filter(M3_coef.x != 0 ) %>%
  select(Chunk_index.x, cpg_index.x, CpG, M3_coef.x, M3_pval.x, M3_lambda.x, M3_coef.y, M3_pval.y, M3_lambda.y ) %>%
  mutate(model = "M3")


plot(id_rerun_merge_M1$M1_pval.x, id_rerun_merge_M1$M1_pval.y)


summary(id_rerun_merge_M1$M1_pval.x < 1e-5)
summary(id_rerun_merge_M1$M1_pval.y < 1e-5)


summary(id_rerun_merge_M3$M3_pval.x < 1e-4)
summary(id_rerun_merge_M3$M3_pval.y < 1e-4)



#for those rows where d_rerun_merge_M1$M1_coef.y is NA, replace these rows original M1_coef.x
id_rerun_merge_M1$M1_coef.y[is.na(id_rerun_merge_M1$M1_coef.y)] <- id_rerun_merge_M1$M1_coef.x[is.na(id_rerun_merge_M1$M1_coef.y)]
id_rerun_merge_M1$M1_pval.y[is.na(id_rerun_merge_M1$M1_pval.y)] <- id_rerun_merge_M1$M1_pval.x[is.na(id_rerun_merge_M1$M1_pval.y)]


#for those rows id_rerun_merge_M1$M1_pval.y = 0, replace the M1_pval.y with the original M1_pval.x
id_rerun_merge_M1$M1_pval.y[id_rerun_merge_M1$M1_pval.y==0] <- id_rerun_merge_M1$M1_pval.x[id_rerun_merge_M1$M1_pval.y==0]
id_rerun_merge_M1$M1_coef.y[id_rerun_merge_M1$M1_pval.y==0] <- id_rerun_merge_M1$M1_coef.x[id_rerun_merge_M1$M1_pval.y==0]


sum(id_rerun_merge_M1$M1_pval.x<1e-5)
sum(id_rerun_merge_M1$M1_pval.y<1e-5)

plot(id_rerun_merge_M1$M1_pval.x, id_rerun_merge_M1$M1_pval.y)

M1_manhattan <- data.frame(CpG = id_rerun_merge_M1$CpG,
                          pval = id_rerun_merge_M1$M1_pval.y,
                          coef = id_rerun_merge_M1$M1_coef.y,
                          pval_old = id_rerun_merge_M1$M1_pval.x)


#id_rerun_merge_M1$M1_pval.y[id_rerun_merge_M1$CpG=="1:3094655-3094656"]


f<-intersect(which(M1_manhattan$pval<5e-8), which(M1_manhattan$pval>0))

M1_manhattan$CpG[f]



#for those rows where d_rerun_merge_M3$M3_coef.y is NA, replace these rows original M3_coef.x
sum(is.na(id_rerun_merge_M3$M3_coef.y))
sum(is.na(id_rerun_merge_M3$M3_pval.y))

#id_rerun_merge_M3$M3_coef.y[is.na(id_rerun_merge_M3$M3_coef.y)] <- id_rerun_merge_M3$M3_coef.x[is.na(id_rerun_merge_M3$M3_coef.y)]
#id_rerun_merge_M3$M3_pval.y[is.na(id_rerun_merge_M3$M3_pval.y)] <- id_rerun_merge_M3$M3_pval.x[is.na(id_rerun_merge_M3$M3_pval.y)]


#for those rows id_rerun_merge_M3$M3_pval.y = 0, replace the M3_pval.y with the original M3_pval.x
sum(id_rerun_merge_M3$M3_pval.y==0)
#id_rerun_merge_M3$M3_pval.y[id_rerun_merge_M3$M3_pval.y==0] <- id_rerun_merge_M3$M3_pval.x[id_rerun_merge_M3$M3_pval.y==0]
#id_rerun_merge_M3$M3_coef.y[id_rerun_merge_M3$M3_pval.y==0] <- id_rerun_merge_M3$M3_coef.x[id_rerun_merge_M3$M3_pval.y==0]




M3_manhattan <- data.frame(CpG = id_rerun_merge_M3$CpG,
                          pval = id_rerun_merge_M3$M3_pval.y,
                          coef = id_rerun_merge_M3$M3_coef.y,
                          plva_old = id_rerun_merge_M3$M3_pval.x)

sum(M3_manhattan$pval==0)
sum(M3_manhattan$pval<1e-5)

save(M1_manhattan,M3_manhattan, file = paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"))

M2_cpgs <- AA_simple_results[which(AA_simple_results$M2_pval<1e-5),]
save(M1_manhattan,M3_manhattan,M2_cpgs,file = paste0(PATH_wk,"/scr/8_rerun/results/4_eval_M123.RData"))

#no need for M2 manhattan plot, since M2 is not rerurn. we just need to reload their results. 

load(paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"),verbose=TRUE)

M1_manhattan$pval[M1_manhattan$CpG=="1:3094655-3094656"]
M1_manhattan$pval[M1_manhattan$CpG=="6:31745741-NA"]

#add all those cpgs in M1_manhattan that are not in M3_manhattan to M3_manhattan, and put their coef as 0 and pval as 1
M3_manhattan <- M3_manhattan %>%
  bind_rows(M1_manhattan[!(M1_manhattan$CpG %in% M3_manhattan$CpG),] %>%
              mutate(coef = 0, pval = 1))



source(paste0(PATH_wk,"scr/7_eval/func_plot_manhattan.R"))


man_plot1 <- man_plot(M1_manhattan, name ="M1_manhattan")
man_plot3 <- man_plot(M3_manhattan, name ="M3_manhattan")
#man_plot2 <- man_plot(M2_manhattan, name ="M2_manhattan")

sum(M1_manhattan$pval < 1e-5)
sum(M1_manhattan$pval < 1e-4)

###########################################################################

sum(M1_manhattan$pval < 1e-5)
sum(M2_manhattan$pval < 1e-5)
sum(M3_manhattan$pval < 1e-5)


###########################################################################
# Assume models_res is a matrix with dimensions 6 x nCpGs,
# where:
# - Row 1: beta estimates for Model 1 (significant if != 0)
# - Row 3: p-values for Model 2 (significant if < 0.05)
# - Row 5: beta estimates for Model 3 (significant if != 0)

# Extract significant CpG indices for each model:
sig_model1 <- M1_manhattan$CpG[which(M1_manhattan$pval<1e-5)]
sig_model2 <- AA_simple_results$CpG[which(AA_simple_results$M2_pval<1e-5)]
sig_model3 <- M3_manhattan$CpG[which(M3_manhattan$pval<1e-5)]

# Create a list of the significant sets
sig_list <- list(
  Model1 = sig_model1,
  Model2 = sig_model2,
  Model3 = sig_model3)

# Plot the Venn diagram.
venn.plot <- venn.diagram(
  x = sig_list,
  filename = NULL,  # if NULL, the diagram is returned as a grid object
  fill = c("red", "blue","green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of Significant CpGs from three models (pvalue < 1e-5)"
)


# Open a JPEG device to save the plot
jpeg(paste0(PATH_wk, "/scr/8_rerun/results/venn_123.jpeg"), width = 800, height = 800)
# Draw the Venn diagram on the device
grid.draw(venn.plot)
# Close the device to finalize the file
dev.off()
###########################################################################


###########################################################################
# veen diagram of the top 1000 significant CpGs from three models

# where:
# - Row 1: beta estimates for Model 1 (significant if != 0)
# - Row 3: p-values for Model 2 (significant if < 0.05)
# - Row 5: beta estimates for Model 3 (significant if != 0)

# Extract significant CpG indices for each model:
sig_model1 <- M1_manhattan$CpG[order(M1_manhattan$pval, decreasing = FALSE)[1:1000]]
sig_model2 <- AA_simple_results$CpG[order(AA_simple_results$M2_pval, decreasing = FALSE)[1:1000]]
sig_model3 <- M3_manhattan$CpG[order(M3_manhattan$pval, decreasing = FALSE)[1:1000]]


# Create a list of the significant sets
sig_list <- list(
  Model1 = sig_model1,
  Model2 = sig_model2,
  Model3 = sig_model3
)

# Plot the Venn diagram.
venn.plot <- venn.diagram(
  x = sig_list,
  filename = NULL,  # if NULL, the diagram is returned as a grid object
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of Top 1000 significant CpGs from three models"
)


# Open a JPEG device to save the plot
jpeg(paste0(PATH_wk, "/scr/8_rerun/results/venn_123_top1k.jpeg"), width = 800, height = 800)

# Draw the Venn diagram on the device
grid.draw(venn.plot)
# Close the device to finalize the file
dev.off()

