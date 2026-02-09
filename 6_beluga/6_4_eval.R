PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


recrd_cpgs <- read.csv(paste0(PATH_wk,"scr/6_beluga/record.csv"), header=TRUE)


sum(recrd_cpgs$last_index != 3000)

file_uncomplete <- recrd_cpgs$filename[recrd_cpgs$last_index != 3000]


file_uncomplete <- gsub("AA_","",file_uncomplete)
file_uncomplete <- gsub("_simple.txt","",file_uncomplete)
file_uncomplete <- as.numeric(file_uncomplete)


#write file_uncomplete to a txt file, one number per line, unquoted
write.table(file_uncomplete, file = paste0(PATH_wk,"scr/6_beluga/file_uncomplete.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)