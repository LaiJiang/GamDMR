#check the region file from 4_merge.R on Beluga. and split into two sets, one with data_chunk_id from a single chunk, the second part with data_chunk_id from 2 chunks.
#we will write script to process these two types seperately.

library(data.table)
library(dplyr)

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
region_file <- fread("C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/dat/merged_regions.csv")


sum(region_file$n_cpgs > 500)

plot(hist(region_file$n_cpgs), main = "Density of n_cpgs in merged regions", xlab = "n_cpgs", ylab = "Density")

table(region_file$n_cpgs)


list_chunks <-unique(region_file$data_chunk_id)

region_file[region_file$data_chunk_id=="0001,0002",]

sum(region_file$data_chunk_id=="0001")

dim(region_file)
region_file[c(45:48),]

#stopped here. how to run mgcv on each of the regions? !!!!!!!!!!!!!
#to assign 1000 jobs, each with 100 regions, take 3 hours.
#when not finnisehd, reverse the region id rank to get it done in total 6 hours.

length(unique(region_file$data_chunk_id))
#split region_file into two parts, one with data_chunk_id from a single chunk. the second part
#with data_chunk_id from 2 chunks.


#the first row of each chunk 
chunk_info <- read.table(file=paste0(PATH_wk,"results/10_sanity/first_row_chunk.txt"))
chunk_info$data_chunk_id <- sprintf("%04d",1:nrow(chunk_info) )
#first extract rows in region_file that are from 2 chunks



region_file_2_chunks <- region_file[grepl(",", region_file$data_chunk_id),]
#then extract the  rows in region_file that are from a single chunk
region_file_1_chunk <- region_file[!grepl(",", region_file$data_chunk_id),]

#for rows in region_file_2_chunks, we check if the data_chunk_id is really from two chnks, or because we splitted the region after record data_chunkId 
region_file_2_chunks$chunk1 <- sapply(strsplit(region_file_2_chunks$data_chunk_id, ","), `[`, 1)
region_file_2_chunks$chunk2 <- sapply(strsplit(region_file_2_chunks$data_chunk_id, ","), `[`, 2)


# Create a named vector for fast lookup: data_chunk_id → V2 (start position of chunk)
chunk_start_map <- setNames(chunk_info$V2, chunk_info$data_chunk_id)

# Add valid_extend column
region_file_2_chunks <- region_file_2_chunks %>%
  mutate(
    chunk2_start = chunk_start_map[chunk2],  # already defined earlier
    correct_data_chunk_id = ifelse(
      region_start < chunk2_start & region_end >= chunk2_start, 
      data_chunk_id,
      ifelse(region_end < chunk2_start, chunk1, chunk2)
    )
  )%>% #remove columns data_chunk_id, chunk1,chunk2, chunk2_start,valid_extend
    select(-data_chunk_id, -chunk1, -chunk2, -chunk2_start) %>% rename(data_chunk_id = correct_data_chunk_id)

#now merge the two data frames
region_file_correct <- rbind(region_file_1_chunk, region_file_2_chunks)



region_file_2_chunks <- region_file_correct[grepl(",", region_file_correct$data_chunk_id),]
#then extract the  rows in region_file that are from a single chunk
region_file_1_chunk <- region_file_correct[!grepl(",", region_file_correct$data_chunk_id),]



#for rows in region_file_2_chunks, we check if the data_chunk_id is really from two chnks, or because we splitted the region after record data_chunkId 
region_file_2_chunks$chunk1 <- sapply(strsplit(region_file_2_chunks$data_chunk_id, ","), `[`, 1)
region_file_2_chunks$chunk2 <- sapply(strsplit(region_file_2_chunks$data_chunk_id, ","), `[`, 2)


# Create a named vector for fast lookup: data_chunk_id → V2 (start position of chunk)
chunk_start_map <- setNames(chunk_info$V2, chunk_info$data_chunk_id)

# Add valid_extend column
region_file_2_chunks <- region_file_2_chunks %>%
  mutate(
    chunk2_start = chunk_start_map[chunk2],  # already defined earlier
    correct_data_chunk_id = ifelse(
      region_start < chunk2_start & region_end >= chunk2_start, 
      data_chunk_id,
      ifelse(region_end < chunk2_start, chunk1, chunk2)
    )
  )

#save two files
fwrite(region_file_1_chunk, "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/dat/region_file_1_chunk.csv"
, row.names = FALSE
, quote = FALSE)

fwrite(region_file_2_chunks, "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/dat/region_file_2_chunks.csv"
, row.names = FALSE
, quote = FALSE)