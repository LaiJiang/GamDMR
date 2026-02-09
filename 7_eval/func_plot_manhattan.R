library(qqman)

man_plot <- function(man_dat, name){# open a JPEG device; adjust width/height/resolution as desired

#pval smaller than 1e-10 forced to 1e-10
man_dat$pval[man_dat$pval < 1e-20] <- 1e-20

# extract everything before ':' as chr
man_dat$chr <- sub(":.*", "", man_dat$CpG)
# extract the number between ':' and '-' as pos, then convert to integer
man_dat$pos <- as.integer(sub(".*:(\\d+)-.*", "\\1", man_dat$CpG))
# make sure chr is numeric
man_dat$chr <- as.numeric(man_dat$chr)



jpeg(
  filename = paste0("C:/Per/LaiJiang/Project/UQAC/meth/scr/7_eval/results/",name, ".jpg"), 
  width    = 1200, 
  height   = 800, 
  units    = "px", 
  res      = 150
)

# draw the Manhattan plot
manhattan(
  man_dat,
  chr            = "chr",
  bp             = "pos",
  p              = "pval",
  snp            = "CpG",
  genomewideline = -log10(5e-8),
  suggestiveline =  -log10(1e-5),
  col            = c("blue", "black")
)


# close the device to write the file
dev.off()

}