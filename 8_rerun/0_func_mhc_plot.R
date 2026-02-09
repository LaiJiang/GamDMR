#input: mhc_cpgs with start and pval two columns 
#input: plot_name


#output: a jpeg file with the manhattan plot



# install.packages("ggplot2")   # if you haven’t already
library(ggplot2)

# Suppose your mhc_cpgs data frame looks like this:
#    Chunk_index               CpG   chr    start    end      pval    M2_coef  M2_F_sdv
# 1: AA_1420_simple.txt 6:28477957-28477958   6 28477957 28477958 0.2258151 -0.019668 0.0148153
# … etc …

# 1) Prepare the data
mhc_cpgs$pos <- mhc_cpgs$start             # genomic coordinate for plotting
mhc_cpgs$logp <- -log10(mhc_cpgs$pval)      # –log10(p-value)



# 2) Basic Manhattan‐style plot with ggplot2
ggplot(mhc_cpgs, aes(x = pos, y = logp)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +  # genome-wide sig.
  scale_x_continuous(
    name = "Genomic position from MHC region",
    labels = scales::comma_format(), 
    breaks = scales::pretty_breaks(10)
  ) +
  scale_y_continuous(name = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste0(plot_name, " : Manhattan plot for CpGs from MHC region chr6"))

#save this plot as a jpeg
ggsave(filename = paste0(PATH_wk,"/scr/8_rerun/results/8_zoom_MHC_",plot_name,".jpeg"),
       width = 10, height = 6, dpi = 300)


