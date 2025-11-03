suppressMessages(suppressWarnings({
  library(Seurat)
  library(tidyverse)
  library(ggpubr)
  library(fgsea)
  library(RColorBrewer)
  library(rstatix)
  library(cowplot)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "7_CNA"))

#### Read the master seurat object ####
seu = readRDS(paste0(base_dir, "5_cNMF_analysis/annotated_seu.rds"))
seu = subset(seu, Cohort %in% c("Aged", "Young"))

cna_results = read.table("cna_results.csv", header = T, sep = ",", row.names = 1)
all(rownames(cna_results) == seu$Cell)
seu = AddMetaData(seu, cna_results)

seu$cohort_coef[seu$cohort_coef_fdr >= 0.1] = NA
p = scCustomize::FeaturePlot_scCustom(seu, pt.size = 0.001,
                                      "cohort_coef", 
                                      na_cutoff = NA, 
                                      order = F) &
  coord_fixed() & 
  theme_pubr(base_size = 5) &
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        legend.justification = "left",
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=5)) & 
  scale_color_gradient2(
    low = "steelblue2",
    high = "firebrick2",
    mid = "white",
    na.value = "gray90"
    
    # breaks = seq(-0.4, 1.2, by = 0.2)  # Uniform break intervals
  )

p = ggrastr::rasterize(p, dpi=1000, dev="ragg_png")
save_plot("Fig.1f.pdf", p, base_height = 90, base_width = 90, units="mm")
