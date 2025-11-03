suppressMessages(suppressWarnings({
  library(RColorBrewer)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "6_SCENIC"))

# Read Z-scores and subset AUC matrices to shared regulons
mtx_files = list.files("output", pattern = "*.auc_mtx.thr_0.01.csv", recursive = TRUE, full.names = T)

auc_mtx_list = lapply(mtx_files, FUN = function(x) { read.csv(x, header = T, row.names = 1) } )
names(auc_mtx_list) = gsub(".*/([a-zA-Z0-9]*)\\..*", "\\1", mtx_files)

regulon_names = lapply(auc_mtx_list, function(x) {colnames(x)})

common_regulons = purrr::reduce(regulon_names, intersect)
common_regulon_auc_mtx_list = lapply(auc_mtx_list, FUN = function(x) { scale(x[,common_regulons]) } )

ap1_complex = c("FOS", "FOSB", "JUN", "JUNB", "JUND", "FOSL1", "FOSL2", "ATF2", "ATF3", "ATF4", "ATF6", "ATF7", "BATF", "BATF2", "BATF3", "MAFA", "MAFB", "MAF", "MAFG", "MAFF", "MAFK")

names(common_regulon_auc_mtx_list) = NULL
common_regulon_auc_mtx = do.call(rbind, common_regulon_auc_mtx_list)

# Make one single long data frame with all AUC Z-scores
auc_long_df = common_regulon_auc_mtx %>%
  as.data.frame() %>%
  rownames_to_column("Cell") %>%
  pivot_longer(-Cell)

# Read metadata and add it to the AUC data frame
meta = readRDS(paste0(base_dir, "5_cNMF_analysis/annotated_seu.rds"))@meta.data

auc_long_df = auc_long_df %>%
  left_join(meta, by="Cell") %>%
  dplyr::rename("Regulon"="name",
                "Regulon_AUC_Z_score"="value")

# Identify top regulons in Young and Aged cohorts
regulons_median_ranks_and_cors = auc_long_df %>%
  mutate(Regulon = gsub("\\..*", "", Regulon),
         Cohort_binary = ifelse(Cohort == "Young", 0, 1)) %>% 
  group_by(Dataset, Regulon) %>%
  summarize(cor_value = cor(Regulon_AUC_Z_score, Cohort_binary)) %>%
  group_by(Dataset) %>%
  arrange(desc(cor_value), .by_group = T) %>%
  mutate(rank = row_number()) 
  
for (dataset in unique(meta$Dataset)) {
specific_regulons = regulons_median_ranks_and_cors %>%
  filter(Dataset==dataset) %>%
  slice_head(n=10) %>%
  pull(Regulon)

df = auc_long_df %>%
  filter(Dataset==dataset) %>%
  mutate(Regulon = gsub("\\..*", "", Regulon)) %>%
  filter(Regulon %in% specific_regulons) %>%
  mutate(Cohort = factor(Cohort, levels=c("Young", "Aged"))) %>%
  arrange(Cohort, Sample) 
  
m = df %>%
  dplyr::select(Cell, Regulon, Regulon_AUC_Z_score) %>%
  pivot_wider(names_from = "Regulon", values_from = "Regulon_AUC_Z_score") %>%
  column_to_rownames("Cell")

cohort_names = df %>%
  filter(Regulon %in% specific_regulons[1]) %>%
  pull(Cohort)

color_palette = c("#c7d1ca", "#617a6a")
color_mapping = setNames(color_palette, unique(cohort_names))
cohort_colors = color_mapping[cohort_names] 

sample_names = df %>%
  filter(Regulon %in% specific_regulons[1]) %>%
  pull(Sample)

color_palette = c(brewer.pal(11, "Spectral"), brewer.pal(9, "Set1"))[1:length(unique(sample_names))]
color_mapping = setNames(color_palette, unique(sample_names))
sample_colors = color_mapping[sample_names] 

# Generate annotation
row_anno = rowAnnotation(
  Cohort = cohort_names,
  Sample = sample_names,
  col = list(
    Cohort = cohort_colors,
    Sample = sample_colors),
  show_legend = T,
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    Cluster = list(nrow = 1),
    grid_height = unit(1, "mm"),
    grid_width = unit(1, "mm"),
    row_gap = unit(0.1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5)),
  simple_anno_size = unit(1.5, "mm")
)
m = as.matrix(m)
col_fuc = colorRamp2(c(quantile(m, 0.001), 0, quantile(m, 0.999)), hcl_palette = "YlGnBu")


# Generate heatmap
hm = Heatmap(m,
             right_annotation = row_anno, 
             col = col_fuc,
             cluster_rows = F,
             cluster_columns = T,
             show_column_dend = F,
             show_row_names = F,
             name = "Z-score",
             row_title_gp = gpar(fontsize = 5),
             row_names_side = "right",
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             heatmap_legend_param = list(direction = "vertical",
                                         labels_gp = gpar(fontsize = 5),
                                         title_gp = gpar(fontsize = 6, fontface = "bold"),
                                         title = "Z-score",
                                         legend_height = unit(7, "mm"),
                                         grid_width = unit(1, "mm")),
             border = T,
             border_gp = gpar(lwd = 0.5)
        
)

hm
pdf(paste0("figures/", dataset, ".top_aging_tf_activity.pdf"), width = 3, height = 3)
draw(hm, merge_legend = TRUE, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()
}
