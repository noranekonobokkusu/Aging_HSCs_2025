suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

meta = readRDS(paste0(base_dir, "1_dataset_preparation/seu_hsc.total.rds"))@meta.data

clusters = read.table(paste0("meta_clusters_starcat/clusters.csv"), sep=",", header = T)
cluster_members_df = clusters %>%
  rowwise() %>% 
  mutate(Dataset = strsplit(Member, "_")[[1]][2]) %>%
  group_by(Cluster) %>%
  filter(n_distinct(Dataset)>1 & n()>3) %>%
  dplyr::select(-Dataset) %>%
  filter(Cluster != "Unclustered") %>%
  dplyr::select(Member, Cluster) %>%
  arrange(Cluster)
cluster_members = cluster_members_df$Member
singletons = read.table(paste0("singletons.txt"),header = F)$V1

cNMF_results_Z_scores = read.table(paste0("cnmf_output/Z_scores.total_cNMF_output.csv"), row.names = 1, header = T, sep=",")
cNMF_results_Z_scores$Sample = NULL
genes = colnames(cNMF_results_Z_scores)
cNMF_results_Z_scores = cNMF_results_Z_scores[! rownames(cNMF_results_Z_scores) %in% singletons, ]
data_type = "Z_scores"
cNMF_results_TPMs = read.table(paste0("cnmf_output/TPMs.total_cNMF_output.csv"), row.names = 1, header = T, sep=",")
cNMF_results_TPMs$Sample = NULL
genes = colnames(cNMF_results_TPMs)
cNMF_results_TPMs = cNMF_results_TPMs[! rownames(cNMF_results_TPMs) %in% singletons, ]

averaged_metaGEP_Z_scores = read.table(paste0("meta_clusters_starcat/metaGEP.Z_scores.tsv"), row.names = 1, header = T, sep = "\t")
averaged_metaGEP_TPMs = read.table(paste0("meta_clusters_starcat/metaGEP.Z_score_based_TPMs.tsv"), row.names = 1, header = T, sep = "\t")

top_genes = averaged_metaGEP_Z_scores %>%
  rownames_to_column("metaGEP") %>%
  pivot_longer(-metaGEP) %>%
  group_by(metaGEP) %>%
  arrange(desc(value), .by_group = T) %>%
  slice_head(n=10) %>%
  pull(name) %>%
  unique()

heatmap_df = cNMF_results_Z_scores[cluster_members, top_genes]

cluster_names = cluster_members_df$Cluster
unique_clusters = sort(unique(cluster_names[cluster_names!="Unclustered"]))
color_palette = c(hcl.colors(length(unique_clusters), "Roma"), "gray")
unique_clusters = c(unique_clusters, "Unclustered")
color_mapping = setNames(color_palette, unique_clusters)
cluster_colors = color_mapping[cluster_names] 

# Define colors of row names
numbers = stringr::str_extract(cluster_members, "_([A-Za-z0-9]+)_", group = 1)
unique_numbers = unique(numbers)
color_palette = brewer.pal(8, "Paired")[c(1,2,8,6)]
color_palette[3] = "#ffcb00"
color_mapping = setNames(color_palette, unique_numbers)
row_colors = color_mapping[numbers] 

ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(0.5, "mm")
ht_opt$HEATMAP_LEGEND_PADDING  = unit(0.5, "mm")

row_anno = rowAnnotation(
  Cluster = cluster_names, 
  Dataset = numbers,
  col = list(
    Cluster = cluster_colors,
    Dataset = row_colors),
  show_legend = T,
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    Cluster = list(direction="vertical"),
    Dataset = list(direction="vertical"),
    grid_height = unit(1, "mm"),
    grid_width = unit(1, "mm"),
    row_gap = unit(0.1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5)
  ),
  simple_anno_size = unit(1.5, "mm")
)

# Define colorscale for the heatmap
col_fun = colorRamp2(c(min(heatmap_df), 0, max(heatmap_df)), c("#5d928d", "white", "#b85023"))
hm = Heatmap(heatmap_df,
             cluster_rows = F, 
             cluster_columns = F,
             row_split = cluster_members_df$Cluster,
             row_title = c("", "", "", ""),
             right_annotation = row_anno,
             col = col_fun,
             show_row_names = F,
             show_heatmap_legend = T,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             heatmap_legend_param = list(direction = "vertical",
                                         labels_gp = gpar(fill = row_colors, fontsize = 5),
                                         title_gp = gpar(fontsize = 6, fontface = "bold"),
                                         title = "Z-score",
                                         legend_height = unit(7, "mm"),
                                         grid_width = unit(1, "mm")),
             row_gap = unit(1, "mm"), column_gap = unit(1, "mm")
)

hm
pdf("meta_clusters_starcat/Suppl.Fig.5b.pdf", width = 120/25.4, height = 120/25.4)
draw(hm, annotation_legend_side = "right", merge_legend = TRUE, heatmap_legend_side = "right", padding = unit(c(0, 0, 2, 0), "mm"))
dev.off()

