suppressMessages(suppressWarnings({
  library(RColorBrewer)
  library(fgsea)
  library(tidyverse)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

if (!dir.exists("meta_clusters_starcat")) dir.create("meta_clusters_starcat")

#### Read input files ####

# Read cNMF programs and filter out singletons
cNMF_results = read.table("cnmf_output/Z_scores.total_cNMF_output.csv", row.names = 1, header = T, sep=",")
cNMF_results$Sample = NULL
genes = colnames(cNMF_results)
singletons = read.table("singletons.txt",header = F)$V1
cNMF_results = cNMF_results[! rownames(cNMF_results) %in% singletons, ]

# Read TPM files to produce averaged TPM meta-programs for *CAT
cGEP_TPMs = read.table("cnmf_output/TPMs.total_cNMF_output.csv", row.names = 1, header = T, sep=",")
cGEP_TPMs = cGEP_TPMs[! rownames(cGEP_TPMs) %in% singletons, ]

# Read cluster assignment
clusters = read.table("meta_clusters_starcat/clusters.csv", sep=",", header = T)

#### Subset representative clusters and compute averaged Z-scores and TPM vectors for them ####

# Compute Pearson correlation for visualization
full_similarity_matrix = cor(scale(t(data.matrix(cNMF_results))))

# Select clusters that have at least four programs coming from at least two datasets
selected_clusters = clusters %>%
  filter(Cluster != "Unclustered") %>%
  rowwise() %>% 
  mutate(Dataset = strsplit(Member, "_")[[1]][2]) %>%
  group_by(Cluster) %>%
  filter(n_distinct(Dataset)>1 & n()>3) %>%
  dplyr::select(-Dataset) %>%
  mutate(Cluster = paste0("metaGEP:", Cluster))

# Make averaged Z-score vectors for GSEA analysis
averaged_metaGEP = cNMF_results %>%
  rownames_to_column("Member") %>%
  left_join(selected_clusters) %>%
  filter(! is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(across(2:(ncol(cNMF_results)+1), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>% 
  column_to_rownames("Cluster")

write.table(averaged_metaGEP, "meta_clusters_starcat/metaGEP.Z_scores.tsv", quote = F, append = F, row.names = T, col.names = T, sep = "\t")

# Make averaged TPM vectors for *CAT
averaged_metaGEP_TPMs = cGEP_TPMs %>%
  dplyr::select(-Sample) %>%
  rownames_to_column("Member") %>%
  left_join(selected_clusters) %>%
  filter(! is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(across(2:(ncol(cNMF_results)+1), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("Cluster")

write.table(averaged_metaGEP_TPMs, "meta_clusters_starcat/metaGEP.Z_score_based_TPMs.tsv", quote = F, append = F, row.names = T, col.names = T, sep = "\t")

# Subset correlation matrix for selected clusters 
cluster_members = selected_clusters$Member
similarity_matrix = full_similarity_matrix[cluster_members, cluster_members]

#### Plot Fig.2b ####

# Define cluster colors
cluster_names = selected_clusters$Cluster
unique_clusters = unique(cluster_names)
color_palette = c(hcl.colors(length(unique_clusters), "Roma"))
color_mapping = setNames(color_palette, unique_clusters)
cluster_colors = color_mapping[cluster_names] 

# Define dataset colors
datasets = stringr::str_extract(rownames(similarity_matrix), "_([A-Za-z0-9]+)_", group = 1)
unique_datasets = unique(datasets)
color_palette = brewer.pal(8, "Paired")[c(1,2,8,6)]
color_palette[3] = "#ffcb00"
color_mapping = setNames(color_palette, unique_datasets)
dataset_colors = color_mapping[datasets] 

# Minor graphic tweaks
ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(0.5, "mm")
ht_opt$HEATMAP_LEGEND_PADDING  = unit(0.5, "mm")

# Define annotation block
anno = columnAnnotation(
  Cluster = cluster_names, 
  Dataset = datasets,
  col = list(
    Cluster = cluster_colors,
    Dataset = dataset_colors),
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
col_fun = colorRamp2(c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), colors = c(hcl.colors(6, "Teal")[2:5], "white", rev(hcl.colors(6, "Oranges")[2:5])))

# Plot heatmap
hm = Heatmap(similarity_matrix,
             width = unit(1, "snpc"), height = unit(1, "snpc"),
             col = col_fun,
             row_order = selected_clusters$Member,
             column_order = selected_clusters$Member,
             show_column_names = F, 
             show_row_names = F,
             show_heatmap_legend = T,
             bottom_annotation = anno,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             show_parent_dend_line = T,
             heatmap_legend_param = list(
               direction = "vertical",
               labels_gp = gpar(fill = dataset_colors, fontsize = 5),
               title_gp = gpar(fontsize = 6, fontface = "bold"),
               title = "Correlation",
               legend_height = unit(7, "mm"),
               grid_width = unit(1, "mm")),
             row_gap = unit(0, "mm"), 
             column_gap = unit(0, "mm"))

hm

# Save as pdf
pdf("meta_clusters_starcat/Fig.2b.pdf", width = 3.3, height = 2.5)
draw(hm, annotation_legend_side = "right", merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

#### Plot all programs with their GSEA annotation for supplements ####

# Reorder correlation matrix entries
full_similarity_matrix = full_similarity_matrix[clusters$Member, clusters$Member]

# Read and reformat GSEA results
gsea_res = data.table::fread(paste0("gsea_results/total_cNMF_fgseaRes.txt"), sep="\t", sep2=c("", " ", ""))

gsea_wide = gsea_res %>%
  dplyr::select(pathway, padj, NES, Program) %>%
  filter(pathway %in% c("Aging_signature", "MEP_signature", "GMP_signature", "CLP_signature")) %>%
  mutate(NES = - sign(NES) * log(padj)) %>%
  dplyr::select(-padj) %>%
  pivot_wider(values_from = NES, names_from = pathway)

gsea_anno = gsea_wide[match(rownames(full_similarity_matrix), gsea_wide$Program), ]

# GSEA colorbar
col_fun_gsea = colorRamp2(c(min(min(gsea_res$NES, na.rm = T), log(0.05)-1), 0, max(max(gsea_res$NES, na.rm = T), -log(0.05)+1)), c("steelblue3", "white","firebrick3"))

# Define GSEA annotation block
row_anno = rowAnnotation(
  GMP_sig = gsea_anno$GMP_signature,
  MEP_sig = gsea_anno$MEP_signature,
  MLP_sig = gsea_anno$CLP_signature,
  Aging_sig = gsea_anno$Aging_signature,
  col = list(
    GMP_sig = col_fun_gsea,
    MEP_sig = col_fun_gsea,
    MLP_sig = col_fun_gsea,
    Aging_sig = col_fun_gsea),
  show_legend = T,
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    GMP_sig = list(direction = "vertical"),
    MEP_sig = list(direction = "vertical"),
    MLP_sig = list(direction = "vertical"),
    Aging_sig = list(direction = "vertical"),
    grid_height = unit(1, "mm"),
    grid_width = unit(1, "mm"),
    row_gap = unit(0.1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5)
  ),
  simple_anno_size = unit(1.5, "mm")
)

# Define cluster colors
cluster_names = clusters$Cluster
unique_clusters = unique(cluster_names)
color_palette = c(hcl.colors(length(unique_clusters) - 1, "Roma"), "grey")
color_mapping = setNames(color_palette, unique_clusters)
cluster_colors = color_mapping[cluster_names] 

# Define dataset colors
datasets = stringr::str_extract(rownames(full_similarity_matrix), "_([A-Za-z0-9]+)_", group = 1)
unique_datasets = unique(datasets)
color_palette = brewer.pal(8, "Paired")[c(1,2,8,6)]
color_palette[3] = "#ffcb00"
color_mapping = setNames(color_palette, unique_datasets)
dataset_colors = color_mapping[datasets] 

# Define top annotation block
column_anno = columnAnnotation(
  Cluster = cluster_names, 
  Dataset = datasets,
  col = list(
    Cluster = cluster_colors,
    Dataset = dataset_colors),
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

# Plot heatmap
hm = Heatmap(full_similarity_matrix,
             width = unit(1, "snpc"), height = unit(1, "snpc"),
             col = col_fun,
             row_order = clusters$Member,
             column_order = clusters$Member,
             show_column_names = F, 
             show_row_names = F,
             show_heatmap_legend = T,
             top_annotation = column_anno,
             right_annotation = row_anno,
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             show_parent_dend_line = T,
             heatmap_legend_param = list(
               direction = "vertical",
               labels_gp = gpar(fill = dataset_colors, fontsize = 5),
               title_gp = gpar(fontsize = 6, fontface = "bold"),
               title = "Correlation",
               legend_height = unit(7, "mm"),
               grid_width = unit(1, "mm")),
             row_gap = unit(0, "mm"), 
             column_gap = unit(0, "mm"))

hm

# Save as pdf
pdf("meta_clusters_starcat/Suppl.Fig.5a.pdf", width = 160/25.4, height = 120/25.4)
draw(hm, annotation_legend_side = "right", merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()
