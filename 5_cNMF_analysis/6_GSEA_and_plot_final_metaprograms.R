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

##### Create a list of pathways #####

# Aging signature
de_results = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep=",", header = T)
aging_sig = de_results %>%
  filter(padj<0.05, log2FoldChange>1) %>%
  pull(gene)

# Human dormancy signature from https://doi.org/10.1016/j.stem.2021.07.003
garcia_table = readxl::read_xlsx(paste0(base_dir, "input_data/signatures/mmc2.xlsx"), sheet = 1)
garcia_signature = garcia_table %>%
  filter(genename!="None") %>%
  filter(logFC < -2,FDR<0.00001) %>%
  pull(genename) %>%
  unique()

graham_signature = read.table(paste0(base_dir, "input_data/signatures/graham_quiescence_up.txt"), header = F)$V1

mep_sig = read.table(paste0(base_dir, "input_data/signatures/mep_sig.tsv"))$x
gmp_sig = read.table(paste0(base_dir, "input_data/signatures/gmp_sig.tsv"))$x
clp_sig = read.table(paste0(base_dir, "input_data/signatures/clp_sig.tsv"))$x

ap1_complex = c("FOS", "FOSB", "JUN", "JUNB", "JUND", "FOSL1", "FOSL2", "ATF2", "ATF3", "ATF4", "ATF6", "ATF7", "BATF", "BATF2", "BATF3", "MAFA", "MAFB", "MAF", "MAFG", "MAFF", "MAFK")

hallmark_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/h.all.v2024.1.Hs.symbols.gmt.txt"))
c2_cp_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/c2.cp.v2024.1.Hs.symbols.gmt.txt"))

main_pathways = list(graham_signature, garcia_signature, aging_sig, mep_sig, gmp_sig, clp_sig, ap1_complex)

names(main_pathways) = c("Graham_quiescence", "Garcia_quiescence", "Aging_signature", "MEP_signature", "GMP_signature", "CLP_signature", "AP-1_complex")

pathways = c(main_pathways, hallmark_pathways, c2_cp_pathways)

#### Run GSEA and visualise enrichment scores for metaGEPs ####

# Read metaGEPs
programs = read.table("meta_clusters_starcat/metaGEP.Z_scores.tsv", row.names = 1, header = T, sep="\t")

fgseaRes_total = data.frame()

# Run GSEA for each metaGEP
for (i in 1:nrow(programs)) {
  entry = programs[i,]
  ranks = as.numeric(entry)
  names(ranks) = colnames(entry)
  
  # GSEA
  fgseaRes = fgsea(pathways = pathways, 
                   stats    = ranks,
                   minSize  = 5,
                   eps      = 0.0,
                   maxSize  = 500,
                   nPermSimple = 10000)
  
  fgseaRes$Program = rownames(entry)
  fgseaRes_total = rbind(fgseaRes_total, fgseaRes)
  
}

# Save results
data.table::fwrite(fgseaRes_total, file=paste0("gsea_results/metaGEP_fgseaRes.txt"), sep="\t", sep2=c("", " ", ""), append = F)

# Select relevant pathways to visualize
relevant_pathways = c("Aging_signature", "Graham_quiescence", "Garcia_quiescence", "MEP_signature", "GMP_signature", "CLP_signature", "AP-1_complex", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_CELL_CYCLE", "HALLMARK_E2F_TARGETS", "HALLMARK_APOPTOSIS", "WP_IL18_SIGNALING", "WP_NUCLEAR_RECEPTORS_METAPATHWAY", "PID_AP1_PATHWAY", "HALLMARK_HYPOXIA", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INFLAMMATORY_RESPONSE")

# Create an intermediate df with both NESs and significance signs, for all combinations of programs and gene sets
temp = fgseaRes_total %>%
  filter(pathway %in% relevant_pathways) %>%
  arrange(pathway = factor(pathway, levels = relevant_pathways)) %>%
  dplyr::select(ID=pathway, padj, NES, Program) %>%
  mutate(sign = case_when(padj < 0.001 ~ "***",
                          padj < 0.01 ~ "**",
                          padj < 0.05 ~ "*",
                          .default = ""))

# Extract NES and significance levels into separate objects
NES_df = temp %>%
  pivot_wider(names_from = Program, values_from = NES, id_cols = ID)
p_val_df = temp %>%
  pivot_wider(names_from = Program, values_from = sign, id_cols = ID)

# Move gene names to rownames
ids = NES_df[,1]$ID
NES_df = as.matrix(NES_df[,-1])
rownames(NES_df) = ids
NES_df = t(NES_df)

ids = p_val_df[,1]$ID
p_val_df = as.matrix(p_val_df[,-1])
rownames(p_val_df) = ids
p_val_df = t(p_val_df)

# Make sure both cell types and gene sets in the same order in two dfs
all(rownames(NES_df) == rownames(p_val_df))
all(colnames(NES_df) == colnames(p_val_df))

# Generate color span for NES
col_fun = colorRamp2(c(min(NES_df, na.rm = T), 0, max(NES_df, na.rm = T)), c("steelblue3", "white", "firebrick3"))

# Define cluster colors
cluster_names = rownames(NES_df)
color_palette = hcl.colors(length(cluster_names), "Roma")
color_mapping = setNames(color_palette, cluster_names)
cluster_colors = color_mapping[cluster_names] 

# Generate annotation
row_anno = rowAnnotation(
  Cluster = cluster_names,
  col = list(
    Cluster = cluster_colors),
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

# Generate heatmap
hm = Heatmap(NES_df,
             right_annotation = row_anno,
             cluster_rows = F,
             cluster_columns = T,
             show_column_dend = F,
             show_row_names = T,
             name = "NES",
             row_title_gp = gpar(fontsize = 5),
             row_names_side = "right",
             row_names_gp = gpar(fontsize = 5),
             column_names_gp = gpar(fontsize = 5),
             heatmap_legend_param = list(direction = "vertical",
                                         labels_gp = gpar(fontsize = 5),
                                         title_gp = gpar(fontsize = 6, fontface = "bold"),
                                         title = "NES",
                                         legend_height = unit(7, "mm"),
                                         grid_width = unit(1, "mm")),
             border = T,
             border_gp = gpar(lwd = 0.5),
             col = col_fun,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%s", p_val_df[i, j]), x, y, gp = gpar(fontsize = 5))
             }
)
hm

# Save as pdf
pdf("meta_clusters_starcat/Fig.2c.pdf", width = 2.5, height = 2.5)
draw(hm, merge_legend = TRUE, heatmap_legend_side = "right", padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()

