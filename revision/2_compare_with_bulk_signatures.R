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
setwd(paste0(base_dir, "revision"))

data_dir = paste0(base_dir, "revision/bulk_studies_results/")

# Load DE results of this study
de_results = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep=",", header = T)

# Load DE results from previous studies
other_studies_results = list()

# Adelman2019 (bulk); PMID: 31085557; genes are extracted from Supplementary Table 6
adelman_results = read.csv(paste0(data_dir, "Adelman2019_results.csv"))
aging_up = adelman_results %>%
  filter(padj<0.05, log2FoldChange > 0.5) %>%
  pull(GeneSym)
aging_down = adelman_results %>%
  filter(padj<0.05, log2FoldChange < -0.5) %>%
  pull(GeneSym)
other_studies_results[["Adelman2019"]] = list(aging_up, aging_down)

# Pang2011; PMID: 22123971; genes are extracted from Appendix (PDF)
aging_up = read.table(paste0(data_dir, "Pang2011_genes_up.txt"))$V1
aging_down = read.table(paste0(data_dir, "Pang2011_genes_down.txt"))$V1
other_studies_results[["Pang2011"]] = list(aging_up, aging_down)

# Nilsson2016; PMID: 27368054; genes are extracted from Supplementary Table 1
results = read.csv(paste0(data_dir, "Nilsson2016_results.csv"))
aging_up = results$Upregulated
aging_down = results$Downregulated
other_studies_results[["Nilsson2016"]] = list(aging_up, aging_down)

# Hennrich2018; PMID: 30275468; genes are extracted from Supplementary Data 5
# Hennrich2018 contains results for HSPCs, not HSCs
results = read.csv(paste0(data_dir, "Hennrich2018_results.csv"))
aging_up = results %>%
  filter(padj<0.05, log2FoldChange > 0.5) %>%
  pull(gene)
aging_down = results %>%
  filter(padj<0.05, log2FoldChange < -0.5) %>%
  pull(gene)
other_studies_results[["Hennrich2018"]] = list(aging_up, aging_down)

# Make volcano plots, highlighting consistently up- or down-regulated genes for each prior study
# Create a directory if it doesn't exist yet
if (!dir.exists("compare_to_bulk_results")) dir.create("compare_to_bulk_results")

for (Dataset in c("Adelman2019", "Pang2011", "Nilsson2016", "Hennrich2018")) {
  
  aging_up = other_studies_results[[Dataset]][[1]]
  aging_down = other_studies_results[[Dataset]][[2]]
  
  p = de_results %>%
    mutate(color = case_when(
      gene %in% aging_up ~ "up",
      gene %in% aging_down ~ "down",
      .default = "unchanged")) %>%
    mutate(label = ifelse( abs(log2FoldChange) > 0.5 & color %in% c("up", "down") & padj < 0.05, gene, "")) %>%
    mutate(color = factor(color, levels = c("unchanged", "down", "up"))) %>%
    arrange(color) %>%
    ggplot(aes(x=log2FoldChange, y=-log10(padj), color=color, label=label)) +
    geom_vline(xintercept = 0.5, colour="grey", linetype="dashed", linewidth=0.2) +
    geom_vline(xintercept = -0.5, colour="grey", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = -log10(0.05), colour="grey", linetype="dashed", linewidth=0.2) +
    geom_point(size=0.5) +
    ggrepel::geom_text_repel(color="black", max.overlaps = 200, size=3, min.segment.length = 0) +
    scale_color_manual(values = c("grey80", "steelblue2", "firebrick2")) +
    theme_pubr() +
    ggtitle(Dataset)
  
  # Individual panels of Suppl. Fig. 3:
  save_plot(paste0("compare_to_bulk_results/volcano.", Dataset, "_de_results_highlighted.png"), p, base_height = 6, base_width = 6, dpi=600)
  
}

# Test whether age-associated signatures from prior studies are enriched in this study's DE results
pathways_for_gsea = list()

for (Dataset in c("Adelman2019", "Pang2011", "Nilsson2016", "Hennrich2018")) {
  pathways_for_gsea[[paste0(Dataset, "_UP")]] = other_studies_results[[Dataset]][[1]]
  pathways_for_gsea[[paste0(Dataset, "_DOWN")]] = other_studies_results[[Dataset]][[2]]
}

ranks = as.numeric(de_results$log2FoldChange)
names(ranks) = de_results$gene

# Run GSEA
fgseaRes = fgsea(pathways = pathways_for_gsea, 
                 stats    = ranks,
                 minSize  = 5,
                 eps      = 0.0,
                 maxSize  = 1000,
                 nPermSimple = 10000)

p = fgseaRes %>%
  ggplot(aes(x=pathway, y=NES, fill=padj, label=round(padj,3), color=padj<0.05)) +
  geom_bar(stat = "identity", linewidth=1) +
  scale_fill_gradient(low="steelblue2", high="grey90") +
  theme_pubr() +
  scale_color_manual(values = c("white", "black"), guide="none") +
  coord_flip() +
  xlab("Gene set")
save_plot("compare_to_bulk_results/Suppl.Fig.4.png", p, base_height = 6, base_width = 6, dpi = 600)
