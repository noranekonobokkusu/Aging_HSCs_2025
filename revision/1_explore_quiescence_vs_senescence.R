suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ggpubr)
  library(fgsea)
  library(cowplot)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"

setwd(paste0(base_dir, "revision"))

#### Address senescence-related questions ####

# Read DE results
de_result = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep = ",", header = T)

# Read all potentially relevant pathways
# Human dormancy signature from https://doi.org/10.1016/j.stem.2021.07.003
garcia_table = readxl::read_xlsx(paste0(base_dir, "input_data/signatures/mmc2.xlsx"), sheet = 1)
garcia_signature = garcia_table %>%
  filter(genename!="None") %>%
  filter(logFC < -2,FDR<0.00001) %>%
  pull(genename) %>%
  unique()

graham_signature = read.table(paste0(base_dir, "input_data/signatures/graham_quiescence_up.txt"), header = F)$V1

hallmark_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/h.all.v2024.1.Hs.symbols.gmt.txt"))
c2_cp_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/c2.cp.v2024.1.Hs.symbols.gmt.txt"))

# Add two senescence signatures outside HALLMARK and C2:CP sets
saul_sen = gmtPathways(paste0(base_dir, "input_data/signatures/SAUL_SEN_MAYO.v2025.1.Hs.gmt"))[[1]]
fridman_sen = gmtPathways(paste0(base_dir, "input_data/signatures/FRIDMAN_SENESCENCE_UP.v2025.1.Hs.gmt"))[[1]]

pathways = list(garcia_signature, graham_signature, saul_sen, fridman_sen)

names(pathways) = c("GARCIA_QUIESCENCE", "GRAHAM_QUIESCENCE", "SAUL_SENESCENCE", "FRIDMAN_SENESCENCE")

pathways = c(pathways, hallmark_pathways, c2_cp_pathways)

# We settled on six different gene sets that sounded diverse and potentially relevant
senescence_pathways = pathways[grep("SENESCENCE", names(pathways), value=T)[c(1,2,3,6,7,8)]]
names(senescence_pathways)

# Are these senescence gene sets consistent between each other? 
common_genes = reduce(senescence_pathways, intersect) # 0 genes are shared across all six gene sets

# We defined two "consensus" senescence-realted gene sets shared across at least 4 (16 genes) or 3 (97 genes) gene sets
shared_senescence_genes_four = stack(senescence_pathways) %>%
  group_by(values) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>%
  filter(n>=4) %>%
  pull(values)

shared_senescence_genes_three = stack(senescence_pathways) %>%
  group_by(values) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>%
  filter(n>=3) %>%
  pull(values)

# We need to add them to the rest
pathways[["SHARED_SENESCENCE_GENES_FOUR"]] = shared_senescence_genes_four

pathways[["SHARED_SENESCENCE_GENES_THREE"]] = shared_senescence_genes_three

# Rename a pathway with a very long name
pathways[["REACTOME_SASP"]] = pathways[["REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"]]
pathways[["REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"]] = NULL

# Plot volcano plots

# Define a function that would highlight and label genes from a gene set on a volcano plot
plot_volcano = function(de_result, pathways, gene_set_name, dir_name) {
  gene_set = pathways[[gene_set_name]]
  # Set thresholds and mark significant hits
  padj_cutoff = 0.05
  log2fc_cutoff = 1 # (100% increase)
  
  genes_to_label = de_result %>%
    filter(padj<padj_cutoff & (gene %in% gene_set) & abs(log2FoldChange)>log2fc_cutoff) %>%
    pull(gene)
  
  p = de_result %>%
    mutate(signif_hit = ifelse(padj<padj_cutoff & abs(log2FoldChange)>log2fc_cutoff, 1, 0),
           is_labeled = ifelse(gene %in% gene_set, 1, 0)) %>%
    arrange(signif_hit) %>%
    ggplot() +
    geom_vline(xintercept = 1, colour="grey", linetype="dashed", linewidth=0.2) +
    geom_vline(xintercept = -1, colour="grey", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept = -log10(0.05), colour="grey", linetype="dashed", linewidth=0.2) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), color = (signif_hit==1)), size=0.5) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), shape = factor(is_labeled)), color = "red", size=0.85) +
    ggrepel::geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(gene %in% genes_to_label, gene, "")), force=5, max.overlaps = 300, direction="both", box.padding = 0.4, min.segment.length = 0, segment.size=0.1, point.padding = 0.1, size = 5/2.845) + 
    scale_shape_manual(
      name = "AP-1 complex",   # Renames the legend title
      values = c("0" = NA, "1" = 16),  # Assigns a valid shape to "1" and hides "0"
      breaks = c("1")  # Ensures only "1" appears in the legend
    ) +
    theme_pubr(base_size = 4) +
    scale_color_manual(values = c("grey", "steelblue"), name = "Significant hit") +
    scale_y_continuous(expand = c(0.01,0.01), name =  expression(-log[10] ~ P[adj])) +
    scale_x_continuous(name =  expression(log[2] ~ "Fold Change")) +
    theme(axis.ticks = element_line(color = "black"),
          axis.line = element_line(linewidth=0.3),
          axis.title = element_text(size=5),
          legend.position = "none",
          title = element_text(size=3.5, face = "bold")) +
    ggtitle(gene_set_name)
  save_plot(paste0(dir_name, "/", Sys.time(), "_", gene_set_name, ".png"), p, base_height = 2, base_width = 2, dpi=600)
}

# In addition to quiescence, senescence, p53 signaling, we also looked at two pathways that potentially indicate stress response - ROS and DNA damage. Additionally, we included JAK-STAT signaling to address Reviewer 2's question
pathways_of_interest = c("SAUL_SENESCENCE", "FRIDMAN_SENESCENCE", "REACTOME_CELLULAR_SENESCENCE", "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE", "REACTOME_ONCOGENE_INDUCED_SENESCENCE", "REACTOME_SASP", "SHARED_SENESCENCE_GENES_FOUR", "SHARED_SENESCENCE_GENES_THREE", "GARCIA_QUIESCENCE", "GRAHAM_QUIESCENCE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "WP_DNA_DAMAGE_RESPONSE", "HALLMARK_P53_PATHWAY", "KEGG_P53_SIGNALING_PATHWAY", "KEGG_JAK_STAT_SIGNALING_PATHWAY")

# Create a directory if it doesn't exist yet
if (!dir.exists("senescence_analysis_results")) dir.create("senescence_analysis_results")

# Generate annotated volcano plots
for (gene_set_name in pathways_of_interest) {
  plot_volcano(de_result, pathways, gene_set_name, "senescence_analysis_results")
}

# Run GSEA and visualize the selected pathways
ranks = as.numeric(de_result$log2FoldChange)
names(ranks) = de_result$gene

fgseaRes = fgsea(pathways = pathways, 
                 stats    = ranks,
                 minSize  = 5,
                 eps      = 0.0,
                 maxSize  = 500,
                 nPermSimple = 10000)

fgseaRes = fgseaRes %>% arrange(padj)

p1 = fgseaRes %>%
  filter(pathway %in% pathways_of_interest) %>%
  dplyr::select(ID=pathway, padj, NES) %>%
  mutate(padj_log = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(x=NES, y=ID, fill = padj_log, label=round(padj, 4))) +
  geom_bar(stat = "identity") +
  geom_text(aes(fontface = ifelse(padj<0.05, "bold", "plain")), nudge_x = -0.6, size = 5/2.845) +
  theme_pubr(base_size = 5) +
  scale_fill_gradient(low = "powderblue", high = "steelblue", trans="log", breaks = scales::pretty_breaks(), name =  expression(-log[10] ~ P[adj])) +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_line(linewidth=0.3), legend.direction = "vertical",
        legend.position = "right",  # Moves the legend above the plot
        legend.justification = "left",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.3, "cm"))
p1 
save_plot("senescence_analysis_results/gsea_pathways_of_interest.png", p1, base_height = 3, base_width = 4, dpi=600)

# p21 (CDKN1A) is moderately significantly upregulated in Aged samples and shows up in several gene sets. How concordant is CDKN1A expression?

seu = readRDS(paste0(base_dir, "5_cNMF_analysis/annotated_seu.rds"))

meta = seu@meta.data
gene_expr = seu@assays$RNA$data["CDKN1A",]
all(names(gene_expr) == rownames(meta)) 
meta$CDKN1A = gene_expr

p = meta %>%
  group_by(Sample, Cohort, Dataset) %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  mutate(Dataset = ifelse(Dataset=="Hourigan2018", "Oetjen2018", Dataset)) %>%
  summarize(n=n(),
            mean_expr = mean(CDKN1A)) %>%
  mutate(sample_used_for_DE = ifelse(n>=20, 1,0)) %>%
  mutate(Cohort = factor(Cohort, levels = c("Cord", "Prenatal", "Infant", "Child", "Young", "Middle", "Aged"))) %>%
  ggplot(aes(x=Dataset,y=mean_expr,fill=Cohort,stroke=sample_used_for_DE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_brewer(palette = "Set1") +
  scale_size_manual(values = c(0.1, 2)) +
  ylab("Average normalized expression, per sample")
p
save_plot("senescence_analysis_results/CDKN1A_average_expression.png", p, base_height = 5, base_width = 6, dpi=600)

# We also assessed enrichment of senescence pathways in metaGEPs

# Read metaGEPs
programs = read.table(paste0(base_dir, "5_cNMF_analysis/meta_clusters_starcat/metaGEP.Z_scores.tsv"), row.names = 1, header = T, sep="\t")

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

p = fgseaRes_total %>% 
  filter(pathway %in% pathways_of_interest) %>%
  mutate(sign = case_when(padj < 0.001 ~ "***",
                          padj < 0.01 ~ "**",
                          padj < 0.1 ~ "*",
                          .default = "")) %>%
  ggplot(aes(x=Program, y=pathway, fill=NES, label=sign)) +
  geom_tile() +
  geom_text(size = 5/2.845) +
  theme_pubr(base_size = 5) +
  scale_fill_gradient2(low = "steelblue2", high = "firebrick2") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
p
save_plot("senescence_analysis_results/metaGEP_pathways_GSEA.png", p, base_height = 4, base_width = 5, dpi=600)

# Plot ranked metaGEPs with gene set members highlighted
for ( pathway_of_interest in pathways_of_interest) {
  p = programs %>%
    rownames_to_column("Program") %>%
    pivot_longer(-Program) %>%
    group_by(Program) %>%
    arrange(desc(value), .by_group = TRUE) %>%
    mutate(r = row_number()) %>%
    mutate(is_in_poi = ifelse(name %in% pathways[[pathway_of_interest]], TRUE, FALSE)) %>%
    arrange(is_in_poi = factor(is_in_poi, levels = c(FALSE, TRUE))) %>%
    ggplot(aes(x=r,y=value, label=ifelse(is_in_poi, name, ""), color=is_in_poi)) +
    geom_point(size=0.2) +
    scale_color_manual(values = c("powderblue", "firebrick2")) +
    ggrepel::geom_text_repel(max.overlaps = 150, show.legend = FALSE, size = 5/2.83465, segment.size=0.2) +
    facet_wrap(.~Program, nrow = 1) +
    theme_pubr(base_size = 5) +
    theme(axis.line = element_line(linewidth=0.4),
          strip.background = element_rect(colour="black", fill="white", 
                                          linewidth=0.8, linetype="solid"),
          strip.text = element_text(size=5),
          legend.key.height = unit(0.25, "cm"),
          legend.spacing.y = unit(0, "mm"),
          legend.box.spacing = unit(0, "mm"),
          legend.box = "horizontal",
          legend.margin = margin(0,0,0,0),
          legend.text = element_text(size=5)) +
    ggtitle(pathway_of_interest)
  
  save_plot(paste0("senescence_analysis_results/", Sys.time(), "_", pathway_of_interest, ".png"), p, base_width = 120, base_height = 30, units = "mm", dpi=600)
  
}
