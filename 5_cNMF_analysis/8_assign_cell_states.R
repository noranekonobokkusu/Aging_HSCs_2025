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
setwd(paste0(base_dir, "5_cNMF_analysis"))

#### Read the master seurat object ####
seu_total = readRDS(paste0(base_dir, "1_dataset_preparation/seu_hsc.total.rds"))

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

mhc_ii = c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "CD74")

hallmark_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/h.all.v2024.1.Hs.symbols.gmt.txt"))
c2_cp_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/c2.cp.v2024.1.Hs.symbols.gmt.txt"))

pathways = list(graham_signature, garcia_signature, aging_sig, mep_sig, gmp_sig, clp_sig, mhc_ii, c2_cp_pathways$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION, hallmark_pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB, hallmark_pathways$HALLMARK_INTERFERON_ALPHA_RESPONSE, hallmark_pathways$HALLMARK_INTERFERON_GAMMA_RESPONSE,
                hallmark_pathways$HALLMARK_E2F_TARGETS,c2_cp_pathways$REACTOME_CELL_CYCLE, c2_cp_pathways$REACTOME_DNA_REPLICATION, c2_cp_pathways$REACTOME_G2_M_CHECKPOINTS, c2_cp_pathways$REACTOME_CELL_CYCLE_MITOTIC)

names(pathways) = c("Graham_quiescence", "Garcia_quiescence", "Aging_signature", "MEP_signature", "GMP_signature", "CLP_signature", "MHC-II", "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "REACTOME_CELL_CYCLE", "REACTOME_DNA_REPLICATION", "REACTOME_G2_M_CHECKPOINTS", "REACTOME_CELL_CYCLE_MITOTIC")

#### Score cells with AddModuleScore ####

seu_total = NormalizeData(seu_total)
seu_total = AddModuleScore(seu_total, pathways, search = T)
colnames(seu_total@meta.data)[(ncol(seu_total@meta.data)-length(pathways) + 1):ncol(seu_total@meta.data)] = names(pathways)
colnames(seu_total@meta.data) # Check

# Read meta-program usages and identify the predominant program

usage_files = list.files(paste0("starcat_usages"), pattern = paste0(".starcat_usage.csv$"), recursive = TRUE, full.names = T)

usage_df = data.frame()
for (i in 1:length(usage_files)) {
  t = read.table(usage_files[i], sep=",", header = T, row.names = 1)
  usage_df = rbind(usage_df, t)
}

programs = colnames(usage_df)

usage_df = usage_df %>%
  rownames_to_column("Cell") %>%
  pivot_longer(-Cell) %>%
  group_by(Cell) %>%
  mutate(Predominant_program = name[which.max(value)]) %>%
  pivot_wider(values_from = value, names_from = name)

# Read UMAP coordinates from scVI-corrected object
seu_scVI = MuDataSeurat::ReadH5AD(paste0(base_dir, "3_scVI_integration/integrated_seu.h5ad"))
seu_scVI@meta.data[c("umap_1", "umap_2")] = NULL
umap_coords = Embeddings(seu_scVI, reduction = "umap")
scvi_coords = Embeddings(seu_scVI, reduction = "scVI")

# Create a copy of the master seurat object and add program usages and UMAP coordinates to it

seu = seu_total
usage_df = usage_df[match(colnames(seu), usage_df$Cell),]
all(colnames(seu) == usage_df$Cell)
all(colnames(seu) == rownames(umap_coords))
usage_df = cbind(usage_df, Leiden_cluster = seu_scVI$leiden_0.5)
seu = AddMetaData(seu, usage_df)
seu[["umap_python"]] = CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(seu))
seu[["scVI"]] = CreateDimReducObject(embeddings = scvi_coords, key = "SCVI_", assay = DefaultAssay(seu))
seu = FindNeighbors(seu, reduction = "scVI", dims = 1:30)
seu = FindClusters(seu, resolution = 0.5)
seu = RunUMAP(seu, reduction = "scVI", dims = 1:30, reduction.name = "umap")

# Save final Seurat object
saveRDS(seu, "annotated_seu.rds")

# Rename Hourigan and Safina datasets
seu@meta.data$Dataset = case_when(seu@meta.data$Dataset == "Hourigan2018" ~ "Oetjen2018",
                                  seu@meta.data$Dataset == "Safina2024" ~ "Safina2025",
                        .default = seu@meta.data$Dataset)

# Plot UMAP with cells colored by the predominant program

p_umap = DimPlot(seu %>% subset(Cohort %in% c("Aged", "Young")),
                    group.by = "Predominant_program",
                    reduction = "umap",
                    shuffle = T,
                    pt.size = 0.05,
                    raster = F) &
  scale_color_manual(values = hcl.colors(length(programs), "Roma")) &
  coord_fixed() &
  guides(color=guide_legend(ncol=2,
                            override.aes = list(size=2))) &
  theme_pubr(base_size = 5) &
  theme(axis.line = element_line(linewidth=0.4),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"))

p_umap = ggrastr::rasterize(p_umap, dpi=1000, dev="ragg_png")
save_plot("Fig.2d.pdf", p_umap, base_height = 72.72, base_width = 72.72, units="mm")


# Plot and compare program fractions between cohorts

fraction_df = seu@meta.data %>%
  dplyr::select(Sample, Dataset, Cohort, Predominant_program) %>%
  group_by(Sample, Dataset, Cohort, Predominant_program) %>%
  summarize(cell_count_per_program = n()) %>%
  group_by(Sample, Dataset, Cohort) %>%
  mutate(cell_count = sum(cell_count_per_program)) %>% 
  filter(cell_count>=20) %>%
  mutate(program_fraction = cell_count_per_program/cell_count) %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  mutate(Cohort = factor(Cohort, levels = c("Young", "Aged")))

stat.test = fraction_df %>%
  group_by(Predominant_program) %>%
  wilcox_test(program_fraction ~ Cohort) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_xy_position(x = "Cohort") 

color_palette = brewer.pal(12, "Paired")[c(4,2,8,12,9,1,6)]
color_palette[3] = "#ffcb00"

p_fractions = ggplot(data = fraction_df, 
                        aes(y=program_fraction, x=Cohort)) +
  geom_boxplot(aes(fill=Cohort), outlier.shape = NA, linewidth=0.25) +
  scale_fill_manual(values = c("#c7d1ca", "#617a6a"),
                    guide = guide_legend(ncol = 1,
                                         title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill=Dataset), position=position_jitterdodge(dodge.width = 0.1), colour="black", pch=21, size=1,stroke=0.25) +
  facet_wrap(Predominant_program~., nrow=1) +
  scale_fill_manual(values = color_palette, 
                    guide = guide_legend(ncol = 2, 
                                         title.position = "top",
                                         override.aes = list(size=2))) +
  theme_pubr(base_size = 5) +
  stat_pvalue_manual(size=5/2.83465, stat.test, label = "p = {scales::pvalue(p.adj)}", bracket.size = 0.25) +
  scale_y_continuous(limits = c(-0.02,1.15), 
                     name = "Fraction of predominant program",
                     breaks = c(0,0.25,0.5,0.75,1),
                     expand = c(0,0)) +
  theme(axis.line = element_line(linewidth=0.4, color="black"),
        axis.ticks = element_line(color="black"),
        strip.background = element_rect(color="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5))

# Save Fig. 2e
save_plot("Fig.2e.pdf", p_fractions, base_height = 50, base_width = 63.63636, units="mm")

