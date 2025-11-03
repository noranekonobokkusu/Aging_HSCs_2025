suppressMessages(suppressWarnings({
  library(tidyverse)
  library(Seurat)
  library(cowplot)
  library(MuDataSeurat)
  library(RColorBrewer)
  library(ggpubr)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "3_scVI_integration/"))

seu = ReadH5AD("hspc/integrated_seu.h5ad")
seu@meta.data[c("umap_1", "umap_2")] = NULL
seu = NormalizeData(seu)
seu = FindNeighbors(seu, reduction = "scVI", dims = 1:30)
seu = FindClusters(seu, resolution = 1.5)
seu = RunUMAP(seu, reduction = "scVI", dims = 1:30, reduction.name = "umap_seurat")

# Explore
DimPlot(seu, group.by = "Dataset", reduction = "umap_seurat") & scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired")))

DimPlot(seu, split.by = "Dataset", ncol=3, reduction = "umap_seurat")

DimPlot(seu, group.by = "seurat_clusters", label=T, reduction = "umap_seurat") & scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"))) & NoLegend()

# Add HSC signature scores
signatures = read.csv(paste0(base_dir, "input_data/signatures/hsc_signatures.csv"))
sig_list = list(signatures$laurenti_HSC, signatures$eppert_CE_HSC_LSC, signatures$jaatinen_HSC_UP)
seu = AddModuleScore(seu, features = sig_list, search = T)
colnames(seu@meta.data)[(ncol(seu@meta.data)-2):ncol(seu@meta.data)] = c("Laurenti", "Eppert", "Jaatinen")

# Plot UMAP
p1 = DimPlot(seu, group.by = "seurat_clusters", label=T, reduction = "umap_seurat", raster = T, raster.dpi = c(1000,1000), label.size = 5/2.83465) + 
  scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"))) +
  NoLegend() + 
  coord_fixed() +
  theme_pubr(base_size = 5) +
  theme(axis.line = element_line(linewidth=0.4),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5), 
        title = element_text(size=6),
        axis.text = element_blank())

# Plot scores
d = seu@meta.data %>%
  dplyr::select(seurat_clusters, Laurenti, Eppert, Jaatinen) %>%
  tidyr::pivot_longer(-seurat_clusters) %>%
  group_by(seurat_clusters, name) %>%
  summarize(m = mean(value)) %>%
  group_by(seurat_clusters) %>%
  mutate(mean_m = mean(m)) %>%
  arrange(desc(mean_m)) 

d$seurat_clusters = factor(d$seurat_clusters, levels = unique(d$seurat_clusters))

p2 = ggplot(d, aes(y=as.factor(seurat_clusters), x=name, fill=m)) +
  geom_tile() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_distiller(palette = "Spectral", 
                       name = "Average score", 
                       guide = guide_colorbar(title.position = "top", barwidth = unit(30, "mm")),
                       breaks = seq(-0.4, 1, by = 0.2)) +
  theme_pubr(base_size = 5) +
  theme(axis.line = element_line(linewidth=0.4),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5), 
        title = element_text(size=6),
        axis.title.x = element_blank()) +
  ylab("Seurat cluster")

# Clusters 26, 3, 9, 0 look like good HSC candidates; 44187 cells in total
seu@meta.data %>%
  filter(seurat_clusters %in% c(26, 3, 9, 0)) %>%
  nrow()

# THe selected clusters include 87% of HSCs annotated by BMM
seu@meta.data %>%
  filter(predicted_CellType=="HSC") %>%
  mutate(hsc_annotated = seurat_clusters %in% c(26, 3, 9, 0)) %>%
  summarize(n=sum(hsc_annotated)/n())

# Plot the fraction of cells annotated as HSCs by BMM, per cluster
d2 = seu@meta.data %>%
  dplyr::select(seurat_clusters, is_hsc) %>%
  group_by(seurat_clusters) %>%
  summarize(n = n(),
            freq = sum(is_hsc=="HSC")/n())

d2$seurat_clusters = factor(d2$seurat_clusters, levels = unique(d$seurat_clusters))

p3 = ggplot(d2, aes(y=as.factor(seurat_clusters), x=1, fill=freq, label=n)) +
  geom_tile() +
  geom_text(size=5/2.83465) +
  scale_x_discrete(expand = c(0,0), name = "") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_distiller(palette = "Spectral", 
                       name = "Fraction", 
                       guide = guide_colorbar(title.position = "top", barwidth = unit(15, "mm")),
                       breaks = seq(-0.4, 1, by = 0.2)) +
  theme_pubr(base_size = 5) +
  theme(axis.line = element_line(linewidth=0.4),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5), 
        axis.text.x = element_blank(),
        title = element_text(size=6)) +
  ylab("Seurat cluster")

# Combine three panels
p_combined = ggarrange(p1, p2, p3, widths = c(0.35, 0.3, 0.1), nrow=1, labels = c("a", "b", "c"), font.label = list(size=8))
save_plot("Suppl.Fig.1.pdf", p_combined, base_height = 100, base_width = 180, units="mm")
