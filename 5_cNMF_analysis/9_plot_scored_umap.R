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

seu = readRDS("annotated_seu.rds")

#### Make Figure 1 panels ####

p_umap = DimPlot(seu, group.by="seurat_clusters", 
                 reduction = "umap",
                 shuffle = T,
                 pt.size = 0.001,
                 raster = F) &
  scale_color_brewer(palette = "Set1") &
  scale_fill_brewer(palette = "Set1") &
  coord_fixed() &
  guides(color=guide_legend(ncol=7,
                            override.aes = list(size=2))) &
  theme_pubr(base_size = 5) &
  theme(axis.line = element_line(linewidth=0.4),
        legend.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.spacing.x = unit(0, 'mm'))

p_umap = ggrastr::rasterize(p_umap, dpi=1000, dev="ragg_png")
save_plot("Fig.1e.pdf", p_umap, base_height = 90, base_width = 90, units="mm")

fraction_df = seu@meta.data %>%
  dplyr::select(Sample, Dataset, Cohort, seurat_clusters) %>%
  group_by(Sample, Dataset, Cohort, seurat_clusters) %>%
  summarize(cell_count_per_cluster = n()) %>%
  group_by(Sample, Dataset, Cohort) %>%
  mutate(cell_count = sum(cell_count_per_cluster)) %>% 
  filter(cell_count>=20) %>%
  mutate(cluster_fraction = cell_count_per_cluster/cell_count) %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  mutate(Cohort = factor(Cohort, levels = c("Young", "Aged")))

stat.test = fraction_df %>%
  group_by(seurat_clusters) %>%
  wilcox_test(cluster_fraction ~ Cohort) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_xy_position(x = "seurat_clusters")

p_fraction = ggplot(data = fraction_df) +
  geom_boxplot(aes(y=cluster_fraction, x=seurat_clusters, fill=Cohort), outlier.shape = NA, linewidth=0.25) +
  geom_point(aes(y=cluster_fraction, x=seurat_clusters, color=Cohort), position = position_jitterdodge(jitter.width = 0.1), size=0.2) +
  scale_fill_manual(values = c("#c7d1ca", "#617a6a"),
                    guide = guide_legend(ncol = 1,
                                         title.position = "top",
                    )) +
  scale_color_manual(values = c("black", "black")) +
  theme_pubr(base_size = 5) +
  stat_pvalue_manual(size=5/2.83465, stat.test %>% filter(p.adj<0.05), label = "p = {scales::pvalue(p.adj)}", bracket.size = 0.25) +
  scale_y_continuous(limits = c(-0.01,0.7),
                     name = "Fraction of predominant program",
                     breaks = c(0,0.25,0.5),
                     expand = c(0,0)) +
  theme(axis.line = element_line(linewidth=0.4),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5))

save_plot("Fig.1f.pdf", p_fraction, base_height = 45, base_width = 50, units="mm")

# Due to the unexpected interaction of scCustomize and ggrastr, I have to add the ninth signature to the original plot object; otherwise, the 8th one is getting dropped later on..

signatures = c("Aging_signature", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_DNA_REPLICATION", "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION")

p = scCustomize::FeaturePlot_scCustom(seu, 
                                      signatures, 
                                      na_cutoff = NA, 
                                      order = F,
                                      num_columns = 3,
                                      pt.size = 0.001) & 
  coord_fixed() & 
  theme_pubr(base_size = 5) &
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        legend.justification = "left",
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=5)) & 
  scale_color_gradientn(
    colours = brewer.pal(9, "YlGnBu"),
    breaks = seq(-0.4, 1.2, by = 0.2)  # Uniform break intervals
  )

plots = p$patches$plots
plots = lapply(plots, function (x) { ggrastr::rasterize(x, dpi=1000, dev="ragg_png") })
p = patchwork::wrap_plots(plots, ncol = 3)
save_plot("Fig.1d.pdf", p, base_height = 120, base_width = 280, units="mm")

#### Make Suppl. Fig. 5 ####

signatures = c("Aging_signature", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Garcia_quiescence", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_DNA_REPLICATION", "REACTOME_G2_M_CHECKPOINTS", "HALLMARK_E2F_TARGETS", "REACTOME_CELL_CYCLE_MITOTIC", "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION")

p = scCustomize::FeaturePlot_scCustom(seu, 
                                      signatures, 
                                      na_cutoff = NA, 
                                      order = F,
                                      num_columns = 4) & 
  coord_fixed() & 
  theme_pubr(base_size = 5) &
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        legend.justification = "left",
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=5)) & 
  scale_color_gradientn(
    colours = brewer.pal(9, "YlGnBu"),
    breaks = seq(-0.4, 1.2, by = 0.2)  # Uniform break intervals
  )

plots = p$patches$plots
plots = lapply(plots, function (x) { ggrastr::rasterize(x, dpi=1000, dev="ragg_png") })
p_a = patchwork::wrap_plots(plots, ncol = 4)

p_b = seu@meta.data %>%
  dplyr::select(Dataset, Sample, Cohort, seurat_clusters, all_of(signatures[1:8])) %>%
  pivot_longer(-c(Dataset, Sample, Cohort, seurat_clusters)) %>%
  group_by(Sample, Dataset) %>%
  filter(n()>=20) %>%
  group_by(Sample, Dataset, seurat_clusters, name) %>%
  summarize(m = mean(value)) %>%
  mutate(name = factor(name, levels = signatures[1:8])) %>%
  ggplot(aes(x=seurat_clusters,  y=m, fill=seurat_clusters)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(linewidth=0.25, outlier.size = 0.3) +
  ylab("Average signature score per sample") +
  facet_wrap(.~name, ncol=4, scales = "free") + 
  theme_pubr(base_size = 5) +
  theme(axis.line = element_line(linewidth=0.4),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5))

seu@meta.data$young_cohort_frac = seu@meta.data %>%
  group_by(seurat_clusters) %>%
  mutate(frac = sum(Cohort == "Young")/ sum(Cohort %in% c("Young", "Aged"))) %>%
  pull(frac)

p_c = scCustomize::FeaturePlot_scCustom(BuenColors::shuf(seu), features = "young_cohort_frac", pt.size = 0.05, order = F) +
  coord_fixed() +
  ggtitle("Fraction of cells from Young\ndonors, per cluster") +
  theme_pubr(base_size = 5) +
  theme(axis.line = element_line(linewidth=0.4, color="black"),
        strip.background = element_rect(color="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        axis.ticks = element_blank(),
        # axis.text = element_blank(),
        legend.text = element_text(size=5)) +
  scale_color_gradientn(
    colours = brewer.pal(9, "YlGnBu"),
    limits = c(0,1)  # Uniform break intervals
  )
p_c = ggrastr::rasterize(p_c, dpi=1000, dev="ragg_png")

p_d = DimPlot(seu %>% subset(Cohort %in% c("Young", "Aged")), 
              group.by="Cohort",
              reduction = "umap",
              shuffle = T,
              pt.size = 0.001,
              raster = F) &
  scale_color_brewer(palette = "Set1", breaks = c("Young", "Aged")) &
  coord_fixed() &
  guides(color=guide_legend(ncol=7,
                            override.aes = list(size=2))) &
  theme_pubr(base_size = 5) &
  theme(axis.line = element_line(linewidth=0.4),
        legend.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        axis.ticks = element_blank(),
        # axis.text = element_blank(),
        legend.spacing.x = unit(0, 'mm'))
p_d = ggrastr::rasterize(p_d, dpi=1000, dev="ragg_png")

p_combined = ggarrange(p_a, ggarrange(p_b, ggarrange(p_c, p_d, ncol=1, labels = c("c", "d"), font.label = list(size=6)), ncol=2, widths = c(0.7, 0.2)), ncol=1, labels = c("a", "b"), font.label = list(size=6))

save_plot("Suppl.Fig.5.pdf", p_combined, base_height = 200, base_width = 180, units="mm")

#### Make boxplot with fractions of young/aged cells in each cluster

fraction_df = seu@meta.data %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  dplyr::select(Sample, Dataset, Cohort, seurat_clusters) %>%
  group_by(Dataset, Cohort, seurat_clusters) %>%
  summarize(cell_count_per_cluster = n()) %>%
  group_by(Dataset, seurat_clusters) %>%
  mutate(cohort_fraction = cell_count_per_cluster/sum(cell_count_per_cluster)) %>%
  mutate(Cohort = factor(Cohort, levels = c("Young", "Aged")))

stat.test = fraction_df %>%
  group_by(seurat_clusters) %>%
  wilcox_test(cohort_fraction ~ Cohort) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_xy_position(x = "seurat_clusters")

colors = brewer.pal(7, "Set1")
colors[6] = "gold"
p_fraction = ggplot(data = fraction_df) +
  geom_boxplot(aes(y=cohort_fraction, x=seurat_clusters, fill=Cohort), outlier.shape = NA, linewidth=0.25) +
  geom_point(aes(y=cohort_fraction, x=seurat_clusters, group=Cohort, color=Dataset), position = position_jitterdodge(jitter.width = 0.1), size=0.2) +
  scale_fill_manual(values = c("#c7d1ca", "#617a6a"),
                    guide = guide_legend(ncol = 1,
                                         title.position = "top",
                    )) +
  scale_color_manual(values = colors,
                     guide = guide_legend(ncol = 3,
                                          title.position = "top",
                     )) +
  theme_pubr(base_size = 5) +
  stat_pvalue_manual(size=5/2.83465, stat.test %>% filter(p.adj<0.05), label = "p = {scales::pvalue(p.adj)}", bracket.size = 0.25) +
  scale_y_continuous(name = "Fraction of cohort",
                     expand = c(0.05,0.05)) +
  theme(axis.line = element_line(linewidth=0.4),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5))
p_fraction

#### Make Suppl. Fig. 8 ####

signatures = c("metaGEP.1", "metaGEP.2", "metaGEP.3", "metaGEP.4")
p = scCustomize::FeaturePlot_scCustom(seu, 
                                      signatures, 
                                      na_cutoff = NA, 
                                      order = F,
                                      num_columns = 2, 
                                      raster = T,
                                      raster.dpi = c(2000,2000), 
                                      pt.size = 5) & 
  coord_fixed() & 
  theme_pubr(base_size = 5) &
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        legend.justification = "left",
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size=5)) & 
  scale_color_gradientn(
    colours = brewer.pal(9, "YlGnBu"),
    breaks = seq(-0.4, 1.2, by = 0.2)  # Uniform break intervals
  )
p
save_plot("Suppl.Fig.8.pdf", p, base_height = 200, base_width = 180, units="mm")

