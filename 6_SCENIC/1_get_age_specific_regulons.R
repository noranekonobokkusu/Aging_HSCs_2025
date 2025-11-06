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

# Read AUC matrices, scale them, and subset to shared regulons
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

# Rename Hourigan and Safina datasets
auc_long_df$Dataset = case_when(auc_long_df$Dataset == "Hourigan2018" ~ "Oetjen2018",
                                auc_long_df$Dataset == "Safina2024" ~ "Safina2025",
                                  .default = auc_long_df$Dataset)

#### Identify and visualize top regulons in Young and Aged cohorts ####

# Get top regulons in Young and Aged cohorts
regulons_mean_ranks_and_cors = auc_long_df %>%
  mutate(Regulon = gsub("\\..*", "", Regulon),
         Cohort_binary = ifelse(Cohort == "Young", 0, 1)) %>% 
  filter(Cohort_binary %in% c(0,1)) %>%
  group_by(Dataset, Regulon) %>%
  summarize(cor_value = cor(Regulon_AUC_Z_score, Cohort_binary)) %>%
  group_by(Dataset) %>%
  arrange(desc(cor_value), .by_group = T) %>%
  mutate(rank = row_number()) %>% 
  group_by(Regulon) %>%
  summarize(mean_rank = mean(rank),
            mean_cor = mean(cor_value)) %>%
  arrange(mean_rank, .by_group = T) %>%
  mutate(final_rank = row_number()) %>%
  arrange(final_rank)

# Get regulons that have high enough correlation with cohort status
top_regulons = regulons_mean_ranks_and_cors %>%
  filter(abs(mean_cor) >= 0.15 & final_rank <=10) %>% 
  pull(Regulon)

# Get top-10 and bottom-10 regulons to illustrate high variability of data
top_10_regulons = regulons_mean_ranks_and_cors %>%
  filter(final_rank <=10 | final_rank >= max(final_rank)-9) %>% 
  pull(Regulon)
  
# Plot ranks vs correlations
p = regulons_mean_ranks_and_cors %>%
  ggplot(aes(x=mean_rank,y=mean_cor, color=Regulon %in% ap1_complex, label=ifelse(Regulon %in% top_10_regulons, Regulon, ""))) +
  geom_point(size=0.5) +
  ggrepel::geom_text_repel(max.overlaps = 100, show.legend = FALSE, size=5/2.83465) +
  theme_pubr(base_size = 5) +
  scale_color_manual(values=c("black", "red"), name = "AP-1 member", labels = c("No", "Yes")) +
  theme(axis.line = element_line(linewidth=0.4),
        legend.key.height = unit(0.25, "cm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5)) +
  xlab("mean rank") +
  ylab("mean correlation")
p
save_plot("Suppl.Fig.6a.pdf", p, base_height = 3, base_width = 3)

# Plot average Z-score heatmaps
p = auc_long_df %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  mutate(Cohort = factor(Cohort, levels=c("Young", "Aged"))) %>%
  group_by(Dataset, Cohort, Regulon) %>%
  summarize(m = mean(Regulon_AUC_Z_score)) %>%
  mutate(Regulon = gsub("\\.\\.\\.", "", Regulon)) %>%
  filter(Regulon %in% top_10_regulons) %>%
  mutate(Regulon = factor(Regulon, levels=top_10_regulons)) %>%
  ggplot(aes(x=Dataset,group=Cohort,y=Regulon,fill=m)) +
  geom_tile() +
  facet_grid(.~Cohort) +
  scale_fill_gradient2(low = "steelblue2", high="firebrick2", name="Average Z-score") +
  scale_y_discrete(limits=rev, position = "right", expand = c(0,0), name = element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  theme_pubr(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5),
        axis.line = element_line(linewidth=0.4),
        axis.ticks = element_line(color="black"),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.7, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(2.5, "mm"),
        legend.key.width = unit(3.5, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5)) +
  xlab("Transcription factor")
p
save_plot("Suppl.Fig.6b.pdf", p, base_height = 3, base_width = 3)

# Plot average Z-score heatmaps (Fig. 1h)
p = auc_long_df %>%
  filter(Cohort %in% c("Young", "Aged")) %>%
  mutate(Cohort = factor(Cohort, levels=c("Young", "Aged"))) %>%
  group_by(Dataset, Cohort, Regulon) %>%
  summarize(m = mean(Regulon_AUC_Z_score)) %>%
  mutate(Regulon = gsub("\\.\\.\\.", "", Regulon)) %>%
  filter(Regulon %in% top_regulons) %>%
  mutate(Regulon = factor(Regulon, levels=top_regulons)) %>%
  ggplot(aes(y=Dataset,group=Cohort,x=Regulon,fill=m)) +
  geom_tile() +
  facet_grid(Cohort ~., switch = "y") +
  scale_fill_gradient2(low = "steelblue2", high="firebrick2", name="Average Z-score") +
  scale_y_discrete(limits=rev, position = "right", expand = c(0,0), name = element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  theme_pubr(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5),
        axis.line = element_line(linewidth=0.4),
        axis.ticks = element_line(color="black"),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.4, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(2.5, "mm"),
        legend.key.width = unit(3.5, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=5)) +
  xlab("Transcription factor")
p
save_plot("Fig.1h.pdf", p, base_height = 110, base_width = 35, units="mm")

#### Identify and visualize top regulons for each meta-program ####

# Plot separate rank vs correlation plots for each Dataset-metaGEP combination
p = auc_long_df %>%
  pivot_longer(cols = starts_with("meta")) %>%
  dplyr::rename("Metaprogram"="name",
                "Meta_score"="value") %>%
  group_by(Dataset, Regulon, Metaprogram) %>%
  mutate(Regulon = gsub("\\..*", "", Regulon)) %>% 
  summarize(cor_value = cor(Regulon_AUC_Z_score, Meta_score)) %>%
  group_by(Metaprogram, Dataset) %>%
  arrange(desc(cor_value), .by_group = T) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x=as.numeric(rank),y=cor_value, color=Regulon %in% ap1_complex, label=ifelse(rank<=10, Regulon, ""))) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 100) +
  facet_grid(Metaprogram ~ Dataset) +
  theme_pubr() +
  scale_color_manual(values=c("black", "red"))
p

# Aggregate the plot above into mean ranks and correlations per metaGEP
regulons_mean_ranks_and_cors_per_program = auc_long_df %>%
  pivot_longer(cols = starts_with("meta")) %>%
  dplyr::rename("Metaprogram"="name",
                "Meta_score"="value") %>%
  group_by(Dataset, Regulon, Metaprogram) %>%
  mutate(Regulon = gsub("\\..*", "", Regulon)) %>% 
  summarize(cor_value = cor(Regulon_AUC_Z_score, Meta_score)) %>%
  group_by(Metaprogram, Dataset) %>%
  arrange(desc(cor_value), .by_group = T) %>%
  mutate(rank = row_number()) %>%
  group_by(Regulon, Metaprogram) %>%
  summarize(mean_rank = mean(rank),
            mean_cor = mean(cor_value)) %>% 
  group_by(Metaprogram) %>%
  arrange(mean_rank, .by_group = T) %>%
  mutate(new_rank = row_number()) 

p = regulons_mean_ranks_and_cors_per_program %>%
  ggplot(aes(x=mean_rank,y=mean_cor, color=Regulon %in% ap1_complex, label=ifelse(new_rank<=10, Regulon, ""))) +
  geom_point(size=0.2) +
  ggrepel::geom_text_repel(max.overlaps = 100, show.legend = FALSE, size = 5/2.83465, segment.size=0.2) +
  facet_wrap(Metaprogram ~ ., ncol=2) +
  scale_color_manual(values=c("black", "red"), name = "AP-1 member", labels = c("No", "Yes")) + 
  scale_x_continuous(name = "mean rank", limits = c(-15, 110)) +
  scale_y_continuous(name = "mean correlation", limits = c(-0.4,0.7)) +
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
        legend.text = element_text(size=5))
p
save_plot("Suppl.Fig.9.pdf", p, base_width = 80, base_height = 60, units = "mm")

# Define top regulons per metaGEP
top_regulons_per_program = regulons_mean_ranks_and_cors_per_program %>%
  group_by(Metaprogram) %>%
  filter(mean_cor >= 0.15 & new_rank <= 10) %>%
  dplyr::select(Metaprogram, Regulon, new_rank)

# Plot heatmap with averaged Z-scores
p = auc_long_df %>%
  mutate(Regulon = gsub("\\..*", "", Regulon)) %>%
  group_by(Dataset, Regulon) %>%
  mutate(Regulon_AUC_Z_score_new = scale(Regulon_AUC_Z_score)[,1]) %>% 
  arrange(desc(Regulon_AUC_Z_score)) %>% 
  left_join(top_regulons_per_program, by = c("Predominant_program" = "Metaprogram", "Regulon" = "Regulon")) %>%
  filter(Regulon %in% top_regulons_per_program$Regulon) %>% 
  group_by(Dataset, Regulon, Predominant_program) %>%
  summarize(mean_Z = mean(Regulon_AUC_Z_score)) %>% 
  mutate(Regulon = factor(Regulon, levels = unique(top_regulons_per_program$Regulon))) %>% left_join(top_regulons_per_program, by = c("Regulon" = "Regulon")) %>%
  ggplot(aes(x=Dataset,y=Regulon,fill=mean_Z)) +
  geom_tile() +
  facet_grid(Metaprogram~Predominant_program, scales = "free", space="free", switch = "y") +
  scale_fill_gradient2(low = "steelblue", high = "firebrick") +
  scale_y_discrete(expand = c(0,0), name = "Top TFs per program", position = "right") +
  scale_x_discrete(expand = c(0,0), name = element_blank()) +
  theme_pubr(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=1),
        axis.line = element_line(linewidth=0.4),
        axis.ticks = element_line(color="black"),
        strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.8, linetype="solid"),
        strip.text = element_text(size=5),
        legend.key.height = unit(2.5, "mm"),
        legend.key.width = unit(1.5, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "mm"),
        legend.box = "vertical",
        legend.position = "right",
        legend.text = element_text(size=5)) +
  ggtitle("Predominant program")
p
save_plot("Fig.2f.pdf", p, base_width = 80, base_height = 60, units = "mm")
