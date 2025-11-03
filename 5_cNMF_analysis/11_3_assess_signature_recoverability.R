suppressMessages(suppressWarnings({
  library(RColorBrewer)
  library(fgsea)
  library(tidyverse)
  library(ggpubr)
}))


base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

#### Prepare gene sets ####

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

main_pathways = list(garcia_signature, aging_sig, mep_sig, gmp_sig, clp_sig, ap1_complex)

names(main_pathways) = c("Garcia_quiescence", "Aging_signature", "MEP_signature", "GMP_signature", "CLP_signature", "AP-1_complex")

pathways = c(main_pathways, hallmark_pathways, c2_cp_pathways)

# Add random gene sets
random_gene_sets = read.table(paste0("random_gene_sets_control/random_gene_sets.txt"), header = TRUE, stringsAsFactors = FALSE)
random_gene_sets = lapply(random_gene_sets, function(column) column[!is.na(column)])

all_pathways = c(pathways, random_gene_sets)

pathways_of_interest = all_pathways[c("Aging_signature", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Garcia_quiescence", "REACTOME_CELL_CYCLE", "MEP_signature", "GMP_signature", "CLP_signature", paste0("Random_set", 1:20))]

#### Collect all GSEA results into a single data frame ####

merged_t_hsc = read.table("random_gene_sets_control/total_cNMF_output.csv", sep=",", header = T, row.names = 1)

fgseaRes_total = data.frame()

pb = txtProgressBar(min = 0, max = nrow(merged_t_hsc), style = 3)

for (i in 1:nrow(merged_t_hsc)) {
  
  entry = merged_t_hsc[i,]
  ranks = as.numeric(entry)
  names(ranks) = colnames(entry)
  
  # GSEA
  fgseaRes = fgsea(pathways = pathways_of_interest, 
                   stats    = ranks,
                   minSize  = 5,
                   eps      = 0.0,
                   maxSize  = 500,
                   nPermSimple = 10000)
  
  fgseaRes$Program = rownames(entry)
  fgseaRes_total = rbind(fgseaRes_total, fgseaRes)
  setTxtProgressBar(pb, i)
  
}

close(pb)

data.table::fwrite(fgseaRes_total, file="random_gene_sets_control/gsea_results_with_random_gene_sets.txt", sep="\t", sep2=c("", " ", ""), append = F)


#### Visualize recoverability of gene sets across stable cNMF runs ####

# Let's look at individual gene sets first

pathways_of_interest = c("Aging_signature", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Garcia_quiescence", "REACTOME_CELL_CYCLE", "MEP_signature", "GMP_signature", "CLP_signature", paste0("Random_set", 1:20))

# 1. Across individual samples (sorted by sample size)

samples = selected_ks$Sample
sample_order = readRDS("annotated_seu.rds")@meta.data %>% filter(Cohort %in% c("Young", "Aged")) %>% group_by(Sample, Dataset) %>% summarise(n=n()) %>% arrange(desc(n)) %>% mutate(Sample = paste0(Sample, "_", Dataset)) %>% filter(n>=100) %>% pull(Sample)

p1 = fgseaRes_total %>%
  group_by(Program) %>%
  arrange(pval, .by_group = T) %>%
  mutate(bonferroni_padj = p.adjust(pval, method="bonferroni")) %>% 
  dplyr::select(pathway, bonferroni_padj, Program, NES) %>%
  mutate(Run = gsub(":.*", "", Program),
         Sample = gsub("_K.*", "", Program),
         sign_pos_enrich = ifelse(NES>0 & bonferroni_padj<0.05, 1, 0),
         pathway = factor(pathway, levels=pathways_of_interest)) %>%
  group_by(Run, Sample, pathway) %>%
  summarize(n = ifelse(sum(sign_pos_enrich)>0, 1, 0)) %>%
  group_by(Sample, pathway) %>%
  summarize(frac = sum(n)/n(),
            nruns = n()) %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Sample) %>%
  mutate(Sample2 = paste0(Sample, " n=", nruns)) %>%
  mutate(Sample2 = factor(Sample2)) %>%
  ggplot(aes(x=pathway,y=Sample2,fill=frac)) +
  geom_tile(color="black") +
  scale_x_discrete(name = "Gene set",  expand = c(0,0)) +
  scale_y_discrete(name = "Sample", expand = c(0,0)) +
  scale_fill_distiller(palette = "Oranges", direction = 1, name = "Fraction recovered") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
p1

# 2. Aggregated across datasets

fgseaRes_total %>%
  group_by(Program) %>%
  arrange(pval, .by_group = T) %>%
  mutate(bonferroni_padj = p.adjust(pval, method="bonferroni")) %>% 
  dplyr::select(pathway, bonferroni_padj, Program, NES) %>%
  mutate(Run = gsub(":.*", "", Program),
         Sample = gsub("_K.*", "", Program),
         Dataset = gsub(".*_", "", Sample),
         sign_pos_enrich = ifelse(NES>0 & bonferroni_padj<0.05, 1, 0),
         pathway = factor(pathway, levels=pathways_of_interest)) %>%
  group_by(Run, Sample, pathway) %>%
  group_by(Run, Sample, Dataset, pathway) %>%
  summarize(n = ifelse(sum(sign_pos_enrich)>0, 1, 0)) %>%
  group_by(Dataset, pathway) %>%
  summarize(frac = sum(n)/n(),
            count = n()) %>%
  mutate(Dataset = paste0(Dataset, " (n=", count, ")")) %>%
  ggplot(aes(x=pathway,y=Dataset,fill=frac)) +
  geom_tile() +
  scale_x_discrete(name = "Gene set",  expand = c(0,0)) +
  scale_y_discrete(name = "Dataset", expand = c(0,0)) +
  scale_fill_distiller(palette = "Oranges", direction = 1, name = "Fraction recovered") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

# Let's add two new "gene sets" - Aging sig & TNF-a, and Aging sig & TNF-a & Quiescence

t = fgseaRes_total %>%
  group_by(Program) %>%
  arrange(pval, .by_group = T) %>%
  mutate(bonferroni_padj = p.adjust(pval, method="bonferroni")) %>% 
  dplyr::select(pathway, bonferroni_padj, Program, NES) %>%
  mutate(sign_pos_enrich = ifelse(NES>0 & bonferroni_padj<0.05, 1, 0)) %>%
  dplyr::select(-c(bonferroni_padj, NES))

new_pathways = t %>%
  mutate(Aging_TNFa_signature = ifelse(sign_pos_enrich[pathway == "Aging_signature"] & sign_pos_enrich[pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"], 1, 0),
         Aging_TNFa_Quiescence_signature = ifelse(sign_pos_enrich[pathway == "Aging_signature"] & sign_pos_enrich[pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] & sign_pos_enrich[pathway == "Garcia_quiescence"], 1, 0)) %>%
  dplyr::select(-c(pathway, sign_pos_enrich)) %>%
  unique() %>%
  pivot_longer(cols = -Program, names_to = "pathway", values_to = "sign_pos_enrich") %>%
  relocate(pathway, .before = Program)

pathways_of_interest = c("Aging_TNFa_Quiescence_signature", "Aging_TNFa_signature", pathways_of_interest)

# Aggregated across datasets, with two additional "gene sets"

t %>%
  bind_rows(new_pathways) %>%
  mutate(Run = gsub(":.*", "", Program),
         Sample = gsub("_K.*", "", Program),
         Dataset = gsub(".*_", "", Sample),
         pathway = factor(pathway, levels=pathways_of_interest)) %>%
  group_by(Run, Sample, Dataset, pathway) %>%
  summarize(n = ifelse(sum(sign_pos_enrich)>0, 1, 0)) %>%
  group_by(Dataset, pathway) %>%
  summarize(frac = sum(n)/n(),
            count = n()) %>%
  mutate(Dataset = paste0(Dataset, " (n=", count, ")")) %>%
  ggplot(aes(x=pathway,y=Dataset,fill=frac, label=count * frac)) +
  geom_tile() +
  scale_x_discrete(name = "Gene set (or combination)",  expand = c(0,0)) +
  scale_y_discrete(name = "Dataset", expand = c(0,0)) +
  geom_text() +
  scale_fill_distiller(palette = "Oranges", direction = 1, name = "Fraction recovered") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

# Let's average random gene set results for a more compact panel

pathways_of_interest = c("Aging_TNFa_Quiescence_signature", "Aging_TNFa_signature", "Aging_signature", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Garcia_quiescence", "REACTOME_CELL_CYCLE", "MEP_signature", "GMP_signature", "CLP_signature", "Random_set")

p2 = t %>%
  bind_rows(new_pathways) %>%
  mutate(Run = gsub(":.*", "", Program),
         Sample = gsub("_K.*", "", Program),
         Dataset = gsub(".*_", "", Sample)) %>%
  group_by(Run, Sample, Dataset, pathway) %>%
  summarize(n = ifelse(sum(sign_pos_enrich)>0, 1, 0)) %>%
  mutate(pathway = gsub("Random_set.*", "Random_set", pathway)) %>%
  group_by(Dataset, pathway) %>%
  summarize(frac = sum(n)/n(),
            count = n()) %>%
  mutate(nruns=count[pathway=="Aging_signature"],
         Dataset = paste0(Dataset, " (n=", nruns, ")"),
         pathway=factor(pathway, levels=pathways_of_interest)) %>%
  ggplot(aes(x=pathway,y=Dataset,fill=frac, label=round(frac,2))) +
  geom_tile() +
  scale_x_discrete(name = "Gene set (or combination)",  expand = c(0,0)) +
  scale_y_discrete(name = "Dataset", expand = c(0,0)) +
  geom_text() +
  scale_fill_distiller(palette = "Oranges", direction = 1, name = "Fraction recovered") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
p2

p = ggarrange(p1, p2, ncol=1, labels = c("a", "b"), heights = c(0.7, 0.4))
cowplot::save_plot("random_gene_sets_control/Suppl.Fig.controls.pdf", p, base_height = 15, base_width = 10)


# This is the same panel as above but transposed, which used to be a panel in Figure 2
t %>%
  bind_rows(new_pathways) %>%
  mutate(Run = gsub(":.*", "", Program),
         Sample = gsub("_K.*", "", Program),
         Dataset = gsub(".*_", "", Sample)) %>%
  group_by(Run, Sample, Dataset, pathway) %>%
  summarize(n = ifelse(sum(sign_pos_enrich)>0, 1, 0)) %>%
  mutate(pathway = gsub("Random_set.*", "Random_set", pathway)) %>%
  group_by(Dataset, pathway) %>%
  summarize(frac = sum(n)/n(),
            count = n()) %>%
  mutate(nruns=count[pathway=="Aging_signature"],
         Dataset = paste0(Dataset, " (n=", nruns, ")"),
         pathway=factor(pathway, levels=rev(pathways_of_interest))) %>%
  ggplot(aes(x=Dataset,y=pathway,fill=frac, label=round(frac,2))) +
  geom_tile(color="black", linewidth=0.1) +
  geom_text(size=5/2.83465) +
  scale_fill_distiller(palette = "Oranges", direction = 1, name="Fraction", guide = guide_colorbar(barwidth = 0.3, barheight = 2.5)) +
  scale_x_discrete(expand = c(0,0), position = "top", name="") +
  scale_y_discrete(expand = c(0,0), name="") +
  theme_pubr(base_size = 5) +
  theme(axis.text.x = element_text(angle=45,vjust=0,hjust=0),
        axis.line = element_line(linewidth=0.4),
        legend.position = "right")

