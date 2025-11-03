suppressMessages(suppressWarnings({
  library(fgsea)
  library(dplyr)
  library(Seurat)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "2_differential_expression_analysis/"))

#### Define a list of pathways ####

# Human dormancy signature from https://doi.org/10.1016/j.stem.2021.07.003
garcia_table = readxl::read_xlsx(paste0(base_dir, "input_data/signatures/mmc2.xlsx"), sheet = 1)
garcia_signature = garcia_table %>%
  filter(genename!="None") %>%
  filter(logFC < -2,FDR<0.00001) %>%
  pull(genename) %>%
  unique()

graham_signature = read.table(paste0(base_dir, "input_data/signatures/graham_quiescence_up.txt"), header = F)$V1

s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
ap1_complex = c("FOS", "FOSB", "JUN", "JUNB", "JUND", "FOSL1", "FOSL2", "ATF2", "ATF3", "ATF4", "ATF6", "ATF7", "BATF", "BATF2", "BATF3", "MAFA", "MAFB", "MAF", "MAFG", "MAFF", "MAFK")
de_results = read.table("DESeq2_results.Aged_vs_Young.csv", sep=",", header = T)
aging_signature = de_results %>%
  filter(padj<0.05, log2FoldChange>1) %>%
  pull(gene)

hallmark_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/h.all.v2024.1.Hs.symbols.gmt.txt"))
c2_cp_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/c2.cp.v2024.1.Hs.symbols.gmt.txt"))

pathways = list(garcia_signature, graham_signature, ap1_complex, aging_signature, s.genes, g2m.genes)

names(pathways) = c("Garcia_quiescence", "Graham_quiescence", "AP-1 complex", "Aging signature", "s.genes", "g2m.genes")

pathways = c(pathways, hallmark_pathways, c2_cp_pathways)

#### Read DE results ####
de_results = read.table("DESeq2_results.Aged_vs_Young.csv", header = T, sep = ",")
ranks = as.numeric(de_results$log2FoldChange)
names(ranks) = de_results$gene

#### Run GSEA ####
fgseaRes = fgsea(pathways = pathways, 
                 stats    = ranks,
                 minSize  = 5,
                 eps      = 0.0,
                 maxSize  = 500,
                 nPermSimple = 10000)

fgseaRes = fgseaRes %>% arrange(padj)

data.table::fwrite(fgseaRes, file="DESeq2_GSEA_results.txt", sep="\t", sep2=c("", " ", ""))

fgseaRes = fgseaRes %>% filter(! pathway %in% c("Aging signature", "AP-1 complex"))
# collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
#                                       pathways, ranks)
# collapsed_fgseaRes <- fgseaRes %>% filter(pathway %in% collapsedPathways$mainPathways)

df = fgseaRes %>%
  arrange(padj) %>%
  slice_head(n=15) %>%
  dplyr::select(ID=pathway, padj, NES) %>%
  mutate(padj = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels = ID))

p1 = df %>%
  ggplot(aes(x=NES, y=ID, fill = padj)) +
  geom_bar(stat = "identity") +
  theme_pubr(base_size = 5) +
  scale_fill_gradient(low = "powderblue", high = "steelblue", trans="log", breaks = scales::pretty_breaks(), name =  expression(-log[10] ~ P[adj])) +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_line(linewidth=0.3), legend.direction = "vertical",
        legend.position = "right",  # Moves the legend above the plot
        legend.justification = "left",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.3, "cm"))
p1 
save_plot("Fig.1c.pdf", p1, base_height = 2, base_width = 3.2)
