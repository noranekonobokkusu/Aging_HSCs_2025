suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ggpubr)
  library(fgsea)
  library(cowplot)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "revision"))

# Create a directory if it doesn't exist yet
if (!dir.exists("compare_to_mouse_aging_signature")) dir.create("compare_to_mouse_aging_signature")

# Read mouse-to-human correspondence
mouse_human_orthologs = read_tsv("HOM_MouseHumanSequence.rpt.txt")

# Get the human-mouse gene correspondence data frame
orthologs = mouse_human_orthologs %>%
  dplyr::select(Number = `DB Class Key`,
                Name = `Common Organism Name`,
                Symbol) %>%
  group_by(Number, Name) %>%
  summarize(Symbol = paste(Symbol, collapse = "|")) %>% 
  pivot_wider(names_from = Name, values_from = Symbol) %>%
  dplyr::rename("mouse"=`mouse, laboratory`) %>% 
  separate_rows(human, sep = "\\|") 

# Read DE results
de_results = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep=",", header = T)

# Define aging signature
aging_sig = de_results %>%
  filter(padj<0.05, log2FoldChange>0.5) %>%
  pull(gene)

# Read mouse aging signature from PMID: 33876187 (Consistency>=4 is a threshold used in the original study)
mouse_aging_signature = read_csv(paste0(base_dir, "input_data/signatures/Svendsen_mouse_aging_signature.csv")) %>% filter(Avg>0, Consistency>=4) %>%
  dplyr::select(Genes)

mouse_aging_signature %>%
  left_join(orthologs, by = c("Genes"="mouse")) %>% 
  mutate(is_shared = human %in% aging_sig) %>% 
  filter(is_shared == TRUE) # 7 mouse genes overlap with human signature

contingency_table_counts = de_results %>%
  left_join(orthologs, by = c("gene"="human")) %>% 
  filter(! is.na(mouse)) %>%
  mutate(human_aging = padj<0.05 & log2FoldChange>0.5,
         mouse_aging = mouse %in% mouse_aging_signature$Genes) %>%
  group_by(human_aging, mouse_aging) %>%
  summarize(n=n()) %>%
  pull(n)

fisher.test(matrix(contingency_table_counts, nrow=2)) # Mouse orthologs are enriched among the genes belonging to the human aging signature; p-value=0.00029

# We can also test it with a binomial test to get an almost identical p-value of 0.00024
p=119/(119+11409)
binom.test(7, 118, p, alternative = "two.sided")

# Gene set over-representation; mouse gene sets were downloaded from https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp
mouse_hallmark_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/mh.all.v2025.1.Mm.symbols.gmt.txt"))
mouse_c2_pathways = gmtPathways(paste0(base_dir, "input_data/signatures/m2.all.v2025.1.Mm.symbols.gmt.txt"))

pathways = c(mouse_hallmark_pathways, mouse_c2_pathways)

# Define the universe of all mouse genes 
background = orthologs %>%
  filter(! is.na(mouse)) %>%
  pull(mouse) %>%
  unique()

# Compute over-representation
foraRes = fora(pathways = pathways, 
                genes    = mouse_aging_signature$Genes,
                universe = background, 
                minSize  = 5)
data.table::fwrite(foraRes, file=paste0("compare_to_mouse_aging_signature/Supplementary_Table_4.tsv"), sep="\t", sep2=c("", " ", ""), append = F)

p_a = foraRes %>%
  filter(grepl("HALLMARK|REACTOME|KEGG", pathway)) %>%
  slice_head(n=20) %>%
  arrange(desc(padj)) %>%
  mutate(pathway=factor(pathway, levels = .data$pathway)) %>% 
  ggplot(aes(y=pathway, x=foldEnrichment, fill=padj)) +
  geom_bar(stat = "identity") +
  theme_pubr() +
  scale_fill_gradient(low="steelblue3", high="powderblue", breaks=scales::pretty_breaks(n = 2), limits=c(0, 0.02))

pathways = foraRes %>%
  filter(grepl("HALLMARK|REACTOME|KEGG", pathway)) %>%
  slice_head(n=20) %>%
  arrange(desc(padj)) %>%
  pull(pathway)

human_gsea_results = data.table::fread(file=paste0(base_dir, "2_differential_expression_analysis/DESeq2_GSEA_results.txt"), sep="\t", sep2=c("", " ", ""))

p_b = human_gsea_results %>%
  filter(pathway %in% pathways) %>%
  mutate(pathway = factor(pathway, levels = pathways)) %>%
  ggplot(aes(x=NES,y=pathway,fill=padj, label=round(padj, 3))) +
  geom_bar(stat="identity") +
  geom_text(aes(x=NES-0.3*sign(NES))) +
  scale_fill_gradientn(colours = c("steelblue2", "grey", "grey"), values = c(0, 0.05, 1)) +
  theme_pubr()

p = ggarrange(p_a, p_b, labels = c("a", "b"), ncol = 1)
save_plot("compare_to_mouse_aging_signature/Suppl.Fig.14.png", p, base_height = 12, base_width = 9, dpi=600)
