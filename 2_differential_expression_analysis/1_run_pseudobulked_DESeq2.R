suppressMessages(suppressWarnings({
  library(Seurat)
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(data.table)
  library(cowplot)
  library(ggpubr)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "2_differential_expression_analysis/"))

# We are interested in deregulation of AP-1 complex members
ap1_complex = c("FOS", "FOSB", "JUN", "JUNB", "JUND", "FOSL1", "FOSL2", "ATF2", "ATF3", "ATF4", "ATF6", "ATF7", "BATF", "BATF2", "BATF3", "MAFA", "MAFB", "MAF", "MAFG", "MAFF", "MAFK")

#### Load the object and pseudobulk it ####
seu_unfiltered = readRDS(paste0(base_dir, "1_dataset_preparation/seu_hsc.total.rds"))

# Mark sample names with their datasets to avoid possible duplicates
seu_unfiltered$Sample = paste0(gsub("_", "-", seu_unfiltered$Sample), "-", seu_unfiltered$Dataset)

# I require a sample to include at least 20 cells to be pseudobulked; I only use datasets that have at least two such samples in both Young and Aged Cohorts
samples_to_use = seu_unfiltered@meta.data %>%
  group_by(Dataset, Sample, Cohort) %>%
  summarize(n=n()) %>%
  filter(n>=20) %>%
  group_by(Dataset) %>%
  mutate(n_young = sum(Cohort=="Young"),
            n_aged = sum(Cohort=="Aged")) %>%
  filter(n_young>=2 & n_aged>=2) %>%
  pull(Sample)

seu = subset(seu_unfiltered, Sample %in% samples_to_use)

# Aggregate counts and metadata per sample
sample_metadata = seu@meta.data %>%
  dplyr::select(Dataset, Cohort, Sex, Sample) %>%
  arrange(Dataset) %>%
  unique()
rownames(sample_metadata) = sample_metadata$Sample

aggregate_expression_result = AggregateExpression(seu, assays = "RNA", return.seurat = F, group.by = c("Sample"))
pseudobulked_counts = aggregate_expression_result$RNA
ribo_genes = grep("^RPS|^RPL", rownames(seu), value=T)
genes_to_keep = setdiff(rownames(seu), ribo_genes)
pseudobulked_counts = pseudobulked_counts[genes_to_keep,rownames(sample_metadata)]

# Check the order of samples is the same in counts and metadata
all(rownames(sample_metadata) == colnames(pseudobulked_counts))

#### Run DESeq2 using Dataset and Sex as covariates ####

# Create DESeq2 object        
ddsColl = DESeqDataSetFromMatrix(pseudobulked_counts, 
                                 colData = sample_metadata, 
                                 design =  ~ Dataset + Sex + Cohort) 

ddsColl$Cohort = relevel(ddsColl$Cohort, ref = "Young")

# Transform counts for data visualization
rld = rlog(ddsColl, blind=TRUE)

# Run DESeq2 differential expression analysis
ddsColl = DESeq(ddsColl, quiet=T)

# Plot dispersion estimates
plotDispEsts(ddsColl)

# Check the coefficients for the comparison
resultsNames(ddsColl)

# Generate results object
res = results(ddsColl, 
              name = "Cohort_Aged_vs_Young",
              alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res = lfcShrink(ddsColl, 
                coef = "Cohort_Aged_vs_Young",
                res=res,
                type = "apeglm",
                quiet = T)

# Turn the DESeq2 results object into a tibble
res_tbl = res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results for AP-1 members
res_tbl %>%
  filter(gene %in% ap1_complex)

#### Make scatterplot with normalized counts per sample for AP-1 members, top-10 up- and down-regulated genes ####

make_scatterplot = function(ddsColl, genes_of_interest, name) {
  
  # Extract normalized counts from dds object
  normalized_counts = counts(ddsColl, normalized = TRUE)
  
  # Subset AP-1 complex
  genes_of_interest_counts = normalized_counts[rownames(normalized_counts) %in% genes_of_interest, ]
  
  # Convert wide matrix to long data frame for ggplot2
  genes_of_interest_counts = data.frame(genes_of_interest_counts)
  genes_of_interest_counts$gene = rownames(genes_of_interest_counts)
  
  # Replace dots in colnames that were introduced by data.frame()
  colnames(genes_of_interest_counts) = gsub("\\.", "-", colnames(genes_of_interest_counts))
  genes_of_interest_counts = melt(setDT(genes_of_interest_counts), 
                    id.vars = c("gene"),
                    variable.name = "Sample") %>% 
    data.frame()
  
  # Join counts data frame with sample_metadata
  genes_of_interest_counts = plyr::join(genes_of_interest_counts, 
                                        as.data.frame(colData(ddsColl)),
                                        by = "Sample")
  
  # Keep only Young and Aged and revelel cohorts
  genes_of_interest_counts = genes_of_interest_counts %>%
    filter(Cohort %in% c("Young", "Aged"))
  genes_of_interest_counts$Cohort = factor(genes_of_interest_counts$Cohort, levels=c("Young", "Aged"))
  
  scatter_p = genes_of_interest_counts %>%
    ggplot(aes(y = value, x = Dataset, fill = Cohort)) +
    geom_boxplot(outlier.shape = NA, linewidth=0.3) +
    geom_point(position=position_jitterdodge(), size=0.3) +
    scale_y_continuous(trans = 'log10') +
    ylab("log10 of normalized expression level") +
    xlab("condition") +
    ggtitle(name) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ gene, scales = "free") +
    theme_pubr(base_size = 5) +
    scale_fill_manual(values = c("#c7d1ca", "#617a6a"),
                      guide = guide_legend(ncol = 2,
                                           title.position = "top",
                      )) +
    scale_y_continuous(name = "Normalized expression", trans = "log10") +
    scale_x_discrete(name = "Dataset") +
    theme(axis.line = element_line(linewidth=0.4),
          axis.text = element_text(angle=90, vjust=0.5),
          strip.background = element_rect(colour="black", fill="white", 
                                          linewidth=0.7, linetype="solid"),
          strip.text = element_text(size=5),
          legend.key.height = unit(0.25, "cm"),
          legend.spacing.y = unit(0, "mm"),
          legend.box.spacing = unit(0, "mm"),
          legend.box = "horizontal",
          legend.margin = margin(0,0,0,0),
          legend.text = element_text(size=5))
  
  scatter_p
  
  save_plot(paste0(gsub(" ", "_", name),  ".normalized_deseq2_expression.pdf"), scatter_p, base_height = 6, base_width = 6)
  
}

top_10_aged = res_tbl %>%
  arrange(padj) %>%
  filter(log2FoldChange > 1) %>%
  slice_head(n=10) %>%
  pull(gene)

top_10_young = res_tbl %>%
  arrange(padj) %>%
  filter(log2FoldChange < -1) %>%
  slice_head(n=10) %>%
  pull(gene)

make_scatterplot(ddsColl, ap1_complex, "AP1 genes")
make_scatterplot(ddsColl, top_10_aged, "Top 10 genes in Aged") # Supplementary Figure 2
make_scatterplot(ddsColl, top_10_young, "Top 10 genes in Young")

#### Save DE results ####

res_table_thres = res_tbl[!is.na(res_tbl$padj), ]
write.table(res_table_thres, "DESeq2_results.Aged_vs_Young.csv", quote = F, sep = ",", col.names = T, row.names = F)

#### Plot Volcano plot ####

# Set thresholds and mark significant hits
padj_cutoff = 0.05
log2fc_cutoff = 1 # (100% increase)

genes_to_label = res_table_thres %>%
  mutate(signif_hit = ifelse(padj<padj_cutoff & abs(log2FoldChange)>log2fc_cutoff, 1, 0),
         sign = ifelse(log2FoldChange<0, "-1", "+1"),
         is_ap1 = ifelse(gene %in% ap1_complex, 1, 0)) %>%
  filter(signif_hit == 1) %>%
  group_by(sign) %>%
  arrange(padj, .by_group = T) %>%
  mutate(rank = ifelse(abs(log2FoldChange) > 1, row_number(), 1000)) %>%  
  mutate(top10 = ifelse(rank <= 10 | is_ap1==1, 1, 0)) %>%
  filter(top10==1) %>%
  pull(gene)

p = res_table_thres %>%
  mutate(signif_hit = ifelse(padj<padj_cutoff & abs(log2FoldChange)>log2fc_cutoff, 1, 0),
         is_ap1 = ifelse(gene %in% ap1_complex, 1, 0)) %>%
  arrange(signif_hit) %>%
  ggplot() +
  geom_vline(xintercept = 1, colour="grey", linetype="dashed", linewidth=0.2) +
  geom_vline(xintercept = -1, colour="grey", linetype="dashed", linewidth=0.2) +
  geom_hline(yintercept = -log10(0.05), colour="grey", linetype="dashed", linewidth=0.2) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = (signif_hit==1)), size=0.5) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), shape = factor(is_ap1)), color = "red", size=0.85) +
  ggrepel::geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(gene %in% genes_to_label, gene, "")), force=5, max.overlaps = 55, direction="both", box.padding = 0.4, min.segment.length = 0, segment.size=0.1, point.padding = 0.1, size = 5/2.845) + 
  scale_shape_manual(
    name = "AP-1 complex",   # Renames the legend title
    values = c("0" = NA, "1" = 16),  # Assigns a valid shape to "1" and hides "0"
    breaks = c("1")  # Ensures only "1" appears in the legend
  ) +
  theme_pubr(base_size = 5) +
  scale_color_manual(values = c("grey", "steelblue"), name = "Significant hit") +
  scale_y_continuous(expand = c(0.01,0.01), name =  expression(-log[10] ~ P[adj])) +
  scale_x_continuous(name =  expression(log[2] ~ "Fold Change")) +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_line(linewidth=0.3),
        legend.position = "none")
p

cowplot::save_plot("Fig.1b.pdf", p, base_height = 2, base_width = 2)
