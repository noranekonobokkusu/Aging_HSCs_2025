suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(MuDataSeurat)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

# Make a folder for random signatures
if (!dir.exists("random_gene_sets_control")) dir.create("random_gene_sets_control")

# Read Seurat object
seu = ReadH5AD(paste0(base_dir, "3_scVI_integration/integrated_seu.h5ad"))

# Read aging signature
de_results = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep=",", header = T)
aging_sig = de_results %>%
  filter(padj<0.05, log2FoldChange>1) %>%
  pull(gene)

# Define expression bins
seu = NormalizeData(seu)
mean_expression = rowMeans(seu@assays$RNA@data)

num_bins = 50
bins = cut(mean_expression, 
           breaks = quantile(mean_expression, 
                             probs = seq(0, 1, length.out = num_bins + 1)), 
           include.lowest = TRUE, 
           labels = FALSE)

gene_bins = data.frame(
  Gene = rownames(seu@assays$RNA@data),
  Expression = mean_expression,
  Bin = bins
)

aging_genes_bins = gene_bins %>%
  filter(Gene %in% aging_sig)

matching_genes = aging_genes_bins %>%
  left_join(gene_bins, by = "Bin", suffix = c("_target", "_match")) %>%
  group_by(Gene_target) %>%
  slice_sample(n = 1) %>%
  select(Gene_target, Bin, MatchedGene = Gene_match)

sampled_genes = matching_genes$MatchedGene

# Visualise real and sampled genes
data = data.frame(Expression = mean_expression)
aging_genes_values = mean_expression[rownames(seu@assays$RNA@data) %in% aging_sig]
sampled_genes_values = mean_expression[rownames(seu@assays$RNA@data) %in% sampled_genes]

p_real = ggplot(data, aes(x = Expression)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = aging_genes_values, color = "red", linetype = "dashed", size = 0.8) +
  labs(
    title = "Genes from my aging signature",
    x = "Averaged Expression",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

p_sampled = ggplot(data, aes(x = Expression)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = sampled_genes_values, color = "red", linetype = "dashed", size = 0.8) +
  labs(
    title = "Randomly sampled genes",
    x = "Averaged Expression",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggarrange(p_real, p_sampled, ncol=1)

# Looks ok, let's use this procedure to generate and save twenty random signatures
ngenes = sum(gene_bins$Gene %in% aging_sig)
nsum = 0
nsignatures = 20
random_gene_sets = data.frame(gene_number = 1:ngenes)
while (nsum < nsignatures) {
  matching_genes = aging_genes_bins %>%
    left_join(gene_bins, by = "Bin", suffix = c("_target", "_match"), relationship = "many-to-many") %>%
    # filter(! Gene_match %in% aging_sig) %>%
    group_by(Gene_target) %>%
    slice_sample(n = 1) %>%
    select(Gene_target, Bin, MatchedGene = Gene_match) %>% pull(MatchedGene)
  if (sum(duplicated(matching_genes))==0) {
    nsum = nsum + 1
    print(length(intersect(aging_sig, matching_genes)))
    random_gene_sets = cbind(random_gene_sets, matching_genes)
  }
}

# Save
random_gene_sets$gene_number = NULL
colnames(random_gene_sets) = paste0("Random_set", 1:nsignatures)
write.table(random_gene_sets, "random_gene_sets_control/random_gene_sets.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)
