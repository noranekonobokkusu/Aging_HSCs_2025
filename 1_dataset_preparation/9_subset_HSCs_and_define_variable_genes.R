suppressMessages(suppressWarnings({
  library(tidyverse)
  library(Seurat)
  library(cowplot)
  library(MuDataSeurat)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "1_dataset_preparation/"))

#### Load the entire annotated HSPCs object and subset HSCs ####
seu_total = readRDS("seu_full_object_hscs_marked.rds")
seu_hsc = subset(seu_total, is_hsc=="HSC")
cat(sprintf("The final dataset is comprised of %d HSCs", ncol(seu_hsc)), "\n")

# Make sure there are no NAs in metadata
seu_hsc@meta.data %>% filter(if_any(everything(), is.na))

#### Save HSCs as an .rds object ####
saveRDS(seu_hsc, file = "seu_hsc.total.rds")

#### Identify variable genes for scVI integration and cNMF, subset seu_hsc, and save it as an .h5ad object ####

# We expect HSC variation can be described by a relatively small number of genes; we settled on selecting genes variable (within top-3000) in at least 4/7 datasets; other approaches for variable genes selection also recover aging-associated metaGEP that becomes more active with age

# Read data
seu = seu_hsc

# Define genes with 0 counts in at least one dataset to exclude them later
zero_count_genes = c()
for (d in unique(seu$Dataset)) {
  seu_temp = subset(seu, Dataset==d)
  zero_count_genes = c(zero_count_genes, names(rowMeans(seu_temp@assays$RNA$counts)[rowMeans(seu_temp@assays$RNA$counts)==0]))
}
zero_count_genes = unique(zero_count_genes)

# Remove ribosomal protein genes and genes with zero counts
ribo_genes = grep("^RPS|^RPL", rownames(seu), value=T)
genes_to_keep = setdiff(rownames(seu), union(zero_count_genes, ribo_genes))

# Define variable genes within each dataset
seu$Sample = paste0(seu$Sample, "_", seu$Dataset)

variable_genes_list = list()
datasets = unique(seu$Dataset)

for (i in 1:length(datasets)) {
  seu_d = subset(seu, Dataset == datasets[i])
  counts = GetAssayData(seu_d, layer="counts", assay="RNA")
  counts = counts[rownames(counts) %in% genes_to_keep,]
  seu_d = CreateSeuratObject(counts=counts, meta.data = seu_d@meta.data)
  DefaultAssay(seu_d) = "RNA"
  seu_d = NormalizeData(seu_d)
  seu_d = FindVariableFeatures(seu_d, nfeatures=3000)
  variable_genes_list[[i]] = VariableFeatures(seu_d)
}

# Filter genes that occur in at least 4/7 datasets
gene_occurrences = table(unlist(variable_genes_list))
common_genes = names(gene_occurrences[gene_occurrences >= 4])
cat(sprintf("There are %d variable genes across 4/%d datasets", length(common_genes), length(datasets)), "\n")

# Make a new Seurat object with variable genes
counts = as.matrix(GetAssayData(seu, layer = "counts"))
counts = counts[rownames(counts) %in% common_genes,]
seu_var = CreateSeuratObject(counts=counts, meta.data = seu@meta.data)

# Save the entire object as an H5AD file
WriteH5AD(seu_var, paste0("seu_hsc.var_genes.h5ad"), assay="RNA")

#### Save HSCs as an .h5ad object for SCENIC ####
counts = as.matrix(GetAssayData(seu_hsc, layer = "counts"))
seu_hsc = CreateSeuratObject(counts=counts, meta.data = seu_hsc@meta.data)
MuDataSeurat::WriteH5AD(seu_hsc, "seu_hsc.total.h5ad", assay="RNA")
