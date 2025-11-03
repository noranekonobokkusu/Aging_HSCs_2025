Sys.setenv(RETICULATE_PYTHON = "~/.venvs/MyEnv/bin/python")

suppressMessages({
  library(reticulate)
  library(tidyverse)
  library(Seurat)
  library(scrubletR)
})

use_virtualenv("~/.venvs/MyEnv", required = TRUE)

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "1_dataset_preparation/"))

# Define thresholds
mito.high = 0.2
count.low = 100
feature.low = 500
upper.quantile = 0.99
doublet_threshold = 0.2

# Define reporting function
summarize_metric = function(data, metric_name) {
  median_val = median(data[[metric_name]], na.rm = TRUE)
  q5 = quantile(data[[metric_name]], 0.05, na.rm = TRUE)
  q95 = quantile(data[[metric_name]], 0.95, na.rm = TRUE)
  
  # Format the message
  message = sprintf(
    "%s: Median = %.2f, 5th = %.2f, 95th = %.2f",
    metric_name, median_val, q5, q95
  )
  
  return(message)
}

# Read expression
folders = list.dirs(paste0(base_dir, "input_data/Safina_data"), full.names = T, recursive = F)[-1] # -1 removes the base folder
data.tib = tibble(Sample = basename(folders), Folder = folders)

# Generate empty list 
seu.ls = vector(mode = "list", length = nrow(data.tib))
ncells_initial_total = 0

# Load Seurat objects and filter low-quality cells and doublets
for (n in 1:nrow(data.tib)) {
  cat(sprintf("Processing sample %s, %d/%d", data.tib$Sample[n], n, nrow(data.tib)), "\n")
  # Load the counts matrix and create Seurat object
  counts.mat = Read10X(data.dir = paste0(data.tib$Folder[n], "/filtered_feature_bc_matrix"))
  seu = CreateSeuratObject(counts = counts.mat, project = data.tib$Sample[n], min.cells = 5, min.features = 100)
  seu = PercentageFeatureSet(seu, "^MT-", col.name = "percent.mito")
  seu$percent.mito = seu$percent.mito/100
  metrics = c("percent.mito", "nCount_RNA", "nFeature_RNA")
  messages = sapply(metrics, function(metric) summarize_metric(seu@meta.data, metric))
  cat(paste(messages, collapse = "\n"), "\n")
  
  # Filter based on metrics
  ncells_initial = ncol(seu)
  ncells_initial_total = ncells_initial_total + ncells_initial
  seu = subset(seu, percent.mito < mito.high & 
                 nCount_RNA > count.low &
                 nCount_RNA < quantile(seu$nCount_RNA, upper.quantile) &
                 nFeature_RNA > feature.low & 
                 nFeature_RNA < quantile(seu$nFeature_RNA, upper.quantile))
  ncells_filtered = ncol(seu)
  cat(sprintf("Removing %d/%d=%.2f low-quality cells", ncells_initial-ncells_filtered, ncells_initial, (ncells_initial-ncells_filtered)/ncells_initial), "\n")
  
  # Remove doublets
  scrublet_obj = get_init_scrublet(seurat_obj = seu, sim_doublet_ratio=6, n_prin_comps=min(ncol(seu)-1, 50), python_home="~/.venvs/MyEnv/bin/python")
  seu[["doublet_scores"]] = scrublet_obj$doublet_scores_obs_
  cat(sprintf("Removing %d/%d=%.2f doublets", sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold), ncells_filtered, (sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold) / ncells_filtered)), "\n\n")
  seu = subset(seu, doublet_scores < doublet_threshold)
  seu.ls[[n]] = seu
  
}

# Merge a list of objects
seu = merge(x = seu.ls[[1]], y = seu.ls[-1], add.cell.ids = data.tib$Sample)
seu = JoinLayers(seu)
seu = RenameCells(seu, new.names = gsub("-", ".", gsub("-1", "", colnames(seu))))

# Perform rough integration and select HSPCs
seu[["RNA"]] = split(seu[["RNA"]], f = seu$orig.ident)
seu = NormalizeData(seu)
seu = FindVariableFeatures(seu)
seu = ScaleData(seu)
seu = RunPCA(seu)
seu = IntegrateLayers(object = seu, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "rpca.integrated", verbose = FALSE)
seu[["RNA"]] = JoinLayers(seu[["RNA"]])
seu = FindNeighbors(seu, reduction = "rpca.integrated", dims = 1:30)
seu = FindClusters(seu, resolution = 1.5)
seu = RunUMAP(seu, dims = 1:30, reduction = "rpca.integrated")
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)
# Subset HSPCs
seu = subset(seu, seurat_clusters %in% c(33,30,36,28))
seu@meta.data[c("RNA_snn_res.1", "seurat_clusters")] = NULL

# Add metadata
metadata = read.table(paste0(base_dir, "input_data/Safina_data/metadata.txt"), header = T)
meta = seu@meta.data
meta = meta %>%
  rownames_to_column("Cell") %>%
  mutate(Sample = gsub("-[0-9]*", "", orig.ident),
         Dataset = "Safina2024") %>%
  left_join(metadata, by="Sample")
rownames(meta) = meta$Cell
seu@meta.data = meta
Idents(seu) = "Sample"

# Drop unnecessary information
counts = GetAssayData(seu, layer = "counts")
new_seu = CreateSeuratObject(counts = counts, meta.data = seu@meta.data)

# Save RDS
saveRDS(new_seu, "Safina_seu_formatted.rds")
