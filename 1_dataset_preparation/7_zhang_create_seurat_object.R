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

# Read count matrix and metadata
m = as.matrix(read.table(file = paste0(base_dir, "input_data/Zhang_data/hspc_raw_matrix.txt")))
meta = read.table(paste0(base_dir, "input_data/Zhang_data/Cell_metadata.txt"), header = T)
extra_meta = read.table(paste0(base_dir, "input_data/Zhang_data/extra_meta.txt"), header = T) # Received from the authors and reformatted
meta = meta %>%
  rename(CellType = RNA_clusters, Cell = cell) %>%
  dplyr::select(-immunophenotypes, -developmental_stages) %>%
  mutate(Sample = gsub("^[A-Z0-9]*_", "", gsub("_L.*_bar.*", "", Cell)),
         Dataset = "Zhang2022") %>%
  left_join(extra_meta, by = "Sample")
meta %>% filter(if_any(everything(), is.na)) # Check for NAs
rownames(meta) = meta$Cell
head(meta)

# Create Seurat Object
seu = CreateSeuratObject(counts = m, meta.data = meta)

# Filter low-quality cells and doublets
seu_split = SplitObject(seu, split.by = "Sample")

# Generate empty list 
seu.ls = vector(mode = "list", length = length(seu_split))
ncells_initial_total = 0

# Load Seurat objects
for (n in 1:length(seu_split)) {
  
  cat(sprintf("Processing sample %s, %d/%d", names(seu_split)[n], n, nrow(data.tib)), "\n")
  # Get the object
  counts = GetAssayData(seu_split[[n]], layer = "counts")
  seu = CreateSeuratObject(counts = counts, min.cells = 5, min.features = 100, meta.data = seu_split[[n]]@meta.data)
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
seu = merge(x = seu.ls[[1]], y = seu.ls[-1])
seu = JoinLayers(seu)

# Save RDS
saveRDS(seu, file = "Zhang_seu_formatted.rds")
