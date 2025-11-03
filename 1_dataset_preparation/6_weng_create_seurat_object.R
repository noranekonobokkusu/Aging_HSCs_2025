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

# Read Seurat object
seu = readRDS(paste0(base_dir, "input_data/Weng_data/Weng2024_MergeAllCellTypes_Seurat_240912.rds"))

# Filter low-quality cells and doublets
seu_split = SplitObject(seu, split.by = "Sample")

# Generate empty list 
seu.ls = vector(mode = "list", length = length(seu_split))
ncells_initial_total = 0

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
  scrublet_obj = get_init_scrublet(seurat_obj = seu, sim_doublet_ratio=6, python_home="~/.venvs/MyEnv/bin/python")
  seu[["doublet_scores"]] = scrublet_obj$doublet_scores_obs_
  
  # HSC-enriched samples are highly homogeneous, so 0.2 doesn't work as a threshold
  if (grepl("HSC", names(seu_split)[n])) {
    thr = 0.5
  } else {
    thr = doublet_threshold
  }
  cat(sprintf("Removing %d/%d=%.2f doublets", sum(scrublet_obj$doublet_scores_obs_>=thr), ncells_filtered, (sum(scrublet_obj$doublet_scores_obs_>=thr) / ncells_filtered)), "\n\n")
  seu = subset(seu, doublet_scores < thr)
  seu.ls[[n]] = seu
  
}
seu = merge(x = seu.ls[[1]], y = seu.ls[-1])
seu = JoinLayers(seu)

# Format metadata and subset HSPCs
meta = seu@meta.data
meta = meta %>%
  mutate(Sample = gsub("_.*", "", (gsub("\\..*", "", Sample))),
         Dataset = "Weng2024",
         Cohort = ifelse(grepl("Young", Sample), "Young", "Aged"),
         CellType = STD.CellType,
         Sex = ifelse(Sample %in% c("Aged1", "Aged2", "Aged3", "Aged5", "Young2", "Young4"), "M", "F"),
         Age = "undefined") %>% # Age can be figured out from the paper, but we won't use it anyway
  dplyr::select(-STD.CellType, -ClonalGroup) %>%
  tibble::rownames_to_column("Cell")

rownames(meta) = meta$Cell
seu@meta.data = meta

seu = subset(seu, CellType %in% c("HSC", "MPP", "CMP", "MKP", "MEP", "CLP", "ProB", "GMP", "EryP", "MDP", "LMPP"))

saveRDS(seu, "Weng_seu_formatted.rds")
