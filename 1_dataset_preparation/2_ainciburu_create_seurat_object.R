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


# Read count files as a Seurat object and filter low-quality cells and doublets
file_list = list.files(path = paste0(base_dir, "input_data/Ainciburu_data/GSE180298_RAW"), pattern = ".*(young|elderly).*\\.h5", full.names = T)

seurat_objects = list()
ncells_initial_total = 0

for (n in 1:length(file_list)) {
  h5_file = file_list[[n]]
  sample_name = gsub(".*_", "", gsub("_filtered.*", "", h5_file))
  cat(sprintf("Processing sample %s, %d/%d", sample_name, n, length(file_list)), "\n")
  h5_data = Read10X_h5(h5_file)
  seu = CreateSeuratObject(counts = h5_data, min.cells = 5, min.features = 100)
  seu = RenameCells(seu, new.names = paste0(gsub("-.*", "", colnames(seu)), "-", sample_name))
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
  cat(sprintf("Removing %d/%d=%.2f doublets", sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold), ncells_filtered, (sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold) / ncells_filtered)), "\n\n")
  seu = subset(seu, doublet_scores < doublet_threshold)
  
  seurat_objects[[sample_name]] = seu
  seurat_objects[[sample_name]]$orig.ident = sample_name
}

seu = merge(seurat_objects[[1]], seurat_objects[-1])
seu = JoinLayers(seu)

# Read and add metadata
young_meta = read.table(paste0(base_dir, "input_data/Ainciburu_data/GSE180298_RAW/metadata_young.txt"), sep = "\t", header = T)
old_meta = read.table(paste0(base_dir, "input_data/Ainciburu_data/GSE180298_RAW/metadata_elderly.txt"), sep = "\t", header = T)
combined_meta = rbind(young_meta, old_meta) %>%
  rownames_to_column("Cell") %>%
  mutate(Cell = gsub("_", "-", Cell)) %>%
  dplyr::select(Cell, CellType)

meta = seu@meta.data
meta$Cell = rownames(meta)
meta = meta %>%
  left_join(combined_meta, by = "Cell") %>%
  mutate(Sample = orig.ident,
         Dataset = "Ainciburu2023",
         Cohort = ifelse(grepl("young", Sample), "Young", "Aged"),
         Age = case_when(orig.ident == "young1" ~ 20,
                         orig.ident == "young2" ~ 23,
                         orig.ident == "young3" ~ 20,
                         orig.ident == "young4" ~ 19,
                         orig.ident == "young5" ~ 23,
                         orig.ident == "elderly1" ~ 61,
                         orig.ident == "elderly2" ~ 74,
                         orig.ident == "elderly3" ~ 72),
         Sex = ifelse(grepl("young", orig.ident), "F", "M"))

rownames(meta) = meta$Cell
seu@meta.data = meta

# Some of the CellTypes are NAs as those cells were excluded by the author filters
saveRDS(seu, "Ainciburu_seu_formatted.rds")
