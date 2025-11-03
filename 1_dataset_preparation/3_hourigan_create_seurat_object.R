# This has changed a lot but was originally based on Peter's script
Sys.setenv(RETICULATE_PYTHON = "~/.venvs/MyEnv/bin/python")

suppressMessages({
  library(reticulate)
  library(tidyverse)
  library(Seurat)
  library(readxl)
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

# Read in the supplemental data
supplemental = read_excel(paste0(base_dir, "input_data/Hourigan_data/20230525_BM_Supplement.xlsx"))

# Read in the cell annotations
annotations = read.csv(paste0(base_dir, "input_data/Hourigan_data/celltype.csv"), row.names = 1)

# Expression and metadata files
folders = list.dirs(paste0(base_dir, "input_data/Hourigan_data"), full.names = T)[-1] # -1 removes the base folder
data.tib = tibble(Sample = basename(folders), Folder = folders)

# Generate empty list 
seu.ls = vector(mode = "list", length = nrow(data.tib))
ncells_initial_total = 0

# Load Seurat objects
for (n in 1:nrow(data.tib)) {
  
  cat(sprintf("Processing sample %s, %d/%d", data.tib$Sample[n], n, nrow(data.tib)), "\n")
  
  # Load the counts matrix and create Seurat object
  counts.mat = Read10X(data.dir = data.tib$Folder[n])
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
  scrublet_obj = get_init_scrublet(seurat_obj = seu, sim_doublet_ratio=6, python_home="~/.venvs/MyEnv/bin/python")
  seu[["doublet_scores"]] = scrublet_obj$doublet_scores_obs_
  cat(sprintf("Removing %d/%d=%.2f doublets", sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold), ncells_filtered, (sum(scrublet_obj$doublet_scores_obs_>=doublet_threshold) / ncells_filtered)), "\n\n")
  seu = subset(seu, doublet_scores < doublet_threshold)
  
  seu.ls[[n]] = seu
}

# Merge a list of objects
seu = merge(x = seu.ls[[1]], y = seu.ls[-1], add.cell.ids = data.tib$Sample)
seu = JoinLayers(seu)
seu = RenameCells(seu, new.names = sub("BM_(.*)-.*", "\\1", colnames(seu)))

# Construct and add metadata
external_meta = annotations %>%
  rownames_to_column("Cell") %>%
  rename(CellType=type) %>%
  mutate(Sample = paste0("BM_", gsub("_.*", "", Cell)),
         Dataset = "Hourigan2018") %>%
  left_join(supplemental, by="Sample") %>%
  mutate(Cohort = case_when(Age<=40 ~ "Young", 
                            Age<60 ~ "Middle",
                            Age>=60 ~ "Aged"),
         Sample = gsub("[0-9]*", "", Sample))

meta = seu@meta.data
meta = meta %>% 
  rownames_to_column("Cell") %>%
  left_join(external_meta, by="Cell")
rownames(meta) = meta$Cell
seu@meta.data = meta

# Keep only stem and progenitor cells
seu = subset(seu, CellType %in% c("Early erythroid progenitors", "Monocyte progenitors", "HSPCs"))

# Save RDS
saveRDS(seu, file = "Hourigan_seu_formatted.rds")
