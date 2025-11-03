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

# Read data
t = read.table(paste0(base_dir, "input_data/Adelman_data/adelman.total.gene.counts.txt"), sep = "\t", header = T)
meta_df = read.table(paste0(base_dir, "input_data/Adelman_data/adelman_cs_metadata.csv"), sep=",", header=T)

# Drop excessive headers
t = t %>%
  filter(GeneID != "GeneID")

colnames(t)[4] = "Cell"
t$Cell %>% unique() %>% length() # Should be 650 cells

# Convert into a matrix-like format
counts.mat = t %>%
  tidyr::pivot_wider(names_from = Cell, values_from = Counts) %>%
  dplyr::select(-Length) %>%
  column_to_rownames("GeneID")
counts.mat[1:3,1:3]

# Create Seurat object and append metadata
seu = CreateSeuratObject(counts = counts.mat, min.cells = 5, min.features = 100)

meta = seu@meta.data
meta$Cell = rownames(meta)
meta = meta %>%
  left_join(meta_df, by=c("Cell" = "SRR_ID")) %>%
  mutate(Dataset = "Adelman2019",
         Cohort = ifelse(Age<=40, "Young", "Aged"),
         Sex = ifelse(grepl("F", Donor), "F", "M"),
         Sample = gsub(".* ", "", gsub("-", "_", Donor)),
         CellType = celltype) %>%
  dplyr::select(-c(cell_title, celltype_general, Sample_description, ID_REF, Donor, celltype))
rownames(meta) = meta$Cell
seu@meta.data = meta
Idents(seu) = "Sample"

# Remove cells with high percentage of mitochondrial genes; Fluidigm data seems to be of a pretty high quality and low-throughput, so I didn't do any filtering based on nCount/nFeature/doublet scores (I tested scrubletR on this dataset and the scores looked unreasonable probably because the data is very homogeneous and cell counts per sample are low)
seu = PercentageFeatureSet(seu, "^MT-", col.name = "percent.mito")
seu$percent.mito = seu$percent.mito/100
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA", "percent.mito"))
seu = subset(seu, percent.mito < mito.high)

# Save formatted object
saveRDS(seu, "Adelman_seu_formatted.rds")
