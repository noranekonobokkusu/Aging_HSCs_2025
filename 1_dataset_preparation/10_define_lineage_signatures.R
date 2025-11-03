suppressMessages(suppressWarnings({
  library(Seurat)
  library(tidyverse)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "1_dataset_preparation"))

# Read HSPCs object and combine HSCs and MPPs
seu = readRDS("seu_full_object_hscs_marked.rds")
seu = NormalizeData(seu)
seu$predicted_CellType[seu$predicted_CellType %in% c("HSC", "MPP-MkEry", "MPP-MyLy")] = "HSC_MPP"

# Contrast HSC+MPP to committed populations
Idents(seu) = "predicted_CellType"
mep_markers = FindMarkers(seu, ident.1 = "MEP", ident.2 = "HSC_MPP")
gmp_markers = FindMarkers(seu, ident.1 = "Early GMP", ident.2 = "HSC_MPP")
clp_markers = FindMarkers(seu, ident.1 = "CLP", ident.2 = "HSC_MPP")

# Define preliminary signatures
mep_sig = mep_markers %>% 
  filter(p_val_adj < 1e-10, avg_log2FC > 1) %>% 
  rownames_to_column("gene") %>%
  slice_head(n=50) %>% 
  pull(gene)

gmp_sig = gmp_markers %>% 
  filter(p_val_adj < 1e-10, avg_log2FC > 1) %>% 
  rownames_to_column("gene") %>%
  slice_head(n=50) %>% 
  pull(gene)

clp_sig = clp_markers %>% 
  filter(p_val_adj < 1e-10, avg_log2FC > 1) %>% 
  rownames_to_column("gene") %>%
  slice_head(n=50) %>% 
  pull(gene)

# To make signatures more specific, let's exclude shared genes
overlap_genes = unique(c(intersect(mep_sig, gmp_sig), intersect(mep_sig, clp_sig), intersect(gmp_sig, clp_sig)))
mep_sig = setdiff(mep_sig, overlap_genes)
gmp_sig = setdiff(gmp_sig, overlap_genes)
clp_sig = setdiff(clp_sig, overlap_genes)

write.table(mep_sig, file = paste0(base_dir, "input_data/signatures/mep_sig.tsv"), append = F, quote = F, sep = "\t", row.names = T, col.names = T)
write.table(gmp_sig, file = paste0(base_dir, "input_data/signatures/gmp_sig.tsv"), append = F, quote = F, sep = "\t", row.names = T, col.names = T)
write.table(clp_sig, file = paste0(base_dir, "input_data/signatures/clp_sig.tsv"), append = F, quote = F, sep = "\t", row.names = T, col.names = T)
