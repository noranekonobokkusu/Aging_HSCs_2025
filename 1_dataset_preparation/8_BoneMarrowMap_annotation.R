suppressMessages(suppressWarnings({
  library(tidyverse)
  library(Seurat)
  library(BoneMarrowMap)
  library(symphony)
  library(cowplot)
  library(MuDataSeurat)
  library(RColorBrewer)
  library(ggpubr)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "1_dataset_preparation/"))

#### Read count matrices and define genes present in all datasets ####
files = list.files(".", pattern = "_seu_formatted\\.rds$", full.names = TRUE)
seu_ls = list()
gene_names_ls = list()

for (i in 1:length(files)) {
  seu = readRDS(files[[i]])
  seu_ls[[i]] = seu
  gene_names_ls[[i]] = rownames(seu_ls[[i]])
}

common_genes = Reduce(intersect, gene_names_ls) #~11k genes

#### Map each of the datasets onto BMM reference; this follows BMM tutorial ####

# Load Symphony reference (pre-downloaded)
ref = readRDS(paste0(base_dir, "input_data/BoneMarrowMap/BoneMarrow_RefMap_SymphonyRef.rds"))
# Set uwot path for UMAP projection
ref$save_uwot_path = paste0(base_dir, "input_data/BoneMarrowMap/BoneMarrow_RefMap_uwot_model.uwot")

# Save annotated datasets into a list
seu_ls_annotated = list()

for (i in 1:length(seu_ls)) {
  
  cat(sprintf("Processing sample %s, %d/%d", unique(seu_ls[[i]]$Dataset), i, length(seu_ls)), "\n")
  
  counts = as.matrix(GetAssayData(seu_ls[[i]], layer = "counts"))
  counts = counts[rownames(counts) %in% common_genes,]
  query = CreateSeuratObject(counts=counts, meta.data = seu_ls[[i]]@meta.data)
  query = NormalizeData(query)
  batchvar = "Sample"
  
  # Map query dataset using Symphony 
  query = map_Query(
    exp_query = query@assays$RNA$counts,
    metadata_query = query@meta.data,
    ref_obj = ref,
    vars = batchvar
  )
  
  # Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
  query = query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 
  
  # Predict Hematopoietic Cell Types by KNN classification
  query = predict_CellTypes(
    query_obj = query, 
    ref_obj = ref, 
    initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
    final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
  ) 
  
  query = AddMetaData(query, query@reductions$umap@cell.embeddings)
  
  seu_ls_annotated[[i]] = query
  gc()
  
}

# Save original BMM transfer output
saveRDS(seu_ls_annotated, file = "seu.bmm_annotation_transferred.rds")

#### Merge the list of Seurat objects into a single object ####

# seu_ls_annotated = readRDS("seu.bmm_annotation_transferred.rds")

new_seu_ls = list()

for (i in 1:length(seu_ls_annotated)) {
  print(i)
  seu = seu_ls_annotated[[i]]
  meta = seu@meta.data[c("Cell", "Cohort", "Sample", "Dataset", "Sex", "mapping_error_score", "predicted_CellType_prob", "predicted_CellType", "Age", "umap_1", "umap_2")]
  counts = as.matrix(GetAssayData(seu, layer = "counts"))
  seu = CreateSeuratObject(counts=counts, meta.data = meta)
  new_seu_ls[[i]] = seu
}

seu_total = merge(x = new_seu_ls[[1]], y = new_seu_ls[-1])
seu_total = JoinLayers(seu_total)

umap_coords = as.matrix(seu_total@meta.data[, c("umap_1", "umap_2")])

# Ensure row names match the Seurat object cell names
rownames(umap_coords) == colnames(seu_total)

# Create a UMAP reduction object
seu_total[["umap"]] = CreateDimReducObject(
  embeddings = umap_coords, 
  key = "umap_", 
  assay = DefaultAssay(seu_total)
)

seu_total@meta.data[c("umap_1", "umap_2")] = seu_total@reductions$umap@cell.embeddings

# Plot BMM-mapped UMAP
p = DimPlot(seu_total, group.by = "predicted_CellType", reduction = "umap", label = T, repel = T, raster.dpi = c(1600,1600), pt.size = 2) &
  NoLegend() &
  scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired"), brewer.pal(12, "Paired")))
save_plot("predicted_CellType_umap.pdf", p, base_height = 10, base_width = 10)

# I want to drop HSC-annotated cells outside the expected HSC area in the UMAP embedding. This would correspond to selecting HSCs with umap_1>0 & umap_1<4 & umap_2>1 & umap_2<5
seu_total$is_hsc = ifelse(( seu_total$predicted_CellType=="HSC" &
                             seu_total$umap_1>0 & 
                             seu_total$umap_1<4 & 
                             seu_total$umap_2>1 & 
                             seu_total$umap_2<5 ), "HSC", "non-HSC")


# Highlight HSC cluster for Fig.1a
p = seu_total@meta.data %>%
  filter(! is.na(is_hsc)) %>%
  mutate(is_hsc = factor(is_hsc, levels = c("HSC", "non-HSC"))) %>%
  arrange(desc(is_hsc)) %>%
  ggplot(aes(x=umap_1,y=umap_2,color=is_hsc)) +
  ggrastr::geom_point_rast(size=0.02,
                  raster.dpi=1000,
                  dev="ragg_png") +
  theme_pubr() + 
  coord_fixed() +
  theme_pubr(base_size = 5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("steelblue", "grey", "grey"))
p
save_plot("Fig.1a.pdf", p, base_height = 2, base_width = 2, dpi=1000)

# Drop unannotated cells
seu_total = subset(seu_total, predicted_CellType %in% unique(seu_total$predicted_CellType)[! is.na(unique(seu_total$predicted_CellType))])

# Drop UMAP coordinates from metadata
seu_total@meta.data[c("umap_1", "umap_2")] = NULL

saveRDS(seu_total, "seu_full_object_hscs_marked.rds")

#### Save as h5ad object for scVI integration ####

# Identify variable genes and create a subsetted object
counts = as.matrix(GetAssayData(seu_total, layer = "counts"))
seu_total$Sample = paste0(seu_total$Sample, "_", seu_total$Dataset)
seu_total[["RNA"]] = split(seu_total[["RNA"]], f = seu_total$Sample)
seu_total = NormalizeData(seu_total)
seu_total = FindVariableFeatures(seu_total, nfeatures = 3000)
variable_genes = VariableFeatures(seu_total)
seu_total = CreateSeuratObject(counts=counts[variable_genes,], meta.data = seu_total@meta.data)

# Check that metadata has no NAs
seu_total@meta.data %>% filter(if_any(everything(), is.na))

# Save the entire object as an H5AD file
WriteH5AD(seu_total, "seu_full_object.h5ad", assay="RNA")
