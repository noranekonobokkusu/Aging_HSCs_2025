# Dataset preparation notes
Except for the Weng dataset, we used either raw sequencing data or count matrices (where available) as input.

For the Weng dataset, we obtained single-cell multi-ome data from two young and two aged donors (https://doi.org/10.6084/m9.figshare.23290004), supplemented with six additional young donors and three additional aged donors (Dr. Vijay Sankaran, personal communication), and combined gene expression counts from all 13 donors into a single object as described in the GitHub repository https://github.com/petervangalen/ReDeeM_2024/ (script 2_MergeSeuratObjectAllCells.Rmd). The object used as input in `6_weng_create_seurat_object.R` is available on [Figshare](https://figshare.com/projects/Aging_HSCs_2025/235781).
