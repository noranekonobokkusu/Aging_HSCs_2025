suppressMessages(suppressWarnings({
  library(tidyverse)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

# Create a directory if it doesn't exist yet
if (!dir.exists("cnmf_output")) dir.create("cnmf_output")

# Read selected Ks
selected_ks = readxl::read_xlsx("stable_single_k.xlsx")
colnames(selected_ks) = c("Sample", "k")
selected_ks$dt = ifelse(grepl("-", selected_ks$k), gsub(".*-", "", selected_ks$k), "0_1")
selected_ks$k = gsub("-.+", "", selected_ks$k)
selected_ks = selected_ks %>%
  filter(k != "-")

gene_names = MuDataSeurat::ReadH5AD(paste0(base_dir, "1_dataset_preparation/seu_hsc.var_genes.h5ad")) %>% rownames()
gene_names = gsub("-", ".", gene_names)

# Read and process cNMF output files for relevant Ks
for (data_type in c("TPMs", "Z_scores")) {
  
  cat(sprintf("Processing: %s", data_type), "\n")
  spectra_score_files = paste0(data_type, "/", selected_ks$Sample, "/", selected_ks$k, ".", selected_ks$dt, ".txt")
  
  merged_t_hsc = data.frame()
  
  pb = txtProgressBar(min = 0, max = length(spectra_score_files), style = 3)
  
  for (i in 1:length(spectra_score_files)) {
    
    f = spectra_score_files[i]
    t = read.table(f, sep="\t", header = T, row.names = 1)
    
    # In some samples, some of the genes were excluded from cNMF; this completes the matrix to the original gene set
    missing_genes = setdiff(gene_names, colnames(t))
    for (gene in missing_genes) {
      t[[gene]] = 0
    }
    t = t[, gene_names]
    
    split_string = strsplit(f, "/")[[1]]
    s = split_string[2]
    k = gsub("\\..*", "", split_string[3])
    
    t[is.na(t)] = 0
    t$Sample = s
    t$K = k
    t$Program = rownames(t)
    t = t %>%
      mutate(Program = paste0(Sample, "_K", K, ":", Program)) %>%
      dplyr::select(-K) %>%
      column_to_rownames("Program")
    
    merged_t_hsc = rbind(merged_t_hsc,t)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  write.table(merged_t_hsc, paste0("cnmf_output/", data_type, ".total_cNMF_output.csv"), quote = F, row.names = T, col.names = T, sep=",", append = F)
  
}
