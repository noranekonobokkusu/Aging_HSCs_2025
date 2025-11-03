suppressMessages(suppressWarnings({
  library(tidyverse)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

# Read multiple selected Ks
selected_ks = readxl::read_xlsx("stable_ks.xlsx")
colnames(selected_ks) = c("Sample", "k")
selected_ks = selected_ks %>% separate_rows(k, sep = ",")
selected_ks$dt = sub("\\.", "_", ifelse(grepl("-", selected_ks$k), gsub(".*-", "", selected_ks$k), "0.1"))
selected_ks$k = gsub("-.+", "", selected_ks$k)
selected_ks = selected_ks %>%
  filter(k != "-")

# Get the full list of gene names, to fill missing values with 0
gene_names = MuDataSeurat::ReadH5AD(paste0(base_dir, "1_dataset_preparation/seu_hsc.var_genes.h5ad")) %>% rownames()
gene_names = gsub("-", ".", gene_names)

# Read Z-score spectra for selected Ks; exclude singletons based on the usage stats (similar to '3_plot_separate_usages_and_find_singletons.R')
spectra_score_files = paste0("Z_scores/", selected_ks$Sample, "/", selected_ks$k, ".", selected_ks$dt, ".txt")
usage_files = paste0("usages/", selected_ks$Sample, "/", selected_ks$k, ".", selected_ks$dt, ".txt")

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
    dplyr::select(-K, -Sample) %>%
    column_to_rownames("Program")
  
  # Identify and drop singletons
  usages = read.table(usage_files[i], sep="\t", header = T, row.names = 1)
  split_string = strsplit(f, "/")[[1]]
  sample = split_string[2]
  k = ncol(usages)
  
  df = usages %>%
    rownames_to_column("cell") %>%
    pivot_longer(-cell) %>%
    group_by(cell) %>%
    mutate(value = value / sum(value)) %>%
    mutate(name = factor(name, levels = paste0("X", 1:k))) %>%
    arrange(desc(name))
  
  singletons = df %>%
    group_by(name) %>%
    dplyr::rename("program"="name") %>%
    summarize(r = max(value)/quantile(value, 0.75)) %>%
    filter(r>10) %>%
    mutate(program = paste0(sample, "_K", k, ":", gsub("X", "", program))) %>%
    pull(program)
  
  t = t[! rownames(t) %in% singletons,]
  merged_t_hsc = rbind(merged_t_hsc,t)
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Write spectra to the file
write.table(merged_t_hsc, paste0("random_gene_sets_control/total_cNMF_output.csv"), quote = F, row.names = T, col.names = T, sep=",", append = F)
