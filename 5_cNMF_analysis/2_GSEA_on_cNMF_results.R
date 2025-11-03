suppressMessages(suppressWarnings({
  library(fgsea)
  library(Seurat)
  library(dplyr)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis"))

# Create a directory if it doesn't exist yet
if (!dir.exists("gsea_results")) dir.create("gsea_results")

##### Create a list of pathways #####

# Aging signature
de_results = read.table(paste0(base_dir, "2_differential_expression_analysis/DESeq2_results.Aged_vs_Young.csv"), sep=",", header = T)
aging_sig = de_results %>%
  filter(padj<0.05, log2FoldChange>1) %>%
  pull(gene)

# Human dormancy signature from https://doi.org/10.1016/j.stem.2021.07.003
garcia_table = readxl::read_xlsx(paste0(base_dir, "input_data/signatures/mmc2.xlsx"), sheet = 1)
garcia_signature = garcia_table %>%
  filter(genename!="None") %>%
  filter(logFC < -2,FDR<0.00001) %>%
  pull(genename) %>%
  unique()

graham_signature = read.table(paste0(base_dir, "input_data/signatures/graham_quiescence_up.txt"), header = F)$V1

mep_sig = read.table(paste0(base_dir, "input_data/signatures/mep_sig.tsv"))$x
gmp_sig = read.table(paste0(base_dir, "input_data/signatures/gmp_sig.tsv"))$x
clp_sig = read.table(paste0(base_dir, "input_data/signatures/clp_sig.tsv"))$x

ap1_complex = c("FOS", "FOSB", "JUN", "JUNB", "JUND", "FOSL1", "FOSL2", "ATF2", "ATF3", "ATF4", "ATF6", "ATF7", "BATF", "BATF2", "BATF3", "MAFA", "MAFB", "MAF", "MAFG", "MAFF", "MAFK")

main_pathways = list(graham_signature, garcia_signature, aging_sig, mep_sig, gmp_sig, clp_sig, ap1_complex)

names(main_pathways) = c("Graham_quiescence", "Garcia_quiescence", "Aging_signature", "MEP_signature", "GMP_signature", "CLP_signature", "AP-1_complex")

##### Run fgsea for each cNMF program separately #####

programs = read.table(paste0("cnmf_output/Z_scores.total_cNMF_output.csv"), row.names = 1, header = T, sep=",")
# Initialize results storage
fgseaRes_list = list()

{
# Set up progress bar
cat("Running GSEA..\n")
pb = txtProgressBar(min = 0, max = nrow(programs), style = 3)

for (i in 1:nrow(programs)) {
  entry = programs[i,]
  s = entry$Sample
  entry[c("Sample")] = NULL
  ranks = as.numeric(entry)
  names(ranks) = colnames(entry)
  
  # GSEA
  fgseaRes = fgsea(pathways = main_pathways, 
                   stats    = ranks,
                   minSize  = 5,
                   eps      = 0.0,
                   maxSize  = 500,
                   nPermSimple = 10000)
  
  fgseaRes$Program = rownames(entry)
  fgseaRes$Sample = s
  
  # Store results in lists
  fgseaRes_list[[i]] = fgseaRes
  
  setTxtProgressBar(pb, i)
}

close(pb)
}

# Combine results from lists
fgseaRes_total = do.call(rbind, fgseaRes_list)
# Save
data.table::fwrite(fgseaRes_total, file=paste0("gsea_results/total_cNMF_fgseaRes.txt"), sep="\t", sep2=c("", " ", ""), append = F)
