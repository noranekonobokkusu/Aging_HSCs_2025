suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
}))

# This code plots usage heatmaps for manual inspection, and identifies likely singletons (programs driven by a single cell)

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"
setwd(paste0(base_dir, "5_cNMF_analysis/"))

if (!dir.exists("usage_plots")) dir.create("usage_plots")

selected_ks = readxl::read_xlsx("stable_single_k.xlsx")
colnames(selected_ks) = c("Sample", "k")
selected_ks$dt = ifelse(grepl("-", selected_ks$k), gsub(".*-", "", selected_ks$k), "0_1")
selected_ks$k = gsub("-.+", "", selected_ks$k)
selected_ks = selected_ks %>%
  filter(k != "-")
usage_files = paste0("usages/", selected_ks$Sample, "/", selected_ks$k, ".", selected_ks$dt, ".txt")
if (!dir.exists("usage_plots/")) dir.create("usage_plots/")

singleton_programs = c()
for (i in 1:length(usage_files)) {
  f = usage_files[i]
  usages = read.table(f, sep="\t", header = T, row.names = 1)
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
  
  p1 = ggplot(df, aes(x=cell,y=name,fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white", high = "red3") +
    theme_pubr() +
    theme(axis.text.x = element_blank(),
          legend.position = "none") +
    ggtitle(sample)
  
  p2 = df %>%
    group_by(name) %>%
    dplyr::rename("program"="name") %>%
    summarize(max = max(value),
              med = median(value),
              q75 = quantile(value, 0.75),
              q90 = quantile(value, 0.9),
              r75 = max/q75,
              r9 = max/q90) %>%
    pivot_longer(-program) %>%
    mutate(value = round(value,2)) %>%
    ggplot(aes(x=name,y=program,label=value,color=value>10)) +
    geom_text() +
    scale_color_manual(values = c("black", "red")) +
    theme_pubr() +
    theme(axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none") +
    ggtitle("")
  
  p = ggarrange(p1, p2, widths = c(0.65, 0.25), vjust = T)
  cowplot::save_plot(paste0("usage_plots/", sample, ".png"), p, base_height = 5, base_width = 12)
  
  singletons = df %>%
    group_by(name) %>%
    dplyr::rename("program"="name") %>%
    summarize(r = max(value)/quantile(value, 0.75)) %>%
    filter(r>10) %>%
    mutate(program = paste0(sample, "_K", k, ":", gsub("X", "", program))) %>%
    pull(program)
  
  singleton_programs = c(singleton_programs, singletons)
}

write.table(singleton_programs, paste0("singletons.txt"), append = F, quote = F, row.names = F, col.names = F)
