suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
}))

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"

setwd(paste0(base_dir, "revision"))

#### Compare percentages/MFIs of Young vs Aged samples ####

df = readxl::read_xlsx("experimental_validation/experimental_validation_young_vs_aged.xlsx")

df$plot = paste0(df$Population, "_", df$Measurement)
df$Population[is.na(df$Population)] = ""
df$Population = gsub("CD38- ", "CD38-\n", df$Population)

for (p in unique(df$plot)) {
  
  df_p = df %>%
    filter(plot == p)
  
  stat.test = df_p %>% 
    filter(to_include == 1) %>%
    wilcox_test(Value ~ Cohort) %>%
    add_xy_position()
  
  df_p = df_p %>%
    mutate(Cohort = factor(Cohort, levels=c("Young", "Aged")),
           to_include = factor(to_include))
  y_label = unique(df_p$Measurement)
  plot_title = unique(df_p$Population)
  y_range = diff(range(df_p$Value))
  
  p = ggplot() +
    geom_boxplot(data = df_p %>% filter(to_include==1),
                 aes(x=Cohort, y=Value),
                 outlier.shape = NA, 
                 linewidth=0.5/2.83465, 
                 width=0.3,
                 color="black") +
    geom_jitter(data = df_p,
                aes(x=Cohort, y=Value, color=to_include),
                width = 0.1, 
                size=0.3) +
    scale_color_manual(values = c("red", "black")) +
    theme_pubr(base_size = 5) +
    scale_y_continuous(name = y_label, limits = c(min(min(df_p$Value)*1.1, 0), max(df_p$Value) * 1.2 + 1), expand = c(0,0)) +
    ggtitle(plot_title) +
    stat_pvalue_manual(stat.test, 
                       label = "p={p}", 
                       tip.length = 0.03, 
                       bracket.nudge.y = max(df_p$Value) * 0.1,
                       size = 5/2.83465, 
                       bracket.size = 0.3/2.83465) +
    theme(legend.position = "none",
          axis.line = element_line(linewidth=0.6/2.83465, color="black"),
          axis.ticks = element_line(linewidth=0.5/2.83465, color="black"),
          plot.title = element_text(size=6)
    )
  
  save_plot(gsub("%", "percent", gsub(" ", "_", paste0("experimental_validation/", y_label, "_", plot_title, ".pdf"))), p, base_height = 40, base_width = 25, units = "mm")
  
}

#### Compare Ki-67+ fractions in p-c-JUN-positive and -negative populations ####

ki_67_in_p_cjun_df = readxl::read_xlsx("experimental_validation/experimental_validation_Ki67_in_p_c_Jun.xlsx")
colnames(ki_67_in_p_cjun_df) = c("Sample", "Q1", "Q2", "Q3", "Q4")

ki_67_in_p_cjun_df$to_include = c(1,1,1,1,0,1,1,1)
ki_67_in_p_cjun_df$Cohort = c("Young", "Young", "Young", "Aged", "Young", "Aged", "Aged", "Aged" )

ki_67_in_p_cjun_df = ki_67_in_p_cjun_df %>%
  mutate(ki_67_in_cJun_pos = Q2/(Q1+Q2),
         ki_67_in_cJun_neg = Q3/(Q3+Q4)) %>%
  pivot_longer(c(ki_67_in_cJun_pos, ki_67_in_cJun_neg))

stat.test = ki_67_in_p_cjun_df %>% 
  filter(to_include == 1) %>%
  wilcox_test(value ~ name, paired = TRUE) %>%
  add_xy_position()

ki_67_in_p_cjun_df = ki_67_in_p_cjun_df %>%
  mutate(to_include = factor(to_include))

p = ggplot() +
  geom_boxplot(data = ki_67_in_p_cjun_df %>% filter(to_include==1),
               aes(x=name, y=value),
               outlier.shape = NA, 
               linewidth=0.5/2.83465, 
               width=0.3,
               color="black") +
  geom_jitter(data = ki_67_in_p_cjun_df,
              aes(x=name, y=value, color=to_include),
              width = 0.1, 
              size=0.3) +
  scale_color_manual(values = c("red", "black")) +
  theme_pubr(base_size = 5) +
  scale_y_continuous(limits = c(-0.01, 0.5), 
                     expand = c(0,0),
                     name = c("%Ki-67+")) +
  scale_x_discrete(labels = c("p-c-JUN-", "p-c-JUN+"),
                   name = "Population") +
  ggtitle("CD34+ CD38-") +
  stat_pvalue_manual(stat.test, 
                     label = "p={p}", 
                     tip.length = 0.03, 
                     bracket.nudge.y = max(ki_67_in_p_cjun_df$value) * 0.1,
                     size = 5/2.83465, 
                     bracket.size = 0.3/2.83465) +
  theme(legend.position = "none",
        axis.line = element_line(linewidth=0.6/2.83465, color="black"),
        axis.ticks = element_line(linewidth=0.5/2.83465, color="black"),
        plot.title = element_text(size=6)
  )
p
save_plot("experimental_validation/KI67_pos_in_p_cJun_populations.pdf", p, base_height = 40, base_width = 25, units = "mm")
