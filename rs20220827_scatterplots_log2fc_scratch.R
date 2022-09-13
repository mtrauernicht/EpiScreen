# For ratios_tib_fc run lines 1-343 in 3_Synergies_Processing.rmd
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)


k27_drugs = CCD_tib %>% filter(!is.na(ratio_slope_plot) & ratio_slope_plot < -0.08 &  feature == "H3K27me3")  %>% distinct(drug_conc) %>% pull(drug_conc)
CCD_tib %>% filter(!is.na(ratio_slope_plot) & feature == "H3K27me3")  %>% distinct(drug_conc, ratio_slope_plot) %>% kable()


ratios_tib_fc %>% filter(drug_conc %in% k27_drugs) %>%
  distinct(IPR, drug_conc, target, log2.fc.ratio, H3K27me3) %>%
ggplot(.,aes(H3K27me3, log2.fc.ratio))  + geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")),
               label.x = "left", label.y = "bottom") +
  facet_wrap(target ~ drug_conc) +
  theme_bw(base_size = 16) 


ggsave("/DATA/projects/DSBrepair/git/EpiScreen/figures/rs20220826_K27me3_CCD_screen_scatterplots_all_tech.pdf", height = 24, width = 20)



m5c_drugs = CCD_tib %>% filter(!is.na(indelrate_slope_plot) &  feature == "m5C")  %>% distinct(drug_conc) %>% pull(drug_conc)
CCD_tib %>% filter(!is.na(indelrate_slope_plot) & feature == "m5C")  %>% distinct(drug_conc, indelrate_slope_plot) %>% kable()


ratios_tib_fc %>% filter(drug_conc %in% m5c_drugs & target == "DNMT") %>%
  distinct(IPR, drug_conc, target, log2.fc.freqCut, m5C) %>%
  ggplot(.,aes(m5C, log2.fc.freqCut))  + geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")),
               label.x = "left", label.y = "bottom") +
  facet_wrap(. ~ drug_conc) +
  theme_bw(base_size = 16) 


ggsave("/DATA/projects/DSBrepair/git/EpiScreen/figures/rs20220826_5mC_CCD_screen_scatterplots.pdf", height = 24, width = 20)

ratios_tib_fc %>% distinct(IPR, drug, conc_char, target, log2.fc.freqCut, log2.fc.ratio) %>%
  pivot_longer(cols = c("log2.fc.freqCut", "log2.fc.ratio")) %>% ggplot(., aes(value, fill = name, alpha = 0.5)) + geom_density()




m5c_drugs = CCD_tib %>% filter(!is.na(indelrate_slope_plot) &  feature == "m5C")  %>% distinct(drug_conc) %>% pull(drug_conc)
CCD_tib %>% filter(!is.na(indelrate_slope_plot) & feature == "m5C")  %>% distinct(drug_conc, indelrate_slope_plot) %>% kable()


ratios_tib_fc %>% filter(drug_conc == "Apicidin 10 µM") %>%
  pivot_longer(cols = chromatin.features.short) %>%
  distinct(IPR, drug_conc, target, log2.fc.freqCut, value, name) %>%
  ggplot(.,aes(value, log2.fc.freqCut))  + geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")),
               label.x = "left", label.y = "bottom") +
  facet_wrap(. ~ name) +
  theme_bw(base_size = 16) 

