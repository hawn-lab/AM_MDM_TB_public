library(tidyverse)

attach("data_clean/MDM-IFN_voom.RData")
attach("results/gene_level/kimma_IFN.RData")

genes.OI <- c("PLA2G3","APOC3","KCNQ1")

dat <- MDM.IFN_kimma_pval$lme.contrast %>% 
  filter(contrast_ref == "untreated") %>% 
  left_join(dat.voom$genes, by=c("gene"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  filter(grepl("APOC", symbol) | grepl("PLA2G3", symbol) | grepl("KCNQ1", symbol)) %>% 
  select(contrast_lvl, gene, symbol, estimate, pval, FDR) %>% 
  filter(FDR < 0.2)
