library(tidyverse)
library(openxlsx)

#### Data ####
load("../AM_MDM_TB/results/gene_set/enrich_cell_TB.RData")

bind_rows(enrich_AM_h, enrich_MDM_h) %>% View
  #Significant pathways
  filter(FDR < 0.1) %>% 
  select(group, pathway, group_in_pathway, size_pathway, `k/K`, FDR) %>% 
  #Clean pathway names
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway),
         pathway = gsub("interferon gamma", "IFNG", pathway),
         pathway = gsub("interferon alpha", "IFNA", pathway),
         pathway = gsub("myc", "MYC", pathway),
         pathway = gsub("mtorc", "MTORC", pathway),
         pathway = gsub("il2 stat5", "IL2 STAT5", pathway),
         pathway = gsub("g2m", "G2M", pathway),
         pathway = gsub("e2f", "E2F", pathway),
         pathway = gsub("p53", "P53", pathway)) %>% 
  #Clean columns and data
  rename(`DEG in pathway (k)`=group_in_pathway,
         `Total genes in pathway (K)`=size_pathway,
         Pathway = pathway) %>% 
  mutate(FDR = signif(FDR, digits = 2),
         `k/K` = signif(`k/K`, digits = 2)) %>% 
  #separate groups
  separate(group, into = c("trash","Cell","fc")) %>% 
  #keep only induced terms
  filter(fc=="up") %>% 
  select(-trash,-fc) %>% 
  arrange(Cell, FDR) %>% 
  write.xlsx("Table1_DEG.enrichment.xlsx")
