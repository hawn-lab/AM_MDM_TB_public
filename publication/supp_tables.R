library(tidyverse)
library(openxlsx)

##### S1: AM/MDM linear models #####
df_list <- list()

load("../AM_MDM_TB/results/gene_level/kimma_lme_cell_TB.RData")
attach("../AM_MDM_TB/data_clean/AM-MDM_voom.RData")
dat.am.mdm <- dat.voom

am.mdm <- dat.am.mdm$genes %>% 
  left_join(AM.MDM_kimma_pval$lme.contrast, by=c("ensembl_gene_id"="gene")) %>% 
  select(ensembl_gene_id, symbol, 
         contrast_ref, contrast_lvl, estimate, pval, FDR) %>% 
  arrange(FDR)

contrast.OI <- data.frame(
  contrast_ref = c("AM Media", "AM Media",  "MDM Media", "AM TB"),
  contrast_lvl = c("AM TB",    "MDM Media", "MDM TB",    "MDM TB")
)

df_list[["TB in AM"]] <- am.mdm %>% 
  inner_join(contrast.OI[1,])
df_list[["TB in MDM"]] <- am.mdm %>% 
  inner_join(contrast.OI[3,])
df_list[["AM vs MDM in media"]] <- am.mdm %>% 
  inner_join(contrast.OI[2,])
df_list[["AM vs MDM in TB"]] <- am.mdm %>% 
  inner_join(contrast.OI[4,])

#Save
write.xlsx(df_list, "TableS1_AM.Mtb.linear.models.xlsx")

#### S4: MDM IFN linear model ####
df_list2 <- list()

load("../MDM_IFN/results/gene_level/kimma_IFN.RData")
attach("../MDM_IFN/data_clean/MDM-IFN_voom.RData")
dat.ifn <- dat.voom

mdm.ifn <- dat.ifn$genes %>% 
  left_join(MDM.IFN_kimma_pval$lme.contrast, by=c("ensembl_gene_id"="gene")) %>% 
  select(ensembl_gene_id, symbol, 
         contrast_ref, contrast_lvl, estimate, pval, FDR) %>% 
  filter(contrast_ref == "untreated") %>% 
  arrange(FDR)

df_list2[["IFNA8 in MDM"]] <- mdm.ifn %>% 
  filter(contrast_lvl == "IFNA8")
df_list2[["IFNb in MDM"]] <- mdm.ifn %>% 
  filter(contrast_lvl == "IFNb")
df_list2[["IFNe in MDM"]] <- mdm.ifn %>% 
  filter(contrast_lvl == "IFNe")
df_list2[["IFNg in MDM"]] <- mdm.ifn %>% 
  filter(contrast_lvl == "IFNg")
# df_list2[["IFNL1 in MDM"]] <- mdm.ifn %>% 
#   filter(contrast_lvl == "IFNL1")

#Save
write.xlsx(df_list2, "TableS5_MDM.IFN.linear.models.xlsx")

###### S2: Enrich ######
load("../AM_MDM_TB/results/gene_set/enrich_cell_TB.RData")

df_list3 <- list()

df_list3[["AM-specific DEG"]] <- enrich_AM_h %>% 
  select(group, size_group, pathway, group_in_pathway, size_pathway,
         `k/K`, pval, FDR, genes) %>% 
  rename(`DEG total`=size_group, `DEG in pathway (k)`=group_in_pathway,
         `Pathway total (K)`=size_pathway) %>% 
  unnest(genes) %>% 
  left_join(dat.am.mdm$genes %>% 
              select(ensembl_gene_id, symbol),
            by = c("genes"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(group, `DEG total`, pathway, `DEG in pathway (k)`, `Pathway total (K)`,
           `k/K`, pval, FDR) %>% 
  summarise(ensembl_gene_id = list(unique(genes)),
            symbol = list(unique(symbol))) %>% 
  mutate(ensembl_gene_id = as.character(ensembl_gene_id),
         symbol = as.character(symbol)) %>% 
  mutate(group = gsub("kimma_","",group))
  
df_list3[["MDM-specific DEG"]] <- enrich_MDM_h %>% 
  select(group, size_group, pathway, group_in_pathway, size_pathway,
         `k/K`, pval, FDR, genes) %>% 
  rename(`DEG total`=size_group, `DEG in pathway (k)`=group_in_pathway,
         `Pathway total (K)`=size_pathway) %>% 
  unnest(genes) %>% 
  left_join(dat.am.mdm$genes %>% 
              select(ensembl_gene_id, symbol),
            by = c("genes"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(group, `DEG total`, pathway, `DEG in pathway (k)`, `Pathway total (K)`,
           `k/K`, pval, FDR) %>% 
  summarise(ensembl_gene_id = list(unique(genes)),
            symbol = list(unique(symbol))) %>% 
  mutate(ensembl_gene_id = as.character(ensembl_gene_id),
         symbol = as.character(symbol)) %>% 
  mutate(group = gsub("kimma_","",group))

#Save
write.xlsx(df_list3, "TableS3_DEG.Hallmark.xlsx")

##### GSEA #####
df_list4 <- list()

load("../AM_MDM_TB/results/gene_set/gsea_cell_TB.RData")

gsea_h_kimma <- gsea_h_kimma %>% 
  select(group, pathway, ES, NES, pval, FDR, leadingEdge) %>% 
  separate(group, into=c("contrast_lvl","contrast_ref"), sep="-") %>% 
  mutate(across(c(contrast_lvl,contrast_ref), ~gsub("_"," ", .))) %>% 
  unnest(leadingEdge) %>% 
  left_join(dat.am.mdm$genes %>% 
              select(ensembl_gene_id, symbol),
            by = c("leadingEdge"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(contrast_lvl, contrast_ref, pathway, ES, NES, pval, FDR) %>% 
  summarise(ensembl_gene_id = list(unique(leadingEdge)),
            symbol = list(unique(symbol))) %>% 
  mutate(ensembl_gene_id = as.character(ensembl_gene_id),
         symbol = as.character(symbol)) %>% 
  arrange(FDR)

df_list4[["TB in AM"]] <- gsea_h_kimma %>% 
  inner_join(contrast.OI[1,])
df_list4[["TB in MDM"]] <- gsea_h_kimma %>% 
  inner_join(contrast.OI[3,])
df_list4[["AM vs MDM in media"]] <- gsea_h_kimma %>% 
  inner_join(contrast.OI[2,])
df_list4[["AM vs MDM in TB"]] <- gsea_h_kimma %>% 
  inner_join(contrast.OI[4,])

#Save
write.xlsx(df_list4, "TableS4_GSEA.xlsx")

#### S5: MDM IFN Enrich ####
load("../MDM_IFN/results/gene_set/enrich_MDM_IFN.RData")
attach("../MDM_IFN/data_clean/MDM-IFN_voom.RData")

# Gene sets specific to only 1 IFN, FDR < 0.05.
venn.ls5 <- list()

for(g in unique(enrich_IFN_h$group)){
  temp.signif <- enrich_IFN_h %>% 
    filter(group==g & FDR < 0.05 & group_in_pathway > 1) %>% 
    pull(pathway) %>% unique()
  venn.ls5[[g]] <- temp.signif
}

#Unique IFNA8 GO terms
IFNA_signif <- 
  venn.ls5[["IFNA8-untreated"]][!venn.ls5[["IFNA8-untreated"]]%in%unlist(venn.ls5[2:length(venn.ls5)])]

temp <- enrich_IFN_h %>% 
  filter(pathway %in% IFNA_signif) %>% 
  unnest(genes) %>% 
  left_join(dat.voom$genes %>% select(ensembl_gene_id,symbol),
            by=c("genes"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(group, size_group, pathway, group_in_pathway, size_pathway, `k/K`,
           pval, FDR) %>% 
  summarise(ensembl_gene_id = list(unique(genes)),
            symbol = list(unique(symbol)),
            .groups="drop") %>% 
  rename(`DEG total`=size_group, `DEG in pathway (k)`=group_in_pathway,	
         `Pathway total (K)`=size_pathway) %>% 
  mutate(group = gsub("-untreated","",group)) 


temp %>% 
  mutate(across(ensembl_gene_id:symbol, ~as.character(.))) %>% 
  write.xlsx(., file = "TableS6_IFNA_enrich_unique.xlsx")

##### AM and MDM specific DEG #####
deg <- read_csv("../AM_MDM_TB/results/gene_level/AM-MDM_kimma.venn.signif.ensembl.csv") %>%
  rownames_to_column() %>%
  pivot_longer(-rowname,names_to="DEG_group", values_to="ensembl_gene_id") %>%
  drop_na(ensembl_gene_id) %>%
  left_join(dat.am.mdm$genes) %>%
  mutate(`AM-Mtb-specific` = ifelse(DEG_group %in% c("AM.down.only","AM.up.only",
                                                     "AM.down_MDM.up","AM.up_MDM.down"),
                                    "Y","N"),
         `MDM-Mtb-specific` = ifelse(DEG_group %in% c("MDM.down.only","MDM.up.only",
                                                      "AM.down_MDM.up","AM.up_MDM.down"),
                                    "Y","N")) %>% 

  select(DEG_group, ensembl_gene_id, symbol, `AM-Mtb-specific`, `MDM-Mtb-specific`) %>%
  arrange(DEG_group) %>%
  unnest(symbol)

#Add FDRs
def_fdr <- am.mdm %>%
  inner_join(contrast.OI) %>% 
  mutate(name = case_when(contrast_ref=="AM Media" & contrast_lvl=="AM TB"~
                            "TB in AM",
                          contrast_ref=="MDM Media" & contrast_lvl=="MDM TB"~
                            "TB in MDM",
                          contrast_ref=="AM Media" & contrast_lvl=="MDM Media"~
                            "AM vs MDM in media",
                          contrast_ref=="AM TB" & contrast_lvl=="MDM TB"~
                            "AM vs MDM in TB")) %>% 
  select(ensembl_gene_id, name, FDR) %>% 
  pivot_wider(values_from = FDR) %>% 
  right_join(deg) %>% 
  select(DEG_group,ensembl_gene_id,symbol, 
         contains("specific"), everything())

df_list5 <- list()
for(group in sort(unique(def_fdr$DEG_group))){
  df_list5[[group]] <- def_fdr %>% 
    filter(DEG_group==group) %>% 
    arrange(symbol)
}

write.xlsx(df_list5, "TableS2_AM.MDM.specific.DEG.xlsx")

deg %>% count(DEG_group)
