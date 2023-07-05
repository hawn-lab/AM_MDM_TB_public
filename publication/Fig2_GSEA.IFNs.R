library(tidyverse)
library(limma)
library(rstatix)
library(ggpubr)
library(patchwork)

#### GSEA Data ####
load("../AM_MDM_TB/results/gene_set/gsea_cell_TB.RData")

fdr.cut <- 0.05
fdr.cut.label <- paste("FDR <", fdr.cut)

#### Define signif terms ####
cell.signif <- gsea_h_kimma %>% 
  filter(FDR <= fdr.cut & group %in% c("MDM_TB-AM_TB", "MDM_Media-AM_Media")) %>%
  pull(pathway) %>% unique()

tb.signif <- gsea_h_kimma %>% 
  filter(FDR <= fdr.cut & group %in% c("MDM_TB-MDM_Media", "AM_TB-AM_Media")) %>%
  pull(pathway)%>% unique()

h.to.plot <- intersect(cell.signif, tb.signif)

#### Format labels ####
dat.to.plot <- gsea_h_kimma %>% 
  #Significant terms
  filter(pathway %in% c(h.to.plot, "HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>% 
  mutate(Significance = ifelse(FDR <= fdr.cut, fdr.cut.label, "NS")) %>% 
  #Beautify labels
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway),
         pathway = gsub("interferon gamma", "IFNG", pathway),
         pathway = gsub("interferon alpha", "IFNA", pathway),
         pathway = gsub("mtorc", "MTORC", pathway),
         pathway = gsub("tnfa", "TNFA", pathway),
         pathway = gsub("nfkb", "NFKB", pathway)) %>% 
  mutate(group2 = recode_factor(factor(group),
                               "MDM_Media-AM_Media"="    AM <--  --> MDM  \nin Media",
                               "MDM_TB-AM_TB"="    AM <--  --> MDM  \nin +Mtb",
                               "AM_TB-AM_Media"="Media <--  --> +Mtb \nin AM",
                               "MDM_TB-MDM_Media"="Media <--  --> +Mtb \nin MDM")) %>% 
  #group and order terms
  mutate(facet_group = case_when(
    #Same direction, difference signif
    pathway %in% c("allograft rejection","MTORC1 signaling","estrogen response early","myogenesis") ~ "iii",
    #different directions
    pathway %in% c("inflammatory response","TNFA signaling via NFKB","epithelial mesenchymal transition","hypoxia") ~ "ii",
    pathway %in% c("IFNA response", "IFNG response") ~ "i"))

#### GSEA plot ####
h.plot1 <- dat.to.plot %>% 
  filter(group %in% c("MDM_Media-AM_Media", "MDM_TB-AM_TB")) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_segment( aes(reorder(pathway, NES), 
                    xend=pathway, y=0, yend=NES)) +
  geom_point(size=2.5, aes(fill = Significance),
             shape=21, stroke=1) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=c("FDR < 0.05" = "#a50026", 
                             "NS"="grey"))+
  lims(y=c(-3.5,3.5)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score") + 
  facet_grid(facet_group ~ group2, scales = "free", space = "free") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank())

# h.plot1

#### Gene Data ####
load("../AM_MDM_TB/data_clean/AM-MDM_voom.RData")
load("../AM_MDM_TB/results/gene_level/kimma_lme_cell_TB.RData")
## Add ABC to panels
col.vec <- c("#117733","#44AA99","#A46199","#882255")

#Format expression data into 1 df
contrast.OI <- data.frame(
  contrast_ref = c("AM Media", "AM Media",  "MDM Media", "AM TB"),
  contrast_lvl = c("AM TB",    "MDM Media", "MDM TB",    "MDM TB")
)

dat.format <- as.data.frame(dat.voom$E) %>% 
  #get gene symbols
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(select(dat.voom$genes, ensembl_gene_id, symbol)) %>% 
  unnest(symbol) %>% 
  #add and format metadata
  pivot_longer(-c(ensembl_gene_id,symbol), names_to = "libID") %>% 
  left_join(dat.voom$targets) %>% 
  mutate(TB = recode(TB, "TB"="+Mtb","Media"="Media")) %>% 
  mutate(x = paste(TB,cell,sep="\n"),
         x = factor(x, levels=c("Media\nAM","+Mtb\nAM",
                                "Media\nMDM","+Mtb\nMDM"))) %>% 
  select(ensembl_gene_id, symbol, libID, x, value, cell, TB)

#maxima to use as y.positions
dat.summ <- dat.format %>% 
  group_by(symbol) %>% 
  summarise(maxE=max(value, na.rm=TRUE),
            minE = min(value, na.rm=TRUE),
            .groups = "drop")

#Format lme results
fdr.format <- AM.MDM_kimma_pval$lme.contrast %>% 
  #Signif contrasts in genes of interest
  inner_join(contrast.OI) %>% 
  filter(FDR < 0.05) %>% 
  #get gene symbols
  left_join(select(dat.voom$genes, ensembl_gene_id, symbol), 
            by=c("gene"="ensembl_gene_id")) %>%
  unnest(symbol) %>% 
  select(symbol, contrast_ref, contrast_lvl, FDR) %>% 
  #Create start/stop groups for FDR bars
  mutate(group1 = recode(contrast_ref, 
                         "AM Media"="Media\nAM", "AM TB"="+Mtb\nAM",
                         "MDM Media"="Media\nMDM", "MDM TB"="+Mtb\nMDM"),
         group2 = recode(contrast_lvl, 
                         "AM Media"="Media\nAM", "AM TB"="+Mtb\nAM",
                         "MDM Media"="Media\nMDM", "MDM TB"="+Mtb\nMDM")) 

#### IFN genes ####
#genes of interest
genes <- c("IFNA1","IFNA8","IFNA13","IFNB1","IFNE","IFNG","IFNK","IFNL1")
dat.subset <- dat.format %>% 
  filter(symbol %in% genes) %>% 
  #Classes
  mutate(class = case_when(symbol %in% c("IFNA1","IFNA8","IFNA13","IFNB1","IFNE","IFNK") ~"Type I IFN",
                           symbol %in% c("IFNG") ~ "Type II IFN",
                           symbol %in% c("IFNL1") ~ "Type III IFN")) 
fdr.subset <- fdr.format %>% filter(symbol %in% genes)

#Format fdr for plotting
stat.dat <- fdr.subset %>%  
  #add y position
  left_join(dat.summ) %>% 
  group_by(symbol) %>% 
  mutate(rnk = order(FDR, decreasing=FALSE)) %>% 
  ungroup() %>% 
  mutate(y.position = maxE+abs(maxE*0.1*rnk)) %>%
  #Format FDR labels
  add_significance("FDR", cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("**", "*", "ns")) %>% 
  left_join(distinct(dat.subset, symbol, class))

type1 <- dat.subset %>% 
  filter(symbol %in% c("IFNA1","IFNA8","IFNA13","IFNB1","IFNE","IFNK")) %>% 
  ggplot(aes(x=x, y=value, color=x)) +
  #Raw data points
  geom_point(color = "grey80", position = position_dodge(width=0.75)) +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5) +
  #Beautify
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        strip.background = element_rect(linewidth=0.5)) +
  scale_color_manual(values = col.vec) +       
  labs(y="Normalized log2 expression", x="", color="",
       subtitle = "Type I IFN") +
  facet_wrap(~symbol, scales="free", nrow=1) +
  #Significant bars
  stat_pvalue_manual(filter(stat.dat, symbol %in% c("IFNA1","IFNA8","IFNA13","IFNB1","IFNE")),
                     label = "FDR.signif", tip.length = 0.02) +
  # clean x axis labels
  scale_x_discrete(labels=c("Media\nAM"="Media",
                            "+Mtb\nAM"="+Mtb",
                            "Media\nMDM"="Media", 
                            "+Mtb\nMDM"="+Mtb")) +
  scale_y_continuous(expand = c(0.15,0))

type2 <- dat.subset %>% 
  filter(symbol %in% c("IFNG")) %>% 
  ggplot(aes(x=x, y=value, color=x)) +
  #Raw data points
  geom_point(color = "grey80", position = position_dodge(width=0.75)) +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5) +
  #Beautify
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        strip.background = element_rect(linewidth=0.5)) +
  scale_color_manual(values = col.vec) +       
  labs(y="Normalized log2 expression", x="", color="",
       subtitle = "Type II IFN") +
  facet_wrap(~symbol, scales="free", nrow=1) +
  #Significant bars
  stat_pvalue_manual(filter(stat.dat, symbol %in% c("IFNG")),
                     label = "FDR.signif", tip.length = 0.02) +
  # clean x axis labels
  scale_x_discrete(labels=c("Media\nAM"="Media",
                            "+Mtb\nAM"="+Mtb",
                            "Media\nMDM"="Media", 
                            "+Mtb\nMDM"="+Mtb")) +
  scale_y_continuous(expand = c(0.15,0))

type3 <- dat.subset %>% 
  filter(symbol %in% c("IFNL1")) %>% 
  ggplot(aes(x=x, y=value, color=x)) +
  #Raw data points
  geom_point(color = "grey80", position = position_dodge(width=0.75)) +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5) +
  #Beautify
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        strip.background = element_rect(linewidth=0.5)) +
  scale_color_manual(values = col.vec) +       
  labs(y="Normalized log2 expression", x="", color="",
       subtitle = "Type III IFN") +
  facet_wrap(~symbol, scales="free", nrow=1) +
  #Significant bars
  stat_pvalue_manual(filter(stat.dat, symbol %in% c("IFNL1")),
                     label = "FDR.signif", tip.length = 0.02) +
  # clean x axis labels
  scale_x_discrete(labels=c("Media\nAM"="Media",
                            "+Mtb\nAM"="+Mtb",
                            "Media\nMDM"="Media", 
                            "+Mtb\nMDM"="+Mtb")) +
  scale_y_continuous(expand = c(0.15,0))

#### Save ####
lo <- "
AAAAACE
BBBBBBD
"

plot_all <- h.plot1+type1+type2+type3+plot_spacer()+
  plot_annotation(tag_levels = "A", tag_suffix = ")") +
  plot_layout(design=lo)
# plot_all

ggsave("Fig2_GSEA.IFN.pdf", plot_all, width = 11, height = 7)

