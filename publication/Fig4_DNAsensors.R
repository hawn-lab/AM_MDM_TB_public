library(tidyverse)
library(limma)
library(rstatix)
library(ggpubr)
library(patchwork)

#### Data ####
load("../AM_MDM_TB/data_clean/AM-MDM_voom.RData")
load("../AM_MDM_TB/results/gene_level/kimma_lme_cell_TB.RData")
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

#### Fig 5: DNA sensors and IRF ####
#genes of interest
genes1 <- c("CGAS","DDX60", "DHX36", "IFI16", "ZBP1", "DHX9", "MRE11", "PYHIN1", "STING1") %>% 
  sort()

dat.subset <- dat.format %>% 
  filter(symbol %in% genes1) %>% 
  #Order and add tags
  mutate(symbol = factor(symbol, levels=genes1)) 
fdr.subset <- fdr.format %>% filter(symbol %in% genes1)

#Format fdr for plotting
stat.dat <- fdr.subset %>%  
  #add y position
  left_join(dat.summ) %>% 
  group_by(symbol) %>% 
  mutate(rnk = order(FDR, decreasing=FALSE)) %>% 
  ungroup() %>% 
  mutate(y.position = maxE+abs(maxE*0.01*rnk)) %>%
  #Format FDR labels
  add_significance("FDR", cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("**", "*", "ns"))

p5 <- dat.subset %>% 
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
        axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        strip.background = element_rect(linewidth=0.5)) +
  scale_color_manual(values = col.vec) +       
  labs(y="Normalized log2 expression", x="", color="",
       title="DNA sensors") +
  facet_wrap(~symbol, scales="free", nrow=2) +
  #Significant bars
  stat_pvalue_manual(stat.dat, 
                     label = "FDR.signif", tip.length = 0.02) +
  # clean x axis labels
  scale_x_discrete(labels=c("Media\nAM"="Media",
                            "+Mtb\nAM"="+Mtb",
                            "Media\nMDM"="Media", 
                            "+Mtb\nMDM"="+Mtb")) +
  scale_y_continuous(expand = c(0.15,0))
p5

#### Save ####
ggsave("Fig4_sensor.genes.pdf", p5, height=6, width=6)

