library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(ggpubr)
library(ggvenn)
library(patchwork)
library(DiagrammeR)

col.vec1 <- c("#DDCC77","#117733","#56B4E9","#D55E00","#BD4B8A")
# col.vec <- c("#648FFF","#FFB000","#DC267F")

#### DNA sensors and LPS ####
dat <- read_excel(
  "../addtl_experiments//DNA_sensor_ligand_paired_AM_MDM_06-21-2023.xlsx",
  sheet="combined") %>% 
  mutate(log10Delta = log10(`2DeltaCT_norm`)) %>% 
  filter(stim %in% c("Media","LPS") & lipo == "-") %>% 
  mutate(stim = fct_relevel(stim, "Media", after = 0))

#### Stats ####
stat_result_int <- data.frame()

for(cy in unique(dat$cyto)){
  print(cy)
  dat_model <- dat %>% filter(cyto == cy)
  
  m1 <- lmer(log10Delta ~ cell*stim+(1|ptID), data=dat_model)
  
  m1.contrast <- emmeans(m1, adjust = "none", 
                         pairwise~cell:stim)$contrasts %>% 
    broom::tidy() %>% 
    separate(contrast, into = c("contrast_ref", 
                                "contrast_lvl"), sep = " - ")
  
  stat_result_int <- m1.contrast %>% 
    mutate(cyto=cy) %>% 
    bind_rows(stat_result_int)
}

contrast.OI <-
  data.frame(contrast_ref = c("AM Media",
                              "AM Media", "MDM Media", "AM LPS"),
             contrast_lvl = c("MDM Media",
                              "AM LPS", "MDM LPS", "MDM LPS"))

dat_summ <- dat %>% 
  group_by(cyto) %>% 
  summarise(maxDelta=max(log10Delta, na.rm=TRUE),
            minDelta = min(log10Delta, na.rm=TRUE),
            .groups = "drop") %>% 
  ungroup()

int_signif <- stat_result_int %>% 
  inner_join(contrast.OI) %>% 
  group_by(contrast_ref, contrast_lvl) %>% 
  mutate(FDR = p.adjust(p.value)) %>% 
  ungroup() %>% 
  filter(FDR<0.05) %>%
  mutate(signif = case_when(FDR < 0.01~"**",
                            FDR < 0.05~"*",
                            TRUE~"")) %>% 
  rename(group1=contrast_ref, group2=contrast_lvl) %>% 
  select(cyto, group1, group2, FDR, p.value, FDR, signif) %>% 
  #add y position
  left_join(dat_summ) %>% 
  group_by(cyto) %>% 
  mutate(rnk = order(FDR, decreasing=FALSE)) %>% 
  ungroup() %>% 
  mutate(y.position = maxDelta+abs(maxDelta*0.3*rnk)) %>% 
  distinct(group1, group2, y.position, p.value, FDR, signif, cyto) 

#### Plot LPS ####
p1 <- dat %>% 
  mutate(cont = paste(cell, stim),
         cont = factor(cont, levels=c("AM Media","AM LPS","MDM Media","MDM LPS"))) %>% 
  
  ggplot(aes(x=cont, y=log10Delta, color=stim)) +
  #Raw data points
  geom_point() +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5) +
  theme_classic(base_size = 12) +
  facet_wrap(~cyto, nrow=1) +
  labs(x="",y="Log10 normalized expression", shape="",
       color="", title="LPS") +
  scale_color_manual("",values=c("grey80", "grey20")) +
  #Significant bars
  stat_pvalue_manual(int_signif,
                     label = "signif", tip.length = 0.02) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill=NA, 
                                    linewidth=1),
        legend.position = "none")
# p1

#### Cytokine data ####
cyto <- read_csv("../addtl_experiments/cyto_elisa.csv") %>% 
  pivot_longer(Mtb_1:Media_3) %>% 
  separate(name, into = c("TB","rep"))

#### Stats ####
stat_result <- data.frame()

for(cy in unique(cyto$cyto)){
  temp <- cyto %>% filter(conc %in% c(0,10) & cyto == cy)
  temp2 <- lm(value ~ stim*TB, temp) %>% broom::tidy()
  stat_result <- temp2 %>% 
    mutate(cyto=cy) %>% 
    bind_rows(stat_result)
}

stat_plot <-stat_result %>% 
  filter(p.value<0.05 & grepl(":",term)) %>% 
  separate(term, into=c("stim"), sep=":", extra = "drop") %>% 
  mutate(stim = gsub("stim","",stim)) %>% 
  mutate(signif = case_when(p.value < 0.01~"**",
                            p.value < 0.05~"*")) %>% 
  left_join(cyto %>% filter(conc%in%c(0,10))) %>% 
  group_by(cyto,stim,p.value,signif,TB) %>% 
  summarise(value = max(value)+200) %>% 
  filter(TB=="Mtb") %>% 
  mutate(x=paste(stim,TB))


p3 <- cyto %>% 
  filter(conc%in%c(0,10)) %>% 
  ggplot(aes(x=stim, y=value, color=TB, fill=TB, group=TB)) +
  #Raw data points
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.75)) +
  #Stdev bar
  stat_summary(fun.data=mean_sdl,
               fun.args = list(mult=1),
               geom="errorbar", width=0.25,
               position = position_dodge(width = 0.75)) +
  theme_classic(base_size = 12) +
  facet_wrap(~cyto, scales = "free") +
  labs(x="",y="Concentration (pg/ml)", shape="",
       color="", title="IFN") +
  scale_color_manual("",values=c("grey80", "grey20"), labels=c("Media","+Mtb")) +
  scale_fill_manual("",values=c("grey80", "grey20"), labels=c("Media","+Mtb")) +
  #signif markers
  geom_text(data=stat_plot, aes(label=signif), color="black",size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, 
                                    linewidth=1))
# p3

#### Flow diagram ####
L <- 0.7
R <- 1.3
C <- 1

p2 <- data.frame(
  label=c(
    'Freeze PBMC',
    'Thaw cells\nin batch',
    'M-CSF x 5 days\nCD14 negative isolation',
    
    'Mock vs. IFN 24 hrs',
    'Mock vs. H37Rv\n(MOl 1:1) 24hr',
    'Supernatants',
    'C) ELISA',
    
    'Mock vs. IFN 6hr',
    'TRIzol\nlysates',
    'RNA isolation',
    'D) RNA sequencing'),
  x = c(C,C,C, L,L,L,L, R,R,R,R),
  y = c(7,6,5, 4,3,2,1, 4,3,2,1)) %>% 
  ggplot() +
  geom_text(aes(x=x,y=y,label=label)) +
  #middle
  geom_segment(x=C, xend=C, y=6.75, yend=6.25,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=C, xend=C, y=5.75, yend=5.25,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  #left
  geom_segment(x=L, xend=L, y=3.75, yend=3.45,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=L, xend=L, y=2.65, yend=2.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=L, xend=L, y=1.65, yend=1.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=0.9, xend=L, y=4.65, yend=4.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  #right
  geom_segment(x=R, xend=R, y=3.75, yend=3.45,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=R, xend=R, y=2.65, yend=2.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=R, xend=R, y=1.65, yend=1.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=1.1, xend=R, y=4.65, yend=4.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  
  theme_void() +
  lims(x=c(0.5,1.5))

# p2

#### RNAseq data ####
load("../MDM_IFN/data_clean/MDM-IFN_voom.RData")
load("../MDM_IFN/results/gene_level/kimma_IFN.RData")
MDM.IFN_kimma_pval$lme.contrast <- MDM.IFN_kimma_pval$lme.contrast %>% 
  mutate(across(c(contrast_ref,contrast_lvl), 
                ~recode(.,"IFNb"="IFNB", "IFNe"="IFNE", "IFNg"="IFNG")))

#### Venn: all IFN ####
venn.ls <- list()
for(ifn in unique(MDM.IFN_kimma_pval$lme.contrast$contrast_lvl)){
  venn.ls[[paste0("MDM +",ifn)]] <-  MDM.IFN_kimma_pval$lme.contrast %>% 
    filter(contrast_ref == "untreated" & contrast_lvl == ifn) %>% 
    filter(FDR < 0.05) %>% 
    select(contrast_lvl, gene, estimate, FDR) %>% 
    pull(gene) %>% unique()
}

venn1 <- ggvenn(venn.ls, show_percentage = FALSE, 
                text_size = 4, set_name_size = 4, 
                stroke_size = 0.5) + 
  scale_fill_manual(values = col.vec1[1:4]) +
  lims(x=c(-2.8,2.8))
# venn1
# IFNL omitted. only has 3 DEG

#### Save ####
lo <- "
AA
BD
CC
"
plot_all <- p1 + p2 + p3 + venn1 +
  plot_layout(heights = c(1,2,1), design = lo)+
  plot_annotation(tag_levels = "A", tag_suffix = ")")
# plot_all

# ggsave("Fig3_IFN.stim.png", plot_all, width=8, height=10)
ggsave("Fig3_IFN.stim.pdf", plot_all, width=8, height=10)
