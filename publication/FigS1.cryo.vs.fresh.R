library(tidyverse)
library(ggpubr)
library(patchwork)

#17 is donor1
#18 is donor2
#21 is donor3

#### B. RNAseq data ####
load("../AM_cyro/data_clean/AM-cryo_voom.RData")

#### Plot ####
p2 <- as.data.frame(dat.voom$E) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "libID") %>% 
  full_join(dat.voom$targets) %>% 
  dplyr::select(rowname:Mtb, -c(libID, group:libID.no)) %>% 
  pivot_wider(names_from = cryo) %>% 
  mutate(Mtb = recode_factor(Mtb, "media"="Media","mtb"="+Mtb")) %>% 
  mutate(donorID2 = case_when(donorID == "AM18"~"Donor 2",
                              donorID == "AM21"~"Donor 3"),
         facet.lab = paste(Mtb, donorID2, sep="\n"),
         facet.lab = factor(facet.lab, 
                            levels=c("Media\nDonor 2",
                                     "Media\nDonor 3",
                                     "+Mtb\nDonor 2",
                                     "+Mtb\nDonor 3"))) %>% 

  ggplot(aes(x=fresh, y=cryo)) +
  geom_point(alpha=0.1) +
  facet_wrap(~facet.lab, nrow = 1) +
  geom_abline(aes(intercept = 0, slope = 1, color="1:1 fit"), 
              size=1, show.legend = FALSE) +
  geom_smooth(method="lm", formula = "y~x", se=FALSE, 
              aes(color="True y~x fit"),
              show.legend = TRUE) +
  scale_color_manual(values=c("grey","red"), name="") +
  theme_classic() +
  coord_fixed() +
  labs(x="Fresh\nNormalized log2 expression",
       y="Cryopreserved\nNormalized log2 expression") +
  stat_cor(method = "pearson", color="red")

# p2

#### Cytokines ####
cyto <- read_csv("../AM_cyro/data_raw/excreted_cytokines_GJPeterson.csv") %>% 
  # None detected group
  group_by(cyto, Media, Stim) %>% 
  mutate(mean_adj_conc = mean(`Adj. Conc.`, na.rm = TRUE),
         ND = ifelse(mean_adj_conc == 0, "N/D",NA)) %>% 
  ungroup() %>% 
  mutate(Mtb = case_when(Stim == "Media"~"Media",
                         Stim == "MOI 2.5"~"+Mtb"),
         Mtb = factor(Mtb, levels=c("Media","+Mtb"))) %>% 
  mutate(group = recode(Media,
                        "10% autologous"="Fresh",
                        "10% FBS"="Cyropreserved\n(10% FBS)",
                        "20% FBS"="Cyropreserved\n(20% FBS)"),
         group = fct_relevel(group, "Fresh", after=0),
         donorID2 = case_when(ID == "AM17"~"Donor 1",
                              ID == "AM18"~"Donor 2",
                              ID == "AM21"~"Donor 3")) %>% 
  mutate(y2 = case_when(cyto=="IL1b"~150,
                        TRUE~500))
 
p3 <- cyto %>% 
  ggplot(aes(x=Mtb, y=`Adj. Conc.`, shape=donorID2, group=group,
             color=group)) +
  geom_point(color = "grey80", position = position_dodge(width=0.75)) +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25,
               position = position_dodge(width = 0.75)) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5,
               position = position_dodge(width = 0.75)) +
  theme_classic(base_size = 12) +
  facet_wrap(~cyto, scales = "free") +
  labs(x="",y="Concentration (pg/ml)", shape="",
       color="") +
  geom_text(aes(label=ND, y=y2),
            position = position_dodge(width = 0.75),
            show.legend = FALSE) +
  scale_color_manual(values=c("#56B4E9","#009E73","#D55E00"))
# p3

#### Save ####

plot_all <- plot_spacer() + p2 + p3 +
  plot_annotation(tag_levels = list(c("B","C")), tag_suffix = ")") +
  plot_layout(ncol=1, heights = c(1,1.2,1))
# plot_all

ggsave("FigS1_cryo.vs.fresh.pdf",plot_all, height=8, width=12)
