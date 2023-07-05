library(tidyverse)
library(limma)
library(ggvenn)
library(patchwork)
library(DiagrammeR)

#### Flow diagram ####
L <- 0.7
R <- 1.3
C <- 1

flowchart <- data.frame(
  label=c(
    'BAL',
    'Pen/Strep 2h\nWashings x6',
    'Rest\novernight',
    'Freeze PBMC',
    'Thaw cells\nin batch',
    'M-CSF x 5 days\nCD14 negative isolation',
    'Mock vs. H37Rv (MOl 1:1) 6h',
    'TRIzol lysates',
    'RNA isolation',
    'RNA sequencing'),
  x = c(L,L,L, R,R,R, C,C,C,C),
  y = c(7,6,5, 7,6,5, 4,3,2,1)) %>% 
  ggplot() +
  geom_text(aes(x=x,y=y,label=label)) +
  #left
  geom_segment(x=L, xend=L, y=6.75, yend=6.45,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=L, xend=L, y=5.65, yend=5.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=L, xend=0.9, y=4.65, yend=4.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  #right
  geom_segment(x=R, xend=R, y=6.75, yend=6.45,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=R, xend=R, y=5.65, yend=5.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=R, xend=1.1, y=4.65, yend=4.35,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  #middle
  geom_segment(x=C, xend=C, y=3.75, yend=3.25,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=C, xend=C, y=2.75, yend=2.25,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +
  geom_segment(x=C, xend=C, y=1.75, yend=1.25,
               arrow=arrow(length = unit(0.1, "inches"),
                           ends = "last", type = "open")) +

  theme_void() +
  lims(x=c(0.5,1.5))

#### Data ####
load("../AM_MDM_TB/data_clean/AM-MDM_voom.RData")
load("../AM_MDM_TB/results/gene_level/kimma_lme_cell_TB.RData")
col.vec <- c("#117733","#44AA99","#A46199","#882255")

#### PCA ####
# Calculate PCA
PCA <- as.data.frame(dat.voom$E) %>% 
  t() %>% 
  prcomp(scale. = TRUE, center = TRUE)

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1]*100,digits=1), "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2]*100,digits=1), "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.voom$targets, by="libID") %>% 
  mutate(TB = recode(TB, "TB"="+Mtb","Media"="Media"),
         col.group = paste(cell,TB),
         col.group = factor(col.group, levels=c("AM Media","AM +Mtb",
                                                "MDM Media","MDM +Mtb")))

PCA1 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_point(aes(color=col.group, shape = col.group),size=3) +
  #Beautify
  theme_classic(base_size = 16) +
  labs(x=PC1.label, y=PC2.label, color="") +
  coord_fixed(ratio=1) +
  theme(panel.border = element_rect(colour = "black", 
                                    fill=NA, linewidth=1),
        # legend.position = "bottom"
        ) +
  scale_color_manual(values = col.vec, name="") +
  scale_shape_manual(values = c(17,17,19,19), name="") #+
  # guides(color = guide_legend(ncol=2),
  #        shape = guide_legend(ncol=2))

# PCA1

#### Venn ####
#List genes significant for cell within a TB condition
signif.wIn.condition <- AM.MDM_kimma_pval$lme.contrast %>% 
  filter((grepl("Media", contrast_ref) & grepl("Media", contrast_lvl)) |
           (grepl("TB", contrast_ref) & grepl("TB", contrast_lvl))) %>% 
  filter(FDR < 0.05) %>% 
  pull(gene) %>% unique()

#Separate up vs down with Mtb within each cell type
for(cell in c("AM","MDM")){
  for(direction in c("up","down")){
    name <- paste("kimma",cell, direction, sep=".")
    
    if(direction == "up"){
      genes <- AM.MDM_kimma_pval$lme.contrast %>% 
        filter(grepl(cell, contrast_ref) & grepl(cell, contrast_lvl) &
                 FDR < 0.05 & estimate > 0 &
                 gene %in% signif.wIn.condition) %>% 
        pull(gene)
    }else{
      genes <- AM.MDM_kimma_pval$lme.contrast %>% 
        filter(grepl(cell, contrast_ref) & grepl(cell, contrast_lvl) &
                 FDR < 0.05 & estimate < 0 &
                 gene %in% signif.wIn.condition) %>% 
        pull(gene)
    }
    
    assign(name, genes, envir = .GlobalEnv)
  }}

venn_dat <- list("AM\nHigher in +Mtb"=kimma.AM.up,
                 "AM\nLower in +Mtb                "=kimma.AM.down,
                 "MDM\n            Higher in +Mtb"=kimma.MDM.up,
                 "MDM\nLower in +Mtb"=kimma.MDM.down)

venn_plot <- ggvenn(venn_dat, show_percentage = FALSE, text_size = 4, set_name_size = 4, 
                    stroke_size = 0.5) + 
  scale_fill_manual(values = col.vec[c(2,2,3,3)]) +
  lims(x=c(-2.5,2.5))
# venn_plot

#### Save ####
lo <- "
ABB
ACC
ACC
"
plot_all <- flowchart + PCA1 + venn_plot + 
  plot_layout(design = lo) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 14))

# plot_all

ggsave("Fig1_flow.PCA.venn.pdf", plot_all, width = 9, height = 6)
