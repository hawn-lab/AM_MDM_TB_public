library(tidyverse)

align.lib <- read_tsv("results/MTB.identity.check/summary.alignment.tsv", 
                      col_names = FALSE,
                      col_types = cols()) %>% 
  filter(grepl("bam",X1)) %>% 
  distinct(X1)

# Align summary
align.summ <- read_tsv("results/MTB.identity.check/summary.alignment.tsv", 
                       col_names = FALSE,
                       col_types = cols()) %>% 
  #Set groups
  mutate(group = ifelse(grepl("_Aligned.sortedByCoord.out.bam", X1),"align",
                        ifelse(grepl("_filter_paired.bam", X1),"filter",
                              NA))) %>% 
  #Get sample name from filename
  mutate(libID = factor(X1, levels=align.lib$X1),
         libID = basename(as.character(libID)),
         libID = gsub("_Aligned.sortedByCoord.out.bam|_filter_paired.bam",
                      "", libID)) %>% 
  fill(libID,group) %>% 
  #Separate data to column
  separate(X1, into=c("h","i"), sep=" \\+ 0 ", fill="right") %>% 
  drop_na(i) %>% 
  #Recode data types (f)
  separate(i, into=c("i"), sep="[(]", extra="drop") %>% 
  mutate(i = fct_recode(factor(i),
                        to.be.aligned="in total ",
                        secondary.align="secondary",
                        chimeric.align="supplementary",
                        PCR.dups="duplicates",
                        align="mapped ",
                        paired="paired in sequencing",
                        R1.paired="read1", R2.paired="read2",
                        align.paired="properly paired ",
                        both.align.paired= "with itself and mate mapped",
                        one.align.paired="singletons " ,
                        both.align.paired.diffCHR="with mate mapped to a different chr",
                        both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>% 
  pivot_wider(names_from = "i", values_from = "h") %>% 
  mutate_at(vars(-libID, -group), as.numeric) %>% 
  select(group, libID, both.align.paired)

#Combine with raw info
metric.all <- read_csv("data_raw/cleaning.metrics.csv") %>% 
  select(libID, to.be.aligned) %>% 
  right_join(align.summ) %>% 
  filter(group=="filter") %>% 
  mutate(perc = both.align.paired/to.be.aligned*100) %>% 
  mutate(mtb= ifelse(grepl("MTB",libID),"Mtb","Media")) %>% 
  separate(libID, into=c("libNO", "ptID", "cryo"), sep="_", remove = FALSE) %>% 
  mutate(cryo = ifelse(cryo != "Cryo","fresh", "cryo"))

#### plot ####
plot <- metric.all %>% 
  
  ggplot(aes(x=paste(ptID,cryo,sep="\n"), y=perc, fill=mtb)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  labs(fill="", x="", y="Percent PF reads aligned to Mtb genome")

ggsave("results/MTB.identity.check/MTB.identity.check.png", plot,
       height=5, width=4)  
  