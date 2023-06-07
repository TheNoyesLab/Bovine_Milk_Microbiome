library(ggplot2)
library(dplyr)
library(phyloseq)
library(microbiome)
library(patchwork)

setwd("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis")

path.rds <- ("RDS")
path.figures <- ("Figures/")

# mastitis pathogens
pathogens <- c("Streptococcus", "Staphylococcus", "Mycoplasma", "Prototheca", "Aerococcus", "Enterococcus", "Lactococcus", "Micrococcus", "Trueperella", "Bacillus", "Corynebacterium", "Listeria", "Acinetobacter", "Aeromonas", "Citrobacter", "Enterobacter", "Escherichia/Shigella", "Flavimonas", "Hafnia", "Klebsiella", "Pantoea", "Plesimonas", "Proteus", "Pseudomonas", "Salmonella", "Serrati", "Serratia", "Stenotrophomonas", "Yersinia", "Yeast", "Nocardia")

ps.sourcetracker <- readRDS(file.path(path.rds, "ps.decontam.st.checked.rds"))


# change ps in relative abundance
ps.rel <- microbiome::transform(ps.sourcetracker, "compositional")
ps.mast.rel <- subset_taxa(ps.rel, Genus %in% pathogens)
ps.mast.gen.rel <-aggregate_rare(ps.mast.rel, level = "Genus", detection = 0.001, prevalence = 0.25)

# plot relative abundance of mastitis per sample
melt.mast <- psmelt(ps.mast.gen.rel)
melt.mast<- melt.mast %>% group_by(Type, Genus) 
p.melt.mast <- ggplot(melt.mast, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(10, "pt")
  )+
  labs(y = "Relative abundance") +
  scale_fill_brewer("Genus", palette = "Paired")

print(p.melt.mast)


ggsave (file.path(path.figures, "Figure 8. Potential mastitis abundance per sample_checked.png"), p.melt.mast, dpi = 600)




# plot grouped relative abundance of mastitis
p.mast.grouped <- plot_composition(ps.mast.gen.rel, sample.sort = NULL,
                                   average_by = "Type",
                                   otu.sort = NULL,
                                   plot.type = "barplot",
                                   verbose = FALSE) +
  xlab("") +
  ylab("Relative abundance")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45),
        axis.ticks.x=element_blank())+
  scale_fill_brewer("Genus", palette = "Paired")
print(p.mast.grouped)

# check Staphylococcus positive samples from CM, get the seq for BLAST
ps.sty.positive <- subset_samples(ps.sourcetracker, X.SampleID %in% c("MS56_S48", "MS35_S3"))
ps.sty.positive <- prune_taxa(taxa_sums(ps.sty.positive)>0, ps.sty.positive)
taxa.sty.positive <- tax_table(ps.sty.positive) %>% data.frame() %>% subset(Genus == ("Staphylococcus"))
seq.sty.positive <- rownames(taxa.sty.positive)

taxa <- readRDS(file.path(path.rds, "taxa.rds")) %>% data.frame()
taxa$seq <- paste0("Seq", seq(dim(taxa)[1]))
taxa.sty.seq <- subset(taxa, taxa$seq %in% seq.sty.positive) %>% rownames()
write.csv (taxa.sty.seq, "results/sequence.staphylococcus.culture.positive.csv")

# get the mean of abundance by Genus
mean.mast.rel <- melt.mast %>% 
  group_by(Genus) %>% 
  summarise (genus_mean = mean(Abundance)*100) %>% 
  arrange(desc(genus_mean))
mean.mast.rel

# get the mean of abundance by Type
mean.mast.rel.type <- melt.mast %>% 
  group_by(Type) %>% 
  summarise (type_mean = mean(Abundance)*100) %>% 
  arrange(desc(type_mean))
mean.mast.rel.type

# get the mean of abundance by Genus and Type
mean.mast.rel.grouped <- melt.mast %>% 
  group_by(Genus,Type) %>% 
  summarise (genus_mean = mean(Abundance)*100)
mean.mast.rel.grouped$genus_mean <- round(mean.mast.rel.grouped$genus_mean, digits=1)
mean.mast.rel.grouped

# merge the plots
p.mast <- (p.melt.mast | p.mast.grouped) + plot_layout(widths = c(3, 1)) + plot_annotation(tag_levels = c("A"))
p.mast
ggsave (file.path(path.figures, "Figure 8. Potential mastitis abundance per sample.png"), p.melt.mast, dpi = 600)
