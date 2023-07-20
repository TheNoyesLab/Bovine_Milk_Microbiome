library(ggplot2)
library(dplyr)
library(phyloseq)
library(microbiome)
library(patchwork)

# set file path

path.rds <- ("Analysis/RDS")
path.figures <- ("Analysis/Figures")

# mastitis pathogens
pathogens <- c("Streptococcus", "Staphylococcus", "Mycoplasma", "Prototheca", "Aerococcus", "Enterococcus", "Lactococcus", "Micrococcus", "Trueperella", "Bacillus", "Corynebacterium", "Listeria", "Acinetobacter", "Aeromonas", "Citrobacter", "Enterobacter", "Escherichia/Shigella", "Flavimonas", "Hafnia", "Klebsiella", "Pantoea", "Plesimonas", "Proteus", "Pseudomonas", "Salmonella", "Serrati", "Serratia", "Stenotrophomonas", "Yersinia", "Yeast", "Nocardia")

ps.sourcetracker <- readRDS(file.path(path.rds, "ps.decontam.st.rds"))


# change ps in relative abundance
ps.rel <- microbiome::transform(ps.sourcetracker, "compositional")
ps.mast.rel <- subset_taxa(ps.rel, Genus %in% pathogens)
ps.mast.gen.rel <-aggregate_rare(ps.mast.rel, level = "Genus", detection = 0.001, prevalence = 0.25)

# plot relative abundance of mastitis per sample
melt.mast <- psmelt(ps.mast.gen.rel)

melt.mast$CowId <- as.factor(melt.mast$CowId)

p.melt.mast <- ggplot(melt.mast, aes(x = CowId, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.5))+
  labs(y = "Relative abundance") +
  scale_fill_brewer("Genus", palette = "Paired")

print(p.melt.mast)

ggsave (file.path(path.figures, "Figure 8. Potential mastitis abundance.png"), p.melt.mast, dpi = 1200, width = 12, height = 6)
