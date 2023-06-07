library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(microbiome)


# set file path

path.rds <("Analysis/RDS")
fig.path <- ("Analysis/Figures")

# set sample type and genus order and genus color
type_order <- c(
  "Teat apex",
  "Teat canal",
  "Stripped milk",
  "Cisternal milk"
)


## 1. before decontam

ps <- readRDS("phyloseq.rds")
ps <- subset_samples(ps, Type %in% type_order)

otu <- otu_table(ps)
tax <- tax_table(ps)

meta <- sample_data(ps) %>% data.frame
meta$Dataset <- rep("Before decontam", 55)
meta <- sample_data(meta)

ps <- phyloseq(otu, tax, meta)
ps

## 2. after decontatm

ps.decontam <- readRDS("ps.decontam.rds")
ps.decontam <- subset_samples(ps.decontam, Type %in% type_order)
ps.decontam <- prune_taxa(taxa_sums(ps.decontam)>0, ps.decontam)
ps.decontam

tax.decontatm <- tax_table(ps.decontam)

otu.decontam <- otu_table(ps.decontam) %>% data.frame()
rownames(otu.decontam) <- paste(rownames(otu.decontam), "After decontam")
otu.decontam <- otu_table(otu.decontam, taxa_are_rows = FALSE)

meta.decontam <- sample_data(ps.decontam) %>% data.frame
meta.decontam$Dataset <- rep("After decontam", 55)
rownames(meta.decontam) <- paste(rownames(meta.decontam), "After decontam")
meta.decontam <- sample_data(meta.decontam)

ps.decontam <- phyloseq(otu.decontam, tax.decontatm, meta.decontam)
ps.decontam


## 3. after sourcetracker

ps.sourcetracker <- readRDS("ps.decontam.st.checked.rds")
ps.sourcetracker <- prune_taxa(taxa_sums(ps.sourcetracker)>0, ps.sourcetracker)
ps.sourcetracker

tax.sourcetracker <- tax_table(ps.sourcetracker)

otu.sourcetracker <- otu_table(ps.sourcetracker) %>% data.frame()
rownames(otu.sourcetracker) <- paste(rownames(otu.sourcetracker), "After SourceTracker")
otu.sourcetracker <- otu_table(otu.sourcetracker, taxa_are_rows = FALSE)

meta.sourcetracker <- sample_data(ps.sourcetracker) %>% data.frame
meta.sourcetracker$Dataset <- rep("After SourceTracker", 55)
rownames(meta.sourcetracker) <- paste(rownames(meta.sourcetracker), "After SourceTracker")
meta.sourcetracker <- sample_data(meta.sourcetracker)

ps.sourcetracker <- phyloseq(otu.sourcetracker, tax.sourcetracker, meta.sourcetracker)
ps.sourcetracker


## 4. merge phyloseq project

ps.all <- merge_phyloseq(ps, ps.decontam, ps.sourcetracker)
ps.all

ps.all.rel <- microbiome::transform(ps.all, "compositional")

## 5.1 aggregate on phylum level
ps.all.rel.phylum <- microbiome::aggregate_rare(ps.all.rel, level = "Phylum", detection = 0.01, prevalence = 0.2)

ps.melt.all.phylum <- psmelt(ps.all.rel.phylum)
ps.melt.all.phylum$Dataset <- factor(ps.melt.all.phylum$Dataset, levels=c("Before decontam", "After decontam", "After SourceTracker"))
ps.melt.all.phylum$Type <- factor(ps.melt.all.phylum$Type, levels=c("Teat apex", "Teat canal", "Stripped milk", "Cisternal milk"))

p.phylum <- ggplot(ps.melt.all.phylum, aes(x = X.SampleID, y = Abundance, fill = Phylum)) +
  geom_bar(stat="identity") +
  facet_grid(Dataset ~ Type, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance")
p.phylum

ggsave(file.path(fig.path, "Figure 7. phylum composition.png"), p.phylum, dpi=600)


## 5.2 aggregate on genus level
genus.order.genus <- c("Acinetobacter",
                 "Bacteroides",
                 "Clostridium sensu stricto 1",
                 "Corynebacterium",
                 "Kocuria",
                 "Methanobrevibacter",
                 "Staphylococcus",
                 "UCG-005",
                 "Unknown",
                 "Other")

scaleColorFillManual.genus <-
  scale_fill_manual(
    values = 
      c(
        Acinetobacter="#9edae5",
        Bacteroides="#bcbd22",
        "Clostridium sensu stricto 1"="#c7c7c7",
        Corynebacterium="#7f7f7f",
        Kocuria="#f7b6d2",
        Methanobrevibacter="#e377c2",
        Staphylococcus="#c5b0d5",
        "UCG-005"="#2ca02c",
        Unknown="#ff7f0e",
        Other="#aec7e8"
      ))



ps.all.rel.genus <- microbiome::aggregate_rare(ps.all.rel, level = "Genus", detection = 0.03, prevalence = 0.2)

ps.melt.all.genus <- psmelt(ps.all.rel.genus)
ps.melt.all.genus$Dataset <- factor(ps.melt.all.genus$Dataset, levels=c("Before decontam", "After decontam", "After SourceTracker"))
ps.melt.all.genus$Type <- factor(ps.melt.all.genus$Type, levels=c("Teat apex", "Teat canal", "Stripped milk", "Cisternal milk"))
ps.melt.all.genus$Genus <- factor(ps.melt.all.genus$Genus, levels = genus.order.genus) 

p.genus <- ggplot(ps.melt.all.genus, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_grid(Dataset ~ Type, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance")+
  scaleColorFillManual.genus 
p.genus

ggsave(file.path(fig.path, "Figure 7. genus composition.png"), p.genus, dpi=600)

