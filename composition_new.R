library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(microbiome)

setwd("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/RDS")

type_order <- c(
  "Teat apex",
  "Teat canal",
  "Stripped milk",
  "Cisternal milk"
)

TA.genus <- unique(melt.ps.all.TA.rel.genus$Genus)
TC.genus <- unique(melt.ps.all.TC.rel.genus$Genus)
SM.genus <- unique(melt.ps.all.QM.rel.genus$Genus)
CM.genus <- unique(melt.ps.all.CM.rel.genus$Genus)

genus.order <- c("Acinetobacter",
                "Aerococcus",
                "Atopostipes",
                "Bacteroides",
                "Clostridium sensu stricto 1",
                 "Corynebacterium",
                "Kocuria",
                "Methanobrevibacter",
                "Rothia",
                "Rikenellaceae RC9 gut group",
                "Staphylococcus",
                "Streptococcus",
                "Succinivibrio",
                "Terrisporobacter",
                "Turicibacter",
                "UCG-005",
                "Weissella",
                "Unknown",
                "Other")

scaleColorFillManual <-
  scale_fill_manual(
    values = 
      c(
        Acinetobacter="#9edae5",
        Aerococcus="#17becf",
        Atopostipes="#dbdb8d",
        Bacteroides="#bcbd22",
        "Clostridium sensu stricto 1"="#c7c7c7",
        Corynebacterium="#7f7f7f",
        Kocuria="#f7b6d2",
        Methanobrevibacter="#e377c2",
        Rothia="#c49c94",
        "Rikenellaceae RC9 gut group"="#8c564b",
        Staphylococcus="#c5b0d5",
        Streptococcus="#9467bd",
        Succinivibrio="#ff9896",
        Terrisporobacter="#d62728",
        Turicibacter="#98df8a",
        "UCG-005"="#2ca02c",
        Weissella="#ffbb78",
        Unknown="#ff7f0e",
        Other="#aec7e8"
      ))

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

fig.path <- ("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/Figures")
ggsave(file.path(fig.path, "Figure 7. phylum composition_sample type_checked_grid.png"), p.phylum, dpi=600)


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

fig.path <- ("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/Figures")
ggsave(file.path(fig.path, "Figure 7. genus composition_sample type_checked_grid.png"), p.genus, dpi=600)


##############################################################

## 5.1 teat apex
ps.all.TA <- subset_samples(ps.all, Type %in% "Teat apex")
ps.all.TA

ps.all.TA.rel <- microbiome::transform(ps.all.TA, "compositional")
ps.all.TA.rel.genus <- microbiome::aggregate_rare(ps.all.TA.rel, level = "Genus", detection = 0.03, prevalence = 0.2)

melt.ps.all.TA.rel.genus <- psmelt(ps.all.TA.rel.genus)
melt.ps.all.TA.rel.genus$Dataset <- factor(melt.ps.all.TA.rel.genus$Dataset, levels = c("Before decontam", "After decontam", "After SourceTracker"))
melt.ps.all.TA.rel.genus$Genus <- factor(melt.ps.all.TA.rel.genus$Genus, levels = genus.order)


p.genus.TA <- ggplot(melt.ps.all.TA.rel.genus, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance",
       title = "Teat apex") +
  scaleColorFillManual
p.genus.TA

## 5.2 teat canal
ps.all.TC <- subset_samples(ps.all, Type %in% "Teat canal")
ps.all.TC

ps.all.TC.rel <- microbiome::transform(ps.all.TC, "compositional")
ps.all.TC.rel.genus <- microbiome::aggregate_rare(ps.all.TC.rel, level = "Genus", detection = 0.03, prevalence = 0.2)

melt.ps.all.TC.rel.genus <- psmelt(ps.all.TC.rel.genus)
melt.ps.all.TC.rel.genus$Dataset <- factor(melt.ps.all.TC.rel.genus$Dataset, levels = c("Before decontam", "After decontam", "After SourceTracker"))
melt.ps.all.TC.rel.genus$Genus <- factor(melt.ps.all.TC.rel.genus$Genus, levels = genus.order)


p.genus.TC <- ggplot(melt.ps.all.TC.rel.genus, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance",
       title = "Teat canal") +
  scaleColorFillManual
p.genus.TC

## 5.3 Cisternal milk
ps.all.CM <- subset_samples(ps.all, Type %in% "Cisternal milk")
ps.all.CM

ps.all.CM.rel <- microbiome::transform(ps.all.CM, "compositional")
ps.all.CM.rel.genus <- microbiome::aggregate_rare(ps.all.CM.rel, level = "Genus", detection = 0.03, prevalence = 0.2)

melt.ps.all.CM.rel.genus <- psmelt(ps.all.CM.rel.genus)
melt.ps.all.CM.rel.genus$Dataset <- factor(melt.ps.all.CM.rel.genus$Dataset, levels = c("Before decontam", "After decontam", "After SourceTracker"))
melt.ps.all.CM.rel.genus$Genus <- factor(melt.ps.all.CM.rel.genus$Genus, levels = genus.order)


p.genus.CM <- ggplot(melt.ps.all.CM.rel.genus, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance",
       title = "Cisternal milk") +
  scaleColorFillManual
p.genus.CM

## 5.4 Stripped milk
ps.all.QM <- subset_samples(ps.all, Type %in% "Stripped milk")
ps.all.QM

ps.all.QM.rel <- microbiome::transform(ps.all.QM, "compositional")
ps.all.QM.rel.genus <- microbiome::aggregate_rare(ps.all.QM.rel, level = "Genus", detection = 0.03, prevalence = 0.2)

melt.ps.all.QM.rel.genus <- psmelt(ps.all.QM.rel.genus)
melt.ps.all.QM.rel.genus$Dataset <- factor(melt.ps.all.QM.rel.genus$Dataset, levels = c("Before decontam", "After decontam", "After SourceTracker"))
melt.ps.all.QM.rel.genus$Genus <- factor(melt.ps.all.QM.rel.genus$Genus, levels = genus.order)


p.genus.QM <- ggplot(melt.ps.all.QM.rel.genus, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(y = "Relative abundance",
       title = "Stripped milk") +
  scaleColorFillManual
p.genus.QM

## merge the plot
p.genus.rel <- (p.genus.TA / p.genus.TC / p.genus.CM / p.genus.QM) +
            plot_annotation(tag_levels = "A")
p.genus.rel

fig.path <- ("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/Figures")
ggsave(file.path(fig.path, "Figure 7. genus composition_sample type_checked.png"), p.genus.rel, width=10, height=12, dpi=600)
