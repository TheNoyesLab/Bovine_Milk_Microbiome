library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(vegan)
library(pairwiseAdonis)
library(metagenomeSeq)

# set file paths
path.rds <- ("Analysis/RDS")
path.figure <- ("Analysis/Figures")

SAMPLE_TYPES <- c(
  "Teat apex",
  "Teat canal",
  "Stripped milk",
  "Cisternal milk",
  "Positive",
  "Air",
  "Blank",
  "Extraction",
  "Library"
)

scaleColorFillManualPanel <-
  scale_fill_manual(
    values = 
      c(
        "Teat apex" = "#1f77b4",
        "Teat canal" = "#ff7f0e",
        "Stripped milk" = "#2ca02c",
        "Cisternal milk" = "#d62728",
        "Library" = "#7f7f7f",
        "Extraction" = "#AB1866",
        "Positive" = "#c63d40",
        "Air" = "#539caf",
        "Blank" = "#ebc850"
      )
  )


## 1. pre-decontam

ps <- readRDS(file.path(path.rds, "phyloseq.rds"))
ps 

metadata <- sample_data(ps) %>% data.frame()
metadata$Type <- factor(metadata$Type, levels = SAMPLE_TYPES)

metadata <- metadata %>% mutate(SampleOrControl = case_when(
  Type %in% c("Teat canal", "Teat apex", "Cisternal milk", "Stripped milk") ~ "sample",
  Type %in% c("Library", "Blank", "Air", "Extraction") ~ "control"
))

sample_data(ps) <- sample_data(metadata)
sample_data(ps)$SampleOrControl <- factor(sample_data(ps)$SampleOrControl)

# 1.1 all samples

# PERMANOVA
bray.dist<-vegdist(otu_table(ps), method='bray') 

beta_div <- adonis2(bray.dist ~ Type, data=metadata, permutations=999)
beta_div

beta_div.pair <-pairwise.adonis2(bray.dist ~ Type, data=metadata, permutations=999)
beta_div.pair

anova(betadisper(bray.dist, metadata$Type))

# NMDS ordination
ord.bray.nmds <- ordinate(ps, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds

ord.pnmds <- plot_ordination(ps, ord.bray.nmds, color = "Type", shape = "SampleOrControl") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("Before decontam")+
  annotate("text", x = -2, y = 1, label = "Stress = 0.125\nPERMANOVA P = 0.001\nVariance P < 0.001",hjust = 0)+
  scale_shape_manual(values = c(17, 16))
  
print(ord.pnmds)


## split milk samples and teat samples

# 1.2 subset milk samples

ps.milk <- subset_samples(ps, Type %in% c("Cisternal milk", "Stripped milk"))
ps.milk <- prune_taxa(taxa_sums(ps.milk)>0, ps.milk)
ps.milk

metadata.milk <- sample_data(ps.milk) %>% data.frame()

# PERMANOVA
bray.dist.milk<-vegdist(otu_table(ps.milk), method='bray') 

beta_div.milk<-pairwise.adonis2(bray.dist.milk ~ Type, data=metadata.milk, permutations=999)
beta_div.milk

anova(betadisper(bray.dist.milk, metadata.milk$Type))

# ordination
ord.bray.nmds.milk <- ordinate(ps.milk, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.milk

ord.pnmds.milk <- plot_ordination(ps.milk, ord.bray.nmds.milk, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("Before decontam (milk microbiome)")+
  annotate("text", x = -0.5, y = 0.13, label = "Stress = 0.101\nPERMANOVA P = 0.193\nVariance P = 0.219", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.milk)

# 1.3 subset teat samples
## sample MS11_S44 is an outlier for ordination

ps.teat <- subset_samples(ps, Type %in% c("Teat canal", "Teat apex") &!X.SampleID == "MS11_S44")
ps.teat <- prune_taxa(taxa_sums(ps.teat)>0, ps.teat)
ps.teat

metadata.teat <- sample_data(ps.teat) %>% data.frame()

# PERMANOVA
bray.dist.teat<-vegdist(otu_table(ps.teat), method='bray') 

beta_div.teat<-pairwise.adonis2(bray.dist.teat ~ Type, data=metadata.teat, permutations=999)
beta_div.teat

anova(betadisper(bray.dist.teat, metadata.teat$Type))

# ordination
ord.bray.nmds.teat <- ordinate(ps.teat, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.teat

ord.pnmds.teat <- plot_ordination(ps.teat, ord.bray.nmds.teat, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("Before decontam (udder microbiome)")+
  annotate("text", x = -1, y = 1.1, label = "Stress = 0.134\nPERMANOVA P = 0.273\nVariance P = 0.518", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.teat)


## 2. After decontam

ps.decontam <- readRDS("ps.decontam.rds")
ps.decontam <- subset_samples(ps.decontam, !Type %in% c("Library", "Extraction"))
ps.decontam <- prune_taxa(taxa_sums(ps.decontam)>0, ps.decontam)
ps.decontam

metadata.decontam <- sample_data(ps.decontam) %>% data.frame()
metadata.decontam$Type <- factor(metadata.decontam$Type, levels = SAMPLE_TYPES)

metadata.decontam <- metadata.decontam %>% mutate(SampleOrControl = case_when(
  Type %in% c("Teat canal", "Teat apex", "Cisternal milk", "Stripped milk") ~ "sample",
  Type %in% c("Library", "Blank", "Air", "Extraction") ~ "control"
))

sample_data(ps.decontam) <- sample_data(metadata.decontam)
sample_data(ps.decontam)$SampleOrControl <- factor(sample_data(ps.decontam)$SampleOrControl)

# 2.1 all samples

# PERMANOVA
bray.dist.decontam<-vegdist(otu_table(ps.decontam), method='bray') 

adonis2(bray.dist.decontam ~ Type, data=metadata.decontam, permutations=999)

beta_div.decontam <-pairwise.adonis2(bray.dist.decontam ~ Type, data=metadata.decontam, permutations=999)
beta_div_ps.decontam <- data.frame(beta_div.decontam)
write.csv(beta_div_ps.decontam, file.path(path.results, "pairwise permanova_decontam.csv"))

anova(betadisper(bray.dist.decontam, metadata.decontam$Type))

# ordination
ord.bray.nmds.decontam <- ordinate(ps.decontam, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.decontam

ord.pnmds.decontam <- plot_ordination (ps.decontam, ord.bray.nmds.decontam, 
                                      color = "Type", shape=("SampleOrControl")) +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("After decontam")+
  annotate("text", x = -1, y = 1.3, label = "Stress = 0.160\nPERMANOVA P = 0.001\nVariance P < 0.001", hjust=0)+
  scale_shape_manual(values = c(17, 16))
print(ord.pnmds.decontam)

## split milk samples and teat samples

# 2.2 subset milk samples

ps.milk.decontam <- subset_samples(ps.decontam, Type %in% c("Cisternal milk", "Stripped milk"))
ps.milk.decontam <- prune_taxa(taxa_sums(ps.milk.decontam)>0, ps.milk.decontam)
ps.milk.decontam

metadata.milk.decontam <- sample_data(ps.milk.decontam) %>% data.frame()

# PERMANOVA
bray.dist.milk.decontam<-vegdist(otu_table(ps.milk.decontam), method='bray') 

beta_div.milk.decontam<-pairwise.adonis2(bray.dist.milk.decontam ~ Type, data=metadata.milk.decontam, permutations=999)
beta_div.milk.decontam

anova(betadisper(bray.dist.milk.decontam, metadata.milk.decontam$Type))

# ordination
ord.bray.nmds.milk.decontam <- ordinate(ps.milk.decontam, "NMDS", "bray", k=3, trymax = 1000)
ord.bray.nmds.milk.decontam

ord.pnmds.milk.decontam <- plot_ordination(ps.milk.decontam, ord.bray.nmds.milk.decontam, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("After decontam (milk microbiome)")+
  annotate("text", x = -0.6, y = 0.2, label = "Stress = 0.157\nPERMANOVA P = 0.287\nVariance P = 0.84", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.milk.decontam)

# 2.3 subset teat samples
## sample MS11_S44 is an outlier for ordination

ps.teat.decontam <- subset_samples(ps.decontam, Type %in% c("Teat canal", "Teat apex") &!X.SampleID == "MS11_S44")
ps.teat.decontam <- prune_taxa(taxa_sums(ps.teat.decontam)>0, ps.teat.decontam)
ps.teat.decontam

metadata.teat.decontam <- sample_data(ps.teat.decontam) %>% data.frame()

# PERMANOVA
bray.dist.teat.decontam<-vegdist(otu_table(ps.teat.decontam), method='bray') 

beta_div.teat.decontam<-pairwise.adonis2(bray.dist.teat.decontam ~ Type, data=metadata.teat.decontam, permutations=999)
beta_div.teat.decontam

anova(betadisper(bray.dist.teat.decontam, metadata.teat.decontam$Type))

# ordination
ord.bray.nmds.teat.decontam <- ordinate(ps.teat.decontam, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.teat.decontam

ord.pnmds.teat.decontam <- plot_ordination(ps.teat.decontam, ord.bray.nmds.teat.decontam, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("After decontam (udder microbiome)")+
  annotate("text", x = -1, y = 1, label = "Stress = 0.133\nPERMANOVA P = 0.461\nVariance P = 0.465", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.teat.decontam)


## 3. After sourcetracker

ps.sourcetracker <- readRDS("ps.decontam.st.checked.rds")
ps.sourcetracker <- prune_taxa(taxa_sums(ps.sourcetracker)>0, ps.sourcetracker)
ps.sourcetracker


metadata.sourcetracker <- sample_data(ps.sourcetracker) %>% data.frame()
metadata.sourcetracker$Type <- factor(metadata.sourcetracker$Type, levels = SAMPLE_TYPES)

sample_data(ps.sourcetracker) <- sample_data(metadata.sourcetracker)

# 3.1 NMDS
ord.bray.nmds.sourcetracker <- ordinate(ps.sourcetracker, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.sourcetracker

ord.pnmds.sourcetracker <- plot_ordination(ps.sourcetracker, ord.bray.nmds.sourcetracker, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
          ggtitle("After SourceTracker")+
  annotate("text", x = -0.1, y = 0.15, label = "Stress = 0.161\nPERMANOVA P = 0.001\nVariance P < 0.001 ",
           hjust = 0)

print(ord.pnmds.sourcetracker)

# PERMANOVA
bray.dist.sourcetracker<-vegdist(otu_table(ps.sourcetracker), method='bray') 

adonis2(bray.dist.sourcetracker ~ Type, data=metadata.sourcetracker, permutations=999)

beta_div.sourcetracker <-pairwise.adonis2(bray.dist.sourcetracker ~ Type, data=metadata.sourcetracker, permutations=999)
beta_div_ps.sourcetracker <- data.frame(beta_div.sourcetracker)
write.csv(beta_div_ps.sourcetracker, file.path(path.results, "pairwise permanova_sourcetracker.checked.csv"))

anova(betadisper(bray.dist.sourcetracker, metadata.sourcetracker$Type))

# 3.2 subset milk samples

ps.milk.sourcetracker <- subset_samples(ps.sourcetracker, Type %in% c("Cisternal milk", "Stripped milk"))
ps.milk.sourcetracker <- prune_taxa(taxa_sums(ps.milk.sourcetracker)>0, ps.milk.sourcetracker)
ps.milk.sourcetracker

metadata.milk.sourcetracker <- sample_data(ps.milk.sourcetracker) %>% data.frame()

# PERMANOVA
bray.dist.milk.sourcetracker<-vegdist(otu_table(ps.milk.sourcetracker), method='bray') 

beta_div.milk.sourcetracker<-pairwise.adonis2(bray.dist.milk.sourcetracker ~ Type, data=metadata.milk.sourcetracker, permutations=999)
beta_div.milk.sourcetracker

anova(betadisper(bray.dist.milk.sourcetracker, metadata.milk.sourcetracker$Type))

# ordination
ord.bray.nmds.milk.sourcetracker <- ordinate(ps.milk.sourcetracker, "NMDS", "bray", k=3, trymax = 1000)
ord.bray.nmds.milk.sourcetracker

ord.pnmds.milk.sourcetracker <- plot_ordination(ps.milk.sourcetracker, ord.bray.nmds.milk.sourcetracker, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("After SourceTracker (milk microbiome)")+
  annotate("text", x = -0.6, y = 0.2, label = "Stress = 0.163\nPERMANOVA P = 0.704\nVariance P = 0.755", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.milk.sourcetracker)


# 3.3 subset teat samples

ps.teat.sourcetracker <- subset_samples(ps.sourcetracker, Type %in% c("Teat canal", "Teat apex") &!X.SampleID == "MS11_S44")
ps.teat.sourcetracker <- prune_taxa(taxa_sums(ps.teat.sourcetracker)>0, ps.teat.sourcetracker)
ps.teat.sourcetracker

metadata.teat.sourcetracker <- sample_data(ps.teat.sourcetracker) %>% data.frame()

# PERMANOVA
bray.dist.teat.sourcetracker<-vegdist(otu_table(ps.teat.sourcetracker), method='bray') 

beta_div.teat.sourcetracker<-pairwise.adonis2(bray.dist.teat.sourcetracker ~ Type, data=metadata.teat.sourcetracker, permutations=999)
beta_div.teat.sourcetracker

anova(betadisper(bray.dist.teat.sourcetracker, metadata.teat.sourcetracker$Type))

# ordination
ord.bray.nmds.teat.sourcetracker <- ordinate(ps.teat.sourcetracker, "NMDS", "bray", k=2, trymax = 1000)
ord.bray.nmds.teat.sourcetracker

ord.pnmds.teat.sourcetracker <- plot_ordination(ps.teat.sourcetracker, ord.bray.nmds.teat.sourcetracker, color = "Type") +
  geom_point(size = 5.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5))+
  ggtitle("After SourceTracker (udder microbiome)")+
  annotate("text", x = -0.9, y = 0.9, label = "Stress = 0.151\nPERMANOVA P = 0.371\nVariance P = 0.125", hjust = 0)+
  scale_shape_manual(values = c(17, 16))

print(ord.pnmds.teat.sourcetracker)


## merge plots
p.nmds.all <- (ord.pnmds + ord.pnmds.milk + ord.pnmds.teat) / (ord.pnmds.decontam + ord.pnmds.milk.decontam + ord.pnmds.teat.decontam) / (ord.pnmds.sourcetracker + ord.pnmds.milk.sourcetracker + ord.pnmds.teat.sourcetracker) +
  plot_annotation(tag_levels = "A")
p.nmds.all

ggsave(file.path(path.figure, "Figure 6. NMDS-all.png"), p.nmds.all, width=15, height=10, dpi=600)
