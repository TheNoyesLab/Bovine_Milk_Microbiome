library(data.table)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)
library(tidyr)
library(microbiome)
library(patchwork)

# set file paths
path.rds <- ("Analysis/RDS")

# read phyloseq object after decontam
ps <- readRDS("ps.decontam.rds")
target <- c("Teat apex", "Teat canal", "Stripped milk", "Cisternal milk")
ps.subset <- subset_samples(ps, Type %in% target)

## 1. remove contaminants from targeted animal samples using control as source

# read source count matrix
setwd("Analysis/results/Sourcetrackerr_contorl_animal_samples_40000/full_results")

air <- fread('control_animal_samples_Air_contributions.txt') %>% as.data.frame()
blank <- fread('control_animal_samples_Blank_contributions.txt') %>% as.data.frame()
ext <- fread('control_animal_samples_Extraction_contributions.txt') %>% as.data.frame()
library <- fread('control_animal_samples_Library_contributions.txt') %>% as.data.frame()
unknown <- fread('control_animal_samples_Unknown_contributions.txt') %>% as.data.frame()

row.names(air) <- air[,1]
air[,1] <- NULL

row.names(blank) <- blank[,1]
blank[,1] <- NULL

row.names(ext) <- ext[,1]
ext[,1] <- NULL

row.names(library) <- library[,1]
library[,1] <- NULL

row.names(unknown) <- unknown[,1]
unknown[,1] <- NULL

# create ps for contaminants identified by sourcetracker
# create the otu table for contaminants by multiplying sample depths
contam_sum <- air + blank + ext + library

sample.depths <- sample_sums(ps.subset)

otu_contam <- sweep(contam_sum,1,sample.depths,'*') %>% round(digits=0)

seq <- otu_table(otu_contam, taxa_are_rows = FALSE)
meta <- sample_data(ps.subset)
tax <- phyloseq::tax_table(ps.subset)
ps.contam.st.checked <- phyloseq(seq, tax, meta)

saveRDS(ps.contam.st, file.path(path.rds,"ps.contam.st.rds"))

# create ps for noncontam after sourcetracker

otu_noncontam <- sweep(unknown,1,sample.depths,'*') %>% round(digits=0)

seq.noncontam <- otu_table(otu_noncontam, taxa_are_rows = FALSE)
ps.decontam.st <- phyloseq(seq.noncontam, tax, meta)

saveRDS(ps.decontam.st, file.path(path.rds,"ps.decontam.st.rds"))

## 2. plot proportion of contaminants

meta <- sample_data(ps.subset) %>% data.frame()

# proportion of contaminants is stored at results$proportions

prop.type <- results$proportions %>% data.frame()
prop.type$Type <- meta$Type
prop.type.long <- prop.type %>% pivot_longer(-Type, names_to = "Source", values_to = "Proportion")

prop <- results$proportions %>% data.frame()
prop$X.SampleID <- meta$X.SampleID 

prop.long <- prop %>% pivot_longer(-X.SampleID,names_to="Source",values_to="Proportion")
prop.long$Type <- prop.type.long$Type

type_colors <-
  c(
    "Teat apex" = "#1f77b4",
    "Teat canal" = "#ff7f0e",
    "Cisternal milk" = "#d62728",
    "Library" = "#7f7f7f",
    "Extraction" = "#AB1866",
    "Air" = "#539caf",
    "Blank" = "#ebc850",
    "Unknown" = "#7AA874"
  )

p.prop.contam <- ggplot(data=prop.long, mapping = aes(x = X.SampleID, y = Proportion, fill = Source))+
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank())+
  scale_fill_manual(values = type_colors)

p.prop.contam


# 3. plot abundance of contaminants
ps.contam.st.melt <- psmelt(ps.contam.st.checked)
detection <- sum(taxa_sums(ps.contam.st.checked))/55

# genera with <5% relative abundance and <25% prevalence were grouped to Others
pseq.contam <- aggregate_rare(ps.contam.st.checked, level="Genus", detection = detection*0.05, prevalence = 0.25)

melt.pseq.contam <- psmelt(pseq.contam)

p.contam <-  ggplot(melt.pseq.contam, aes(x = X.SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(10, "pt"),
    panel.grid.major.x = element_blank()
  )+
  scale_fill_brewer("Genus", palette = "Paired")
p.contam

# 4. plot NMDS ordination 
type.colors <- c(
  "Teat apex" = "#1f77b4",
  "Teat canal" = "#ff7f0e",
  "Stripped milk" = "#2ca02c",
  "Cisternal milk" = "#d62728"
)

ord.bray.nmds <- ordinate(ps.contam.st.checked, "NMDS", "bray", k=2, trymax = 1000)
ord.pnmds.contam.st <- plot_ordination(ps.contam.st.checked, ord.bray.nmds, color = "Type") +
  geom_point(size = 4.0) +   
  theme_bw() +  
  scale_color_manual(values = type.colors) +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 14)) +
  labs(color = "Contaminants")+
  annotate("Text", x=-1, y=1.5, label = "Stress = 0.133\nPERMANOVA P = 0.001", hjust=0)
print(ord.pnmds.contam.st)

# PERMANOVA
otu.contam.st <- otu_table(ps.contam.st.checked)
meta <- meta %>% data.frame

permanova.contam.st.adonis <- adonis2(otu.contam.st ~ Type,
                               data = meta, permutations=999, method = "bray")
permanova.contam.st <- pairwise.adonis2(otu.contam.st ~ Type,
                                      data = meta, permutations=999, method = "bray")
permanova.contam.st

## 5. check contaminants in milk sample using control, teat apex and teat cannal as sources

results.qm <- readRDS(file.path(path.rds, "sourcetracker_quarter_milk.rds"))
results.cm <- readRDS(file.path(path.rds, "sourcetracker_cisternal_milk.rds"))

# add sample type to the stripped milk proportion table
qm<-as.data.frame(results.qm$proportions)
qm$X.SampleID <- rownames(qm)
qm.long <- pivot_longer(qm, cols = -X.SampleID, names_to = "Source", values_to = "Proportion")
qm.long$Type<- rep("Stripped milk", times =104)

# add sample type to the cisternal milk proportion table
cm <- as.data.frame(results.cm$proportions)
cm$X.SampleID <- rownames(cm)
cm.long <- pivot_longer(cm, cols = -X.SampleID, names_to = "Source", values_to = "Proportion")
cm.long$Type<- rep("Cisternal milk", times =98)

# merge cm and qm proportion table
st_milk_prop <- rbind(qm.long, cm.long)

st_milk_prop$Source <- factor(st_milk_prop$Source, levels=c("Air", "Blank", "Extraction", "Library",  "Teat apex", "Teat canal", "Cisternal milk", "Unknown"))

p.milk <- ggplot(data=st_milk_prop, mapping = aes(x = X.SampleID, y = Proportion, fill = Source))+
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank())+
  scale_fill_manual(values = type_colors)

p.milk

# merge plots
p.contam.st <- (((p.prop.contam) + (p.milk)) / ((p.contam) + (ord.pnmds.contam.st))) + plot_annotation(tag_levels = "A")
p.contam.st

ggsave(file.path(path.figure, "Figure 4. sourcetracker_contaminants_per_sample.png"), p.contam.st, width = 12, height = 10, dpi=600)
