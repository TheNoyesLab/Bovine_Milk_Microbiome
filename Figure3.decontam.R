library(decontam)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(pairwiseAdonis)
library(patchwork)
library(microbiome)

TYPE <- c(
  "Cisternal milk", 
  "Teat canal", 
  "Teat apex", 
  "Stripped milk"
)

# Set up paths
path.rds <- "RDS/"
path.fig <- "Figures/"

# Load phyloseq object
ps <- readRDS(file.path(path.rds, "phyloseq.rds"))
ps

sample_data(ps)$SampleOrControl[sample_data(ps)$SampleOrControl == "Neither"] <- "Sample"
sample_data(ps)$CopyNumber <- as.numeric(sample_data(ps)$CopyNumber)

### Prepare data for 'frequency' method
exclude_controls <- c("Water", "Extraction")

ps.samples <- subset_samples(ps, !Matrix %in% exclude_controls) # remove control samples
ps.samples <- prune_taxa(taxa_sums(ps.samples) > 0, ps.samples) # remove taxa with no counts
ps.samples

# Run frequency method
# DNA contrentration is stored in the CopyNumber column
contamdf.freq <- isContaminant(ps.samples, method="frequency", conc="CopyNumber", threshold = 0.5)

# Show number of sequence features classified as contaminants
# FALSE: No. of sequence features that were not classified as contaminants (14,120)
# TRUE:  No. of sequence features that were classified as contaminants (422)
# Note: Most of the sequence features classified as contaminants will be very low prevalence
table(contamdf.freq$contaminant)

# Vector of sequece feature indices in otu table that are contaminants
# Note: Some of the most abundant sequence features were classified as contaminants
contam_feature_indices <- which(contamdf.freq$contaminant == TRUE)
head(contam_feature_indices)

#select contaminants from otu table
contam_seq_freq <- rownames(contamdf.freq[contamdf.freq$contaminant == "TRUE",])

contamdf.freqMod <- contamdf.freq %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+")))
        )

contamdf.freqMod$Prevalence <- factor(contamdf.freqMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

theme_set(theme_bw())
# Part of the Area Green Tableau color palette (Version 9.x)
# Similar to Davis et al., (2018)
scaleColorFillManualFrequency <-
  scale_fill_manual(
    values = 
      c(
        "2"    = "#dbe8b4",
        "3-5"  = "#9ad26d",
        "6-10" = "#6cae59",
        "11+"  = "#4a8c1c"
      )
  )

# Plot score statistics output from 'frequency' method
freqScorePlot <- ggplot(contamdf.freqMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method")                                    +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

freqScorePlot


### Prepare data for 'prevalence' method
sample_data(ps)$is.neg <- sample_data(ps)$SampleOrControl == "Control"

# Run 'prevalence' method
# Information about whether the sample ir a control or true sample is stored in 'is.neg' variable
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)

# Show number of sequence features classified as contaminants
# FALSE: No. of sequence features that were not classified as contaminants (14,559)
# TRUE:  No. of sequence features that were classified as contaminants (95)
table(contamdf.prev$contaminant)

# Vector of sequece feature indices in otu table that are contaminants
# Note: Some of the most abundant sequence features were classified as contaminants
contam_feature_indices_prev <- which(contamdf.prev$contaminant == TRUE)
head(contam_feature_indices_prev)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleOrControl == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleOrControl == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
pa.plot <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  theme(legend.position = "bottom")+
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)")
pa.plot

### Prepare data for 'prevalence' method, H2O is control
ps.no.extraction <- subset_samples(ps, !Type %in% "Extraction") # remove extraction samples
ps.no.extraction <- prune_taxa(taxa_sums(ps.no.extraction) > 0, ps.no.extraction) # remove taxa with no counts
ps.no.extraction
sample_data(ps.no.extraction)$is.neg <- sample_data(ps.no.extraction)$SampleOrControl == "Control"

# Run 'prevalence' method
# Information about whether the sample ir a control or true sample is stored in 'is.neg' variable
contamdf.prev.water <- isContaminant(ps.no.extraction, method="prevalence", neg="is.neg", threshold=0.5)

# Show number of sequence features classified as contaminants
# FALSE: No. of sequence features that were not classified as contaminants (14,547)
# TRUE:  No. of sequence features that were classified as contaminants (11)
table(contamdf.prev.water$contaminant)

# Vector of sequece feature indices in otu table that are contaminants
# Note: Some of the most abundant sequence features were classified as contaminants
contam_feature_indices_prev_water <- which(contamdf.prev.water$contaminant == TRUE)
head(contam_feature_indices_prev_water)

contamdf.prev.waterMod <- contamdf.prev.water %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                             ifelse(prev > 2 & prev <= 5, "3-5",
                                    ifelse(prev >= 6 & prev <= 10, "6-10", "11+")))
  )
contamdf.prev.waterMod$Prevalence <- factor(contamdf.prev.waterMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))


prev.waterScorePlot <- ggplot(contamdf.prev.waterMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method (library)")                         +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank()
  )
prev.waterScorePlot

#select seq name of contaminants identified by prevalence method using water as control
contam_seq_prev_water <- rownames(contamdf.prev.water[contamdf.prev.water$contaminant == "TRUE",])

### Prepare data for 'prevalence' method, Extraction is control
ps.no.library <- subset_samples(ps, !Type %in% "Library") # remove library samples
ps.no.library <- prune_taxa(taxa_sums(ps.no.library) > 0, ps.no.library) # remove taxa with no counts
ps.no.library
sample_data(ps.no.library)$is.neg <- sample_data(ps.no.library)$SampleOrControl == "Control"

# Run 'prevalence' method
# Information about whether the sample ir a control or true sample is stored in 'is.neg' variable
contamdf.prev.extraction <- isContaminant(ps.no.library, method="prevalence", neg="is.neg", threshold=0.5)

# Show number of sequence features classified as contaminants
# FALSE: No. of sequence features that were not classified as contaminants (14,527)
# TRUE:  No. of sequence features that were classified as contaminants (111)
table(contamdf.prev.extraction$contaminant)

# Vector of sequece feature indices in otu table that are contaminants
# Note: Some of the most abundant sequence features were classified as contaminants
contam_feature_indices_prev_extraction <- which(contamdf.prev.extraction$contaminant == TRUE)
head(contam_feature_indices_prev_extraction)

contamdf.prev.extractionMod <- contamdf.prev.extraction %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                             ifelse(prev > 2 & prev <= 5, "3-5",
                                    ifelse(prev >= 6 & prev <= 10, "6-10", "11+")))
  )
contamdf.prev.extractionMod$Prevalence <- factor(contamdf.prev.extractionMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))


prev.extractionScorePlot <- ggplot(contamdf.prev.extractionMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method (extraction)")                      +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.title.y = element_blank()
  )
prev.extractionScorePlot

#select seq names of contaminants identified by prevalence method using extraction as control
contam_seq_prev_extraction <- rownames(contamdf.prev.extraction[contamdf.prev.extraction$contaminant == "TRUE",])

### create noncontam ps project
#select unique seq names of contaminants identified by both frequency and prevalence methods separately 
contam_seq <- c(contam_seq_freq, contam_seq_prev_extraction, contam_seq_prev_water)
contam_seq <- unique(contam_seq)

# create ps project after removing contaminants
allTaxa = taxa_names(ps)
Taxa.decontam <- allTaxa[!(allTaxa %in% contam_seq)]
ps.decontam <- prune_taxa(Taxa.decontam, ps)
saveRDS(ps.decontam, file.path(path.rds, "ps.decontam.rds"))

ps.contam <- prune_taxa(contam_seq, ps)
saveRDS(ps.contam, file.path(path.rds, "ps.contam.rds"))

# plot 16S gene copy vs abundance of contaminants

SAMPLE_TYPES <- c(
  "Teat apex",
  "Teat canal",
  "Stripped milk",
  "Cisternal milk",
  "Air",
  "Blank",
  "Extraction",
  "Library"
)


my_palette<- c(
  "Teat apex" = "#1f77b4",
  "Teat canal" = "#ff7f0e",
  "Stripped milk" = "#2ca02c",
  "Cisternal milk" = "#d62728",
  "Library" = "#7f7f7f",
  "Extraction" = "#bcbd22",
  "Air" = "#539caf",
  "Blank" = "#ebc850"
)

metadata <- sample_data(ps.contam) %>% data.frame
counts <- sample_sums(ps.contam) %>% data.frame
colnames(counts) <- "counts"
counts$X.SampleID <- rownames(counts)
metadata <- left_join(metadata, counts, by = "X.SampleID")
metadata <- metadata[metadata$Type != "Library",]
metadata$Type <- factor(metadata$Type, levels = SAMPLE_TYPES)

p.contam.counts <- ggplot(data=metadata, aes(x=log(CopyNumber), y=counts/1000))+
                  geom_point(aes(color = Type), size = 2)+
                  geom_smooth(aes(group =1), method = "lm", se = TRUE) +
                   stat_cor(label.x = 7, label.y = 110) +
                   stat_regline_equation(label.x = 7, label.y = 100)+
                  xlab("Log (16S gene copy number)")+
                  ylab("Counts / 1000")+
                 scale_color_manual(values = my_palette)
p.contam.counts


# plot the abundance of contaminants
detection <- sum(taxa_sums(ps.contam))/70*0.05
pseq.contam <- aggregate_rare(ps.contam, level="Genus", detection = detection, prevalence = 0.2)

melt.contam <- psmelt(pseq.contam)
melt.contam$Type <- factor(melt.contam$Type, levels = SAMPLE_TYPES)

p.contam <- ggplot(melt.contam, aes(x=X.SampleID, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x") +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(10, "pt"),
    panel.grid.major.x = element_blank() )+
          scale_fill_brewer("Genus", palette = "Paired")
p.contam

# NMDS of contaminants
ps.contam.sample <- subset_samples(ps.contam, Type %in% c("Cisternal milk", "Stripped milk", "Teat apex", "Teat canal"))
ps.contam.sample <- prune_taxa(taxa_sums(ps.contam)>0, ps.contam.sample)

type_colors <- c(
  "Teat apex" = "#1f77b4",
  "Teat canal" = "#ff7f0e",
  "Stripped milk" = "#2ca02c",
  "Cisternal milk" = "#d62728"
)

ord.bray.nmds.contam <- ordinate(ps.contam.sample, "NMDS", "bray", k=3, trymax = 1000)
ord.pnmds.contam <- plot_ordination(ps.contam.sample, ord.bray.nmds.contam, color = "Type") +
  geom_point(size = 3.0) +   
  theme_bw() +  
  scale_color_manual(values = type_colors) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust=0.5))+
  annotate("text", x = -1.5, y = -1.5, label = "Stress = 0.124\nPERMANOVA P = 0.001",
           hjust = 0)
print(ord.pnmds.contam)

# PERMANOVA
bray.dist.contam<-vegdist(otu_table(ps.contam.sample), method='bray') 
bray.dist.contam

metadata.sample <- sample_data(ps.contam.sample) %>% data.frame()

beta_div.contam <-pairwise.adonis2(bray.dist.contam ~ Type, data=metadata.sample, permutations=999)
beta_div.contam


# merge plots

fig.decontam <- (freqScorePlot + prev.waterScorePlot + prev.extractionScorePlot) / (pa.plot + p.contam.counts + ord.pnmds.contam) / (p.contam ) + 
  plot_annotation(tag_levels = c("A"))

ggsave(file.path (path.fig, "Figure 3. Decontam_Contaminants.png"), fig.decontam, width = 16, height = 12, dpi=600)

