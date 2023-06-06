library(dada2); packageVersion("dada2")
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("/Users/deng0291/Desktop/Projects/01.Cow_Milk_Microbiome")

path.fastq <- "Data/NonPMA_fastq"
path.mock.fasta <- "/Users/deng0291/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/reference/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/zymo_logdist_update.fasta"
path.meta <- "Analysis/Meta"
path.rds <- "Analysis/RDS"
path.fig <- "Analysis/Figures"
path.results <- "Analysis/results"

list.files(path.fastq)
# Forward and reverse fastq filenames have format: MS67_S42_R1_001.fastq.gz and MS67_S42_R2_001.fastq.gz
fnFs <- sort(list.files(path.fastq, pattern="MS67_S42_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path.fastq, pattern="MS67_S42_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:2); sample.names <- apply(sample.names, 2, paste0, collapse = "_")

plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.fastq, "filtered.positive", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.fastq, "filtered.positive", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,17), truncLen=c(270,220),
                     maxN=0, maxEE=c(3,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, randomize = TRUE, multithread=TRUE)
errR <- learnErrors(filtRs, randomize = TRUE, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:258]
sum(seqtab2)/sum(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), sum(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

taxa <- assignSpecies(seqtab.nochim, path.mock.fasta)
counts <- sum(seqtab.nochim) #185968

meta <- read.csv(file.path(path.meta, 'meta.csv'))
row.names(meta) <- meta$X.SampleID
meta <- subset(meta, rownames(meta) %in% "MS67_S42")

saveRDS(taxa, file.path(path.rds, "taxa_positive.rds"))
saveRDS(seqtab.nochim, file.path(path.rds, "asv_positive.rds"))

# Make a phyloseq project for positive 
library(phyloseq)
seq <- otu_table(seqtab.nochim, taxa_are_rows = TRUE)
tax <- tax_table(taxa)
meta <- sample_data(meta)
ps.positive <- phyloseq(seq, tax, meta)

saveRDS(ps.positive, file.path(path.rds, "ps.positive.rds"))

# check positive vs. standard
ps.positive <- readRDS (file.path(path.rds, "ps.positive.rds"))

mock <- c("Pseudomonas", "Escherichia/Shigella", "Salmonella", "Lactobacillus", 
          "Enterococcus", "Staphylococcus", "Listeria", "Bacillus", "Saccharomyces", "Cryptococcus")

mock.melt <- psmelt(subset_taxa(ps.positive, Genus %in% mock))
mock.melt$Genus <- as.character(mock.melt$Genus)
mock.melt$Genus <- ifelse(mock.melt$Genus == "Escherichia/Shigella", "Escherichia", mock.melt$Genus)
mock.melt$Genus <- factor(mock.melt$Genus, levels = c("Listeria", "Pseudomonas", "Bacillus", "Escherichia", "Salmonella", "Lactobacillus", "Staphylococcus"))
write.csv(mock.melt, file.path(path.results, "positive_mock.csv"))

# read the csv table added with mock community
mock_pos <- read.csv(file.path(path.results, "positive_mock_standard.csv"))

mock_pos$Genus <- factor(mock_pos$Genus, levels=c("Listeria", "Pseudomonas", "Bacillus", "Escherichia", "Salmonella", "Lactobacillus", "Enterococcus", "Staphylococcus", "Unknown"))

genusPalette <- c(Bacillus="#009e73", Enterococcus="#332d2d", Escherichia="#0072b2", Lactobacillus="#f0e442", Listeria="#ecb333", Pseudomonas="#d62728", Salmonella="#ec6e0b", Staphylococcus="#76c6f3", Unknown="#bababa") 

p.positive <- ggplot(data=mock_pos, aes(fill=Genus, x=Sample, y=Relative.abundance, width=.4)) +             theme_bw()+
            geom_bar(position = "stack", stat = "identity", alpha=0.8)+
            scale_fill_manual(values=genusPalette)+
            ylab("Relative abundance")+
            xlab("")+
            coord_flip()+
             theme(
               axis.text.x = element_text(size = 10),
               axis.title.y = element_text(size = 10),
               legend.position = "bottom",
               legend.key.size = unit(0.5, "cm")
                   )
p.positive <- print(p.positive)

## check sequencing depth and 16S copy numbers

ps <- readRDS(file.path(path.rds, "phyloseq.rds"))
ps

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

theme_set(theme_bw())
# Part of the Tableau 10 palette (Version 9.x)

      
ColorFillManual <- c(
        "Teat apex" = "#1f77b4",
        "Teat canal" = "#ff7f0e",
        "Stripped milk" = "#2ca02c",
        "Cisternal milk" = "#d62728",
        "Library" = "#7f7f7f",
        "Extraction" = "#bcbd22",
        "Air" = "#539caf",
        "Blank" = "#ebc850"
      )


sdata_df <- as(sample_data(ps), "data.frame")
sdata_df$Type <- factor(sdata_df$Type, levels = SAMPLE_TYPES)
sdata_df$PF_Clusters <- as.numeric(sdata_df$PF_Clusters)

## plot No. Reads by Sample Type
library(gghalves)

p.reads <- sdata_df %>% 
  ggplot(aes(x = Type, y = PF_Clusters/1000, fill = Type, color=Type))+
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
               color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  geom_half_point(position = position_nudge(x=-0.4), 
                      size=2, range_scale=0.4, alpha = 0.6) +
  ylab("Number of Reads / 1000") +
  ylim(0,320)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))
p.reads <- print(p.reads)


# statistical analysis suing glm, ANOVA and emmeans
## we cannot put CowID in the model, because some samples are not from cows
glm.sample.depth <- glm(PF_Clusters/1000 ~ Type, data = sdata_df)
summary(glm.sample.depth)
car::Anova(glm.sample.depth)
emmeans::emmeans(glm.sample.depth, pairwise~Type)


## plot copy No. by Sample Type
p.copy<- sdata_df %>% 
  ggplot(aes(x = Type, y = log(CopyNumber), fill = Type, color=Type))+
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F)+
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  geom_half_point(position = position_nudge(x=-0.4), 
                  size=2, range_scale=0.4, alpha = 0.6) +
  ylab("Log (16S copy number)") +
  ylim(0,16)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))
p.copy <- print(p.copy)


# statistical analysis suing glm, ANOVA and emmeans
glm.sample.CopyN <- glm(log(CopyNumber) ~ Type, data = sdata_df)
summary(glm.sample.CopyN)
car::Anova(glm.sample.CopyN)
emmeans::emmeans(glm.sample.CopyN, pairwise~Type)


## merge panel figure
p.quality <- ((p.reads | p.copy) / (p.positive)) + plot_annotation(tag_levels = "A") + plot_layout(height=c(2,1))
p.quality

ggsave(file.path(path.fig, "Figure 2. Depth_16S copy_Positive.png"), width=10, height=8, p.quality, dpi=600)


save.image("~/Desktop/Projects/01.Cow_Milk_Microbiome/Analysis/RData/Check positive.RData")
