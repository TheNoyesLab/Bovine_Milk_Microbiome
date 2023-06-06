library(dada2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gghalves)

## check positive mock sample
ps.positive <- readRDS ("ps.positive.rds"))

# mock community species
mock <- c("Pseudomonas", "Escherichia/Shigella", "Salmonella", "Lactobacillus", 
          "Enterococcus", "Staphylococcus", "Listeria", "Bacillus", "Saccharomyces", "Cryptococcus")

# subset the taxa from positive phyloseq object coming from mock community
mock.melt <- psmelt(subset_taxa(ps.positive, Genus %in% mock))
mock.melt$Genus <- as.character(mock.melt$Genus)
mock.melt$Genus <- ifelse(mock.melt$Genus == "Escherichia/Shigella", "Escherichia", mock.melt$Genus)
mock.melt$Genus <- factor(mock.melt$Genus, levels = c("Listeria", "Pseudomonas", "Bacillus", "Escherichia", "Salmonella", "Lactobacillus", "Staphylococcus"))

# plot the composition of positive control
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

## check sequencing depth and 16S copy numbers

ps <- readRDS("phyloseq.rds"))

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


# statistical analysis using glm, ANOVA and emmeans
glm.sample.CopyN <- glm(log(CopyNumber) ~ Type, data = sdata_df)
summary(glm.sample.CopyN)
car::Anova(glm.sample.CopyN)
emmeans::emmeans(glm.sample.CopyN, pairwise~Type)


## merge panel figure
p.quality <- ((p.reads | p.copy) / (p.positive)) + plot_annotation(tag_levels = "A") + plot_layout(height=c(2,1))
p.quality

ggsave("Figure 2. Depth_16S copy_Positive.png"), width=10, height=8, p.quality, dpi=600)
