library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(patchwork)
library(gghalves)

# set file path
path.rds <- "Analysis/RDS"
path.figure <- ("Analysis/Figures")

# set order for sample type
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

## 1. loapd ps pre-decontam

ps <- readRDS("phyloseq.rds")
ps

metadata <- sample_data(ps) %>% data.frame()
metadata$Type <- factor(metadata$Type, levels = SAMPLE_TYPES)

depth <- sample_sums(ps) %>% data.frame 
depth$X.SampleID <- metadata$X.SampleID

# alpha diversity
tabAlphaDiv <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "InvSimpson"))
tabAlphaDiv$X.SampleID <- row.names(tabAlphaDiv)
metadata <- metadata %>% left_join(tabAlphaDiv, by = "X.SampleID")

### 1.1 richness
p.richness <- metadata %>% 
  ggplot(aes(x = Type, y = Observed, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("Richness") +
  labs(title="Before decontam")+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,2000)
p.richness

# statistical analysis 
glm.richness <- glm (Observed ~ Type , data = metadata)
summary(glm.richness)
car::Anova(glm.richness, Type=III)
emmeans.richness <- emmeans::emmeans(glm.richness, pairwise~Type)

richness.contra <- emmeans.richness$contrasts %>% data.frame
richness.mean <- emmeans.richness$emmeans %>% data.frame


### 1.2 shannon
p.shannon <- metadata %>% 
  ggplot(aes(x = Type, y = Shannon, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("Shannon") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0.9,7.1)
p.shannon

# statistical analysis using glm, ANOVA and emmeans
glm.shannon <- lmer(Shannon ~ Type, data = metadata)
summary(glm.shannon)
car::Anova(glm.shannon)
emmeans.shannon <- emmeans::emmeans(glm.shannon, pairwise~Type)

shannon.contra <- emmeans.shannon$contrasts %>% data.frame
shannon.mean <- emmeans.shannon$emmeans %>% data.frame



### 1.3 InvSimpson
p.InvSimpson <- metadata %>% 
  ggplot(aes(x = Type, y = InvSimpson, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("Inver Simpson") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,200)
p.InvSimpson

# statistical analysis suing glm, ANOVA and emmeans
glm.InvSimpson <- glm(InvSimpson ~ Type, data = metadata)
summary(glm.InvSimpson)
car::Anova(glm.InvSimpson)
emmeans::emmeans(glm.InvSimpson, pairwise~Type)


## 2. after-decontam

ps.decontam <- readRDS(file.path(path.rds, "ps.decontam.rds"))
ps.decontam <- subset_samples(ps.decontam, !Type %in% c("Library", "Extraction"))
ps.decontam <- prune_taxa(taxa_sums(ps.decontam)>0, ps.decontam)

metadata.decontam <- sample_data(ps.decontam) %>% data.frame()
metadata.decontam$Type <- factor(metadata.decontam$Type, levels = SAMPLE_TYPES)

depth.decontam <- sample_sums(ps.decontam) %>% data.frame
depth.decontam$X.SampleID <- metadata.decontam$X.SampleID

# alpha diversity
tabAlphaDiv.decontam <- phyloseq::estimate_richness(ps.decontam, measures = c("Observed", "Shannon", "InvSimpson"))
tabAlphaDiv.decontam$X.SampleID <- row.names(tabAlphaDiv.decontam)
metadata.decontam <- metadata.decontam %>% left_join(tabAlphaDiv.decontam, by = "X.SampleID")

### 2.1 richness
p.richness.decontam <- metadata.decontam %>% 
  ggplot(aes(x = Type, y = Observed, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("")+
  labs(title="After decontam")+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,2000)
p.richness.decontam

# statistical analysis suing glm, ANOVA and emmeans
glm.richness.decontam <- glm(Observed ~ Type, data = metadata.decontam)
summary(glm.richness.decontam)
car::Anova(glm.richness.decontam)
emmeans::emmeans(glm.richness.decontam, pairwise~Type)


### 2.2 shannon
p.shannon.decontam <- metadata.decontam %>% 
  ggplot(aes(x = Type, y = Shannon, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0.9,7.1)
p.shannon.decontam

# statistical analysis suing glm, ANOVA and emmeans
glm.shannon.decontam <- glm(Shannon ~ Type, data = metadata.decontam)
summary(glm.shannon.decontam)
car::Anova(glm.shannon.decontam)
emmeans::emmeans(glm.shannon.decontam, pairwise~Type)


### 2.3 InvSimpson
p.InvSimpson.decontam <- metadata.decontam %>% 
  ggplot(aes(x = Type, y = InvSimpson, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,200)
p.InvSimpson.decontam

# statistical analysis suing glm, ANOVA and emmeans
glm.InvSimpson.decontam <- glm(InvSimpson ~ Type, data = metadata.decontam)
summary(glm.InvSimpson.decontam)
car::Anova(glm.InvSimpson.decontam)
emmeans::emmeans(glm.InvSimpson.decontam, pairwise~Type)

## 3. after-sourcetracker

ps.sourcetracker <- readRDS(file.path(path.rds, "ps.decontam.st.rds"))
ps.sourcetracker <- prune_taxa(taxa_sums(ps.sourcetracker)>0, ps.sourcetracker)
ps.sourcetracker

metadata.sourcetracker <- sample_data(ps.sourcetracker) %>% data.frame()
metadata.sourcetracker$Type <- factor(metadata.sourcetracker$Type, levels = SAMPLE_TYPES)

depth.sourcetracker <- sample_sums(ps.sourcetracker) %>% data.frame
depth.sourcetracker$X.SampleID <- metadata.sourcetracker$X.SampleID

# alpha diversity
tabAlphaDiv.sourcetracker <- phyloseq::estimate_richness(ps.sourcetracker, measures = c("Observed", "Shannon", "InvSimpson"))
tabAlphaDiv.sourcetracker$X.SampleID <- row.names(tabAlphaDiv.sourcetracker)
metadata.sourcetracker <- metadata.sourcetracker %>% left_join(tabAlphaDiv.sourcetracker, by = "X.SampleID")

### 3.1 richness
p.richness.sourcetracker <- metadata.sourcetracker %>% 
  ggplot(aes(x = Type, y = Observed, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("")+
  labs(title="After SourceTracker")+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,2000)
p.richness.sourcetracker

# statistical analysis suing glm, ANOVA and emmeans
glm.richness.sourcetracker <- glm(Observed ~ Type, data = metadata.sourcetracker)
summary(glm.richness.sourcetracker)
car::Anova(glm.richness.sourcetracker)
emmeans::emmeans(glm.richness.sourcetracker, pairwise~Type)


### 3.2 shannon
p.shannon.sourcetracker <- metadata.sourcetracker %>% 
  ggplot(aes(x = Type, y = Shannon, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0.9,7.1)
p.shannon.sourcetracker

# statistical analysis suing glm, ANOVA and emmeans
glm.shannon.sourcetracker <- glm(Shannon ~ Type, data = metadata.sourcetracker)
summary(glm.shannon.sourcetracker)
car::Anova(glm.shannon.sourcetracker)
emmeans::emmeans(glm.shannon.sourcetracker, pairwise~Type)


### 3.3 InvSimpson
p.InvSimpson.sourcetracker <- metadata.sourcetracker %>% 
  ggplot(aes(x = Type, y = InvSimpson, fill = Type, color = Type)) +
  scale_color_manual (values = ColorFillManual)+
  scale_fill_manual (values = ColorFillManual)+
  geom_half_violin(position = position_nudge(x = 0.15, y = 0),
                   color=NA, alpha = 0.6, scale = "width", trim = F, side = "R")+
  geom_half_point(position = position_nudge(x=-0.4),
                  size=2, range_scale=0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha=0.6) +
  xlab("Sample Type") +
  ylab("") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  ylim(0,200)
p.InvSimpson.sourcetracker

# statistical analysis suing glm, ANOVA and emmeans
glm.InvSimpson.sourcetracker <- glm(InvSimpson ~ Type, data = metadata.sourcetracker)
summary(glm.InvSimpson.sourcetracker)
car::Anova(glm.InvSimpson.sourcetracker)
emmeans::emmeans(glm.InvSimpson.sourcetracker, pairwise~Type)


### merge figures

p.alpha.diversity <- ((p.richness / p.shannon / p.InvSimpson) | (p.richness.decontam / p.shannon.decontam / p.InvSimpson.decontam) | (p.richness.sourcetracker / p.shannon.sourcetracker / p.InvSimpson.sourcetracker)) + 
  plot_layout(widths=c(8, 7, 6))
p.alpha.diversity


ggsave(file.path(path.figure, "Figure 5. Alpha diversity.checked.png"), p.alpha.diversity, width = 9, height = 6, dpi=600)   
