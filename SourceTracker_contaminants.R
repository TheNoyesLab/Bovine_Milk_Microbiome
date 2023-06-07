library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(MicrobiotaProcess)
library(tidyr)

# set up paths

path.rds <- ("/Analysis/RDS")

# read phyloseq object after decontam
ps.decontam <- readRDS("RDS/ps.decontam.rds")
metadata <- sample_data(ps.decontam)
metadata <- data.frame(metadata)

# plot rarefaction curve
ps.sample <- subset_samples(ps.decontam, Type %in% c("Stripped milk", "Cisternal milk", "Teat canal", "Teat apex"))

rareres <- get_rarecurve(obj=ps.sample, chunks=400)
p_rare <- ggrarecurve(obj=rareres,
                      factorNames="Type",
                      shadow=FALSE,
                      indexNames=c("Observe")) +
  theme_bw()+
  theme(axis.text=element_text(size=8),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks = seq(0, 175000, 5000))
p_rare

# rarefaction was set at 40,000 when most samples reach the plateau
rarefaction = 40000

## sourcetracker was run for three times, we only provide the code for the 1st run below
## 1. sources are "Library", "Blank", "Air", "Extraction"; sink are "Cisternal milk", "Stripped milk", "Teat canal", "Teat apex"
## 2. sources are "Library", "Blank", "Air", "Extraction", "Teat canal", "Teat apex"; sink are "Cisternal milk", "Stripped milk"
## 3. sources are "Library", "Blank", "Air", "Extraction", "Teat canal", "Teat apex", "Cisternal milk"; sink are "Stripped milk" 
## The codes were sourced from <https://github.com/danknights/sourcetracker> and <https://github.com/lakarstens/ControllingContaminants16S/blob/master/Analyses/ControllingContaminants16S_SourceTrackerPrep.Rmd> 
##################################################
# 1. assign control as source, samples as sink
metadata <- metadata %>%
  mutate(SourceSink = case_when(
    Type %in% c("Library", "Blank", "Air", "Extraction") ~ "source",
    Type %in% c("Cisternal milk", "Stripped milk", "Teat canal", "Teat apex") ~ "sink"
  ))

# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header
otus <- otu_table(ps.decontam)
otus <- data.frame(otus)
otus <- as.matrix(otus)

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
st_otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

# create directory to store outcomes

outdir='Sourcetracker_contorl_animal_sampels_40000'
filebase='control_animal_samples'

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Type
if(is.element('CowId',colnames(metadata))) desc <- metadata$CowId

# load SourceTracker package
source('SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, full.results = TRUE, rarefaction_depth  =rarefaction)

## Export results
# get average of full results across all runs of sourcetracker
res.mean <- apply(results$full.results,c(2,3,4),mean)

# Get depth of each sample for relative abundance calculation
sample.depths <- apply(results$full.results[1,,,,drop=F],4,sum)

# create directory to store the results
subdir <- paste(outdir,'full_results',sep='/')
dir.create(subdir,showWarnings=FALSE, recursive=TRUE)

# write each environment as a separate file
for(i in 1:length(results$train.envs)){
  env.name <- results$train.envs[i]
  filename.fractions <- sprintf('%s/%s_%s_contributions.txt', subdir, filebase, env.name)
  res.mean.i <- res.mean[i,,]
  # handle the case where there is only one sink sample
  if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)
  
  # make rows be samples, columns be features
  res.mean.i <- t(res.mean.i)
  
  # ensure proper names are retained
  colnames(res.mean.i) <- colnames(st_otus)
  rownames(res.mean.i) <- results$samplenames
  
  # calculate and save relative abundance
  res.mean.i.ra <- sweep(res.mean.i,1,sample.depths,'/')
  sink(filename.fractions)
  cat('SampleID\t')
  write.table(res.mean.i.ra,quote=F,sep='\t')
  sink(NULL)
}

#generate summary plots
if(dim(results$draws)[2] > 1) {
  plot.types <- c('pie')
} else plot.types <- c('pie', 'bar')
envs<-metadata[rownames(results$proportions),'Type']
envs<-unlist(envs)
envs <- as.factor(envs)
labels = sprintf('%s_%s',envs, rownames(results$proportions))
plotixs <- sort(as.numeric(envs),index=TRUE)$ix
for(plot.type in plot.types){
  # plot each environment separately
  for(env in unique(envs)){
    plotixs <- which(envs == env)
    png(sprintf('%s/%s_%s_%s.png',outdir,filebase,plot.type,env),width=500,height=500)
    plot(results, type=plot.type, labels=labels, include.legend=TRUE, indices=plotixs)
    dev.off()    
  }
}

saveRDS(results, file.path(path.rds, "sourcetracker_contaminants_results.rds"))

############################################################
# 2. define quarter milk as sink and the rests are source
metadata<- metadata %>%
  mutate(SourceSink = case_when(
    Type %in% c("Teat canal", "Teat apex", "Cisternal milk", "Library", "Blank", "Air", "Extraction") ~ "source",
    Type==c("Quarter milk") ~ "sink"
  ))
# repeat the code above get the results.qm
saveRDS(results.qm, file.path("RDS/", "sourcetracker_results_quartermilk.rds"))
###################################################
#3. define cisternal milk as sink and the rests (exclude of quarter milk) are source
metadata<- metadata %>%
  mutate(SourceSink = case_when(
    Type %in% c("Teat canal", "Teat apex", "Library", "Blank", "Air", "Extraction") ~ "source",
    Type==c("Cisternal milk") ~ "sink"
  ))
# repeat the code above to the results.cm
saveRDS(results.cm, file.path("RDS/", "sourcetracker_results_cisternalmilk.rds"))
