library(dada2); packageVersion("dada2")

setwd("/Users/deng0291/Desktop/Projects/01.Cow_Milk_Microbiome/")

path.fastq <- "Data/NonPMA_fastq"
path.meta <- "Analysis/Meta"
path.rdp <- "Analysis/reference/silva_nr99_v138.1_train_set.fa.gz"
path.rdp.species <- "Analysis/reference/silva_species_assignment_v138.1.fa.gz"
path.rds <- "Analysis/RDS"
path.results <- "Analysis/results"

list.files(path.fastq)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path.fastq, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path.fastq, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:2); sample.names <- apply(sample.names, 2, paste0, collapse = "_")

plotQualityProfile(fnFs[1:8])
plotQualityProfile(fnRs[1:8])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.fastq, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.fastq, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
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
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(path.results, "track.csv"))

#remove UMGC internal standard, positive and one quarter milk sample with low reads
seqtab.nochim <- subset(seqtab.nochim, !(rownames(seqtab.nochim) %in% c("MS55_S18", "UMGC_304", "MS67_S42")))
dim(seqtab.nochim)


taxa <- assignTaxonomy(seqtab.nochim, path.rdp, multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

meta <- read.csv(file.path(path.meta, 'meta.csv'))
row.names(meta) <- meta$X.SampleID

saveRDS(taxa, file.path(path.rds, "taxa.rds"))
saveRDS(seqtab.nochim, file.path(path.rds, "asv.rds"))
saveRDS(meta, file.path(path.rds, "meta.rds"))
