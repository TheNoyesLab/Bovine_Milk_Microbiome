library(dada2)
library(dplyr)

##Set up working directories

#Define the working directory and file paths.
path.fastq <- "Data" # parent directory for raw data
path.meta <- "Analysis/Meta" # directory for metadata
path.rds <- "Analysis/RDS"
path.rdp <- "Analysis/reference/silva_nr99_v138.1_train_set.fa.gz" # directory for reference database
path.mock.fasta <- "Analysis/reference/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/zymo_logdist_update.fasta" # directory for mock genome database

##Filtering and Trimming

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path.fastq, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path.fastq, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:2); sample.names <- apply(sample.names, 2, paste0, collapse = "_")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.fastq, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.fastq, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,17), truncLen=c(270,220),
                     maxN=0, maxEE=c(3,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


## Learn erros 

errF <- learnErrors(filtFs, randomize = TRUE, multithread=TRUE)
errR <- learnErrors(filtRs, randomize = TRUE, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## Inference of sequence variants

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

## Merging of paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

##Create a sequence table

seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:258]
sum(seqtab2)/sum(seqtab)

## Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

##Track read retention through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


#remove UMGC internal standard, positive and one quarter milk sample with low reads
seqtab.nochim <- subset(seqtab.nochim, !(rownames(seqtab.nochim) %in% c("MS55_S18", "UMGC_304", "MS67_S42")))
dim(seqtab.nochim)

## Assign taxonomy without positive control

taxa <- assignTaxonomy(seqtab.nochim, path.rdp, multithread=TRUE)

## Create phyloseq object without positive control
seq <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)

tax <- tax_table(taxa)

# read metadata
meta <- read.csv(file.path(path.meta, 'meta.csv'))
row.names(meta) <- meta$X.SampleID
meta <- subset(meta, !X.SampleID %in% c("MS55_S18","UMGC_304", "MS67_S42"))
# Check sampleNames of seq and meta
identical(sample_names(seq), sort(sample_names(meta))) # TURE

# Generate phyloseq object without positive control
ps <- phyloseq(seq, tax, meta)

# Change sequence headers
taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

saveRDS(ps, file.path(path.rds,'phyloseq.rds'))

# Assign taxonomy to positive control
seq.positive <- subset(seqtab.nochim, rownames(seqtab.nochim) %in% "MS67_S42")

taxa.positive <- assignSpecies(seqtab.positive, path.mock.fasta)

meta.positive <- subset(meta, rownames(meta) %in% "MS67_S42")

# Generate phyloseq object for positive control
ps.postive <- phyloseq(seq.postive, tax.postive, meta.postive)
saveRDS(ps.positive, file.path(path.rds,"ps.positive.rds"))


