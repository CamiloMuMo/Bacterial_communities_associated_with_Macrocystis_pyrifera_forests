### Data processing Illumina myseq paired end with DADA2
# Big data workflow. Each step is self consistent and it is possible to stop after each step
# R1 and R2 fastq reads must be in two different directories
# Create an R project in the working directory with the following subdirectories:
#                                         - R1 (contains R1)
#                                         - R2 (contains R2)
#                                         - dada.scripts (contains R.script for dada2)
#                                         - taxonomy (contains the databases to assign taxonomy)
message("Upload libraries")
library(dada2)
library(data.table)

message("STEP1: QUALITY CHECK")
message("File parsing")
pathF <- "R1"
pathR <- "R2"
fastqFsq <- sort(list.files(pathF, pattern="R1.fastq", full.names=TRUE))
fastqRsq <- sort(list.files(pathR, pattern="R2.fastq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fastqFsq), ".R1."), `[`, 1) # CHANGE WITH SEQUENCE NAMES. Extract sample names

# Reading quality profiles (F)
plotQualityProfile(fastqFsq[1:2])
# Reading quality profiles (R)
plotQualityProfile(fastqRsq[1:2])

message("STEP2: FILTERING")

message("File parsing")
pathF <- "R1"
pathR <- "R2"
fastqFs <- sort(list.files(pathF, pattern="R1.fastq"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq"))
dir.create("track_reads")

dir.create(file.path(pathF, "filtered1"))
filtpathF1 <- file.path(pathF, "filtered1")

dir.create(file.path(pathR, "filtered1"))
filtpathR1 <- file.path(pathR, "filtered1")

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

message("Filtering")
#THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
out1 <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF1, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR1, fastqRs),
              truncLen=c(250,150),
              maxN=0, maxEE=c(2,2), trimLeft=c(19,20), rm.phix=TRUE,
              compress=FALSE, multithread=TRUE) # Adjust filtering parameters if needed. On Windows set multithread=FALSE
write.table(out1, "track_reads/filter1.txt")

message("STEP3: INFER SEQUENCE VARIANTS")
message("File parsing")
pathF <- "R1"
pathR <- "R2"
filtpathF <- file.path(pathF, "filtered1")
filtpathR <- file.path(pathR, "filtered1")
filtFs <- list.files(filtpathF, pattern="R1.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R2.fastq", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), ".R1."), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), ".R2."), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

message("Learn the error rates and plot them")
dir.create(file.path(filtpathF, "error.rates_filter1"))
dir.create(file.path(filtpathR, "error.rates_filter1"))

errF <- learnErrors(filtFs, nbases=1e10, multithread=TRUE, randomize=TRUE)
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs, nbases=1e10, multithread=TRUE, randomize=TRUE, MAX_CONSIST = 20)
plotErrors(errR, nominalQ=TRUE)

message("Sample inference and merger of paired-end reads")
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR,minOverlap=10, maxMismatch=1)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table
dir.create("out_filter1")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "out_filter1/seqtab_filter1.rds")
getN <- function(x) sum(getUniques(x))
merged <- sapply(mergers, getN)
write.table(merged, "track_reads/merged_filter1.txt")

message("STEP4: MERGE RUNS, REMOVE CHIMERAS, ASSIGN TAXONOMY ###")
# Merge runs
 st1 <- readRDS(" container folder") <- insert your container folder

message("Remove chimeras")
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
saveRDS(seqtab, "seqtab_filter_nochim.rds")
nochim <- rowSums(seqtab)
write.table(nochim, "nochim_filter.txt")

message("Assign taxonomy with Silva")
taxa.silva <- assignTaxonomy(seqtab, file.path("silva_nr99_v138.2_train_set.fa.gz"), minBoot = 80, multithread=TRUE, verbose=TRUE)
species.silva <- addSpecies(taxa.silva, file.path("silva_species_assignment_v138.2.fa.gz"), verbose=TRUE) # Assign species
saveRDS(taxa.silva, file.path("taxa.silva.rds"))
saveRDS(species.silva, file.path("species.silva.rds"))

message("Save files, tack, taxa, species")
dir.create(file.path("final.data_filter1"))

seq.abundances <- as.data.frame(t(seqtab))
seq.abundances <- setDT(as.data.frame(seq.abundances), keep.rownames = TRUE)
colnames(seq.abundances)[colnames(seq.abundances) == c("rn")] <- c("ASV_seq")
seq.abundances <- setDT(as.data.frame(seq.abundances), keep.rownames = TRUE)
seq.abundances$rn <- paste0("ASV_", seq.abundances$rn)
colnames(seq.abundances)[colnames(seq.abundances) == c("rn")] <- c("ASV")

message("Silva")
taxa1.silva <- setDT(as.data.frame(taxa.silva), keep.rownames = TRUE)
colnames(taxa1.silva)[colnames(taxa1.silva) == c("rn")] <- c("ASV_seq")
taxa.abundance.silva <- merge(seq.abundances, taxa1.silva, by = "ASV_seq")
write.table(taxa.abundance.silva, file = file.path("taxa.silva.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

species1.silva <- setDT(as.data.frame(species.silva), keep.rownames = TRUE)
colnames(species1.silva)[colnames(species1.silva) == c("rn")] <- c("ASV_seq")
species.abundance.silva <- merge(seq.abundances, species1.silva, by = "ASV_seq")
write.table(species.abundance.silva, file = file.path("species.silva.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
