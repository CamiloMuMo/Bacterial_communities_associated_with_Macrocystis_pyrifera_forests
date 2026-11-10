# 01. Paired-end V4 reads to taxa.silva.xlsx ----
# Run these sections in order, with the folder containing this script as the working directory.
# Input: original primer-bearing FASTQs in R1/ and R2/, samples.tsv, and two SILVA FASTAs.
# References: Callahan et al. (2016), doi:10.1038/nmeth.3869; https://benjjneb.github.io/dada2/
library(dada2)       # * 1.32.0
library(ShortRead)   # * 1.62.0 (stream the original FASTQs for the primer audit)
library(writexl)     # * 1.5.4 (export the two worksheets)
library(ggplot2)     # * 4.0.3 (save DADA2 quality and error plots)

# 1. Match the forward and reverse reads ----
# One row per sequenced library, in the original DADA2 order. File names may be changed
# to the downloaded ENA names without changing ID, Library or row order.
samples <- read.delim("samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(samples) == 34L, !anyDuplicated(samples$ID), !anyDuplicated(samples$Library),
          !anyDuplicated(samples$R1), !anyDuplicated(samples$R2),
          sum(samples$Type == "Sediment") == 25L, sum(samples$Type == "Sea_Water") == 6L)
control_ids <- c("Blank1", "Blank2", "Blank3")
environmental <- samples$Type %in% c("Sediment", "Sea_Water")
stopifnot(setequal(samples$ID[!environmental], control_ids))
fastqFs <- file.path("R1", samples$R1)
fastqRs <- file.path("R2", samples$R2)
# The same files may be stored uncompressed or gzip-compressed.
fastqFs <- ifelse(!file.exists(fastqFs) & file.exists(paste0(fastqFs, ".gz")), paste0(fastqFs, ".gz"), fastqFs)
fastqRs <- ifelse(!file.exists(fastqRs) & file.exists(paste0(fastqRs, ".gz")), paste0(fastqRs, ".gz"), fastqRs)
# For archive-style names, the author's ENA run column can be filled once before release.
# Original service filenames still work when ENA_run is blank (including all three controls).
if ("ENA_run" %in% names(samples)) {
  known <- !is.na(samples$ENA_run) & nzchar(samples$ENA_run)
  stopifnot(!anyDuplicated(samples$ENA_run[known]))
  enaF <- file.path("R1", paste0(samples$ENA_run, "_1.fastq.gz"))
  enaR <- file.path("R2", paste0(samples$ENA_run, "_2.fastq.gz"))
  useF <- known & !file.exists(fastqFs) & file.exists(enaF)
  useR <- known & !file.exists(fastqRs) & file.exists(enaR)
  fastqFs[useF] <- enaF[useF]
  fastqRs[useR] <- enaR[useR]
}
missing_files <- c(fastqFs, fastqRs)[!file.exists(c(fastqFs, fastqRs))]
if (length(missing_files)) stop("FASTQ not found. Check the file names in samples.tsv:\n", paste(missing_files, collapse = "\n"))
stopifnot(!anyDuplicated(fastqFs), !anyDuplicated(fastqRs))
names(fastqFs) <- names(fastqRs) <- samples$Library
threads <- .Platform$OS.type != "windows"  # FALSE also runs in a single thread
dir.create("DADA2_output", showWarnings = FALSE)
dir.create("DADA2_output/filtered", showWarnings = FALSE)
dir.create("results/supplement", recursive = TRUE, showWarnings = FALSE)
# The author used SILVA 138.2. These are the public DADA2-formatted 138.2 filenames.
# If the locally archived originals have different names, change only these two paths.
training_file <- "Taxonomy/silva_nr99_v138.2_toGenus_trainset.fa.gz"
species_file  <- "Taxonomy/silva_v138.2_assignSpecies.fa.gz"
if (!all(file.exists(c(training_file, species_file)))) stop("Check the two SILVA 138.2 reference paths in section 1.")
write.table(data.frame(Reference = c(training_file, species_file),
                       MD5 = unname(tools::md5sum(c(training_file, species_file)))),
            "DADA2_output/reference_files.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
ggsave("DADA2_output/quality_R1.pdf", plotQualityProfile(unname(fastqFs[1:2])), width = 9, height = 5)
ggsave("DADA2_output/quality_R2.pdf", plotQualityProfile(unname(fastqRs[1:2])), width = 9, height = 5)

write.table(data.frame(ID = samples$ID, Library = samples$Library, Type = samples$Type,
                       R1 = unname(fastqFs), R2 = unname(fastqRs)),
            "DADA2_output/file_pairing.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 2. Locate primers BEFORE trimming: Table S2A ----
# Find the first full, exact IUPAC-compatible primer in the first 60 nucleotides.
# Positions are one-based. A mismatch or a partial primer is not an exact full match.
primer_patterns <- c(R1 = "GTG[CT]CAGC[AC]GCCGCGGTAA",               # 515F: GTGYCAGCMGCCGCGGTAA
                     R2 = "GGACTAC[ACGT][ACG]GGGT[AT]TCTAAT")        # 806R-B: GGACTACNVGGGTWTCTAAT
trim_bases <- c(R1 = 19L, R2 = 20L)
first_read_id <- matrix(NA_character_, nrow(samples), 2, dimnames = list(samples$ID, c("R1", "R2")))
primer_rows <- list()
shifted_reads <- list()
for (i in seq_len(nrow(samples))) {
  for (mate in c("R1", "R2")) {
    file <- if (mate == "R1") fastqFs[i] else fastqRs[i]
    key <- paste(samples$ID[i], mate, sep = ":")
    stream <- FastqStreamer(unname(file), n = 100000L)
    bins <- integer(61L)  # first bin = no match; start at nucleotide 1 = second bin
    shifted <- list()
    tryCatch({
      repeat {
        chunk <- yield(stream)
        if (!length(chunk)) break
        if (is.na(first_read_id[i, mate])) {
          header <- as.character(ShortRead::id(chunk)[1L])
          first_read_id[i, mate] <- sub("/[12]$", "", sub("[[:space:]].*$", "", header))
        }
        bases <- substr(as.character(sread(chunk)), 1L, 60L)
        pos <- as.integer(regexpr(primer_patterns[[mate]], bases, perl = TRUE))
        bins <- bins + tabulate(ifelse(pos < 0L, 1L, pos + 1L), nbins = 61L)
        select <- which(pos > 1L)
        if (length(select)) {
          # Save only these rare reads, to check which survived filtering and their actual prefix.
          identifiers <- sub("[[:space:]].*$", "", as.character(ShortRead::id(chunk))[select])
          substring_after_cut <- substring(bases[select], trim_bases[[mate]] + 1L,
                                           pos[select] + trim_bases[[mate]] - 1L)
          shifted[[length(shifted) + 1L]] <- data.frame(Read_ID = identifiers, Prefix = substring_after_cut)
        }
      }
    }, finally = close(stream))
    nonzero <- which(bins > 0L)
    primer_rows[[key]] <- data.frame(ID = samples$ID[i], Type = samples$Type[i], Mate = mate,
      Primer_start = ifelse(nonzero == 1L, NA_integer_, nonzero - 1L), Reads = bins[nonzero])
    shifted_reads[[key]] <- if (length(shifted)) do.call(rbind, shifted) else
      data.frame(Read_ID = character(), Prefix = character())
    stopifnot(!anyDuplicated(shifted_reads[[key]]$Read_ID))
  }
}
different <- is.na(first_read_id[, 1]) | is.na(first_read_id[, 2]) | first_read_id[, 1] != first_read_id[, 2]
if (any(different)) stop("Check paired FASTQ identifiers before denoising: ", paste(samples$ID[different], collapse = ", "))
primer_by_library <- do.call(rbind, primer_rows)
primer_by_library$Position <- ifelse(is.na(primer_by_library$Primer_start), "No exact match within first 60 nt",
                              ifelse(primer_by_library$Primer_start == 1L, "First nucleotide", "After first nucleotide"))
write.table(primer_by_library, "DADA2_output/primer_positions_by_library.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
primer_environment <- subset(primer_by_library, ID %in% samples$ID[environmental])
S2A <- aggregate(Reads ~ Mate + Position, primer_environment, sum)
S2A$Percent <- 100 * S2A$Reads / ave(S2A$Reads, S2A$Mate, FUN = sum)
S2A <- S2A[order(S2A$Mate, match(S2A$Position, c("First nucleotide", "After first nucleotide", "No exact match within first 60 nt"))), ]
write.table(S2A, "results/supplement/S2A_primer_positions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(S2A, row.names = FALSE)

# 3. Filter and trim the paired reads ----
filtFs <- file.path("DADA2_output/filtered", paste0(samples$ID, "_R1.fastq"))
filtRs <- file.path("DADA2_output/filtered", paste0(samples$ID, "_R2.fastq"))
names(filtFs) <- names(filtRs) <- samples$Library
filtered <- filterAndTrim(fwd = unname(fastqFs), filt = unname(filtFs),
                          rev = unname(fastqRs), filt.rev = unname(filtRs),
                          truncLen = c(250, 150), trimLeft = c(19, 20),
                          maxN = 0, maxEE = c(2, 2), rm.phix = TRUE,
                          compress = FALSE, multithread = threads)
stopifnot(identical(basename(rownames(filtered)), basename(unname(fastqFs))), all(filtered[, "reads.out"] > 0))
# No new trimming rule is introduced. Full primers starting at nucleotide 1 lie entirely
# within trimLeft. For exact shifted hits the cut can leave a suffix; audit it in surviving reads.
posttrim_rows <- list()
for (i in seq_len(nrow(samples))) {
  for (mate in c("R1", "R2")) {
    key <- paste(samples$ID[i], mate, sep = ":")
    lookup <- shifted_reads[[key]]
    surviving <- logical(nrow(lookup))
    prefix_matches <- logical(nrow(lookup))
    n_filtered <- full_after <- 0L
    filtered_bases <- 0
    file <- if (mate == "R1") filtFs[i] else filtRs[i]
    stream <- FastqStreamer(unname(file), n = 100000L)
    tryCatch({
      repeat {
        chunk <- yield(stream)
        if (!length(chunk)) break
        bases <- substr(as.character(sread(chunk)), 1L, 60L)
        n_filtered <- n_filtered + length(bases)
        filtered_bases <- filtered_bases + sum(width(sread(chunk)))
        full_after <- full_after + sum(grepl(primer_patterns[[mate]], bases, perl = TRUE))
        if (nrow(lookup)) {
          ids <- sub("[[:space:]].*$", "", as.character(ShortRead::id(chunk)))
          ix <- match(ids, lookup$Read_ID)
          k <- which(!is.na(ix))
          if (length(k)) {
            surviving[ix[k]] <- TRUE
            prefix_matches[ix[k]] <- startsWith(bases[k], lookup$Prefix[ix[k]])
          }
        }
      }
    }, finally = close(stream))
    stopifnot(n_filtered == unname(filtered[i, "reads.out"]))
    posttrim_rows[[key]] <- data.frame(ID = samples$ID[i], Type = samples$Type[i], Mate = mate,
      Raw_exact_shifted = nrow(lookup), Shifted_retained = sum(surviving),
      Expected_prefix_observed = sum(prefix_matches), Filtered_reads = n_filtered, Filtered_bases = filtered_bases,
      Filtered_with_full_primer_first_60_nt = full_after)
  }
}
posttrim <- do.call(rbind, posttrim_rows)
raw_totals <- aggregate(Reads ~ ID + Mate, primer_by_library, sum)
stopifnot(all(raw_totals$Reads == filtered[match(raw_totals$ID, samples$ID), "reads.in"]))
write.table(posttrim, "DADA2_output/after_trim_by_library.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# A negative full-motif scan cannot prove the absence of a partial or mismatched primer.
rm(shifted_reads)

# 4. Learn errors, infer ASVs and merge paired ends ----
# August, March, water, sediments and all three controls were processed together.
# Keep the original seed once, before BOTH calls; do not reset it between directions.
set.seed(100)
errF_log <- capture.output(errF <- learnErrors(filtFs, nbases = 1e10, randomize = TRUE, multithread = threads, verbose = TRUE))
writeLines(errF_log, "DADA2_output/error_learning_R1.txt")
ggsave("DADA2_output/errors_R1.pdf", plotErrors(errF, nominalQ = TRUE), width = 8, height = 7)
errR_log <- capture.output(errR <- learnErrors(filtRs, nbases = 1e10, randomize = TRUE, MAX_CONSIST = 20, multithread = threads, verbose = TRUE))
writeLines(errR_log, "DADA2_output/error_learning_R2.txt")
# If all filtered bases together fall below nbases, learnErrors must read all supplied libraries.
error_pool <- posttrim[, c("ID", "Type", "Mate", "Filtered_reads", "Filtered_bases")]
error_pool$All_supplied_libraries_needed <- ave(error_pool$Filtered_bases, error_pool$Mate, FUN = sum) < 1e10
write.table(error_pool, "DADA2_output/error_learning_libraries.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
ggsave("DADA2_output/errors_R2.pdf", plotErrors(errR, nominalQ = TRUE), width = 8, height = 7)
mergers <- setNames(vector("list", nrow(samples)), samples$Library)
denoised_F <- denoised_R <- numeric(nrow(samples))
for (i in seq_len(nrow(samples))) {
  message("Processing ", samples$ID[i])
  derepF <- derepFastq(filtFs[[i]])
  ddF <- dada(derepF, err = errF, multithread = threads)
  derepR <- derepFastq(filtRs[[i]])
  ddR <- dada(derepR, err = errR, multithread = threads)
  denoised_F[i] <- sum(getUniques(ddF))
  denoised_R[i] <- sum(getUniques(ddR))
  mergers[[i]] <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = 10, maxMismatch = 1)
}
seqtab <- makeSequenceTable(mergers)
stopifnot(identical(rownames(seqtab), samples$Library))

# 5. Merged lengths BEFORE chimera removal: Table S2B ----
# Count only sequences with at least one merged read in the 31 environmental libraries.
env_reads <- colSums(seqtab[environmental, , drop = FALSE])
env_reads <- env_reads[env_reads > 0]
lengths_merged <- aggregate(list(unique_sequences = rep(1L, length(env_reads)), merged_reads = unname(env_reads)),
                            list(length_nt = nchar(names(env_reads))), sum)
lengths_merged$pct_reads <- 100 * lengths_merged$merged_reads / sum(lengths_merged$merged_reads)
write.table(lengths_merged, "DADA2_output/merged_lengths_by_nt.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# Fixed manuscript intervals; show any unexpected lengths instead of silently dropping them.
lengths_merged$interval <- cut(lengths_merged$length_nt, c(-Inf, 230, 251, 252, 253, 254, 351, Inf),
                               c("<231", "231-251", "252", "253", "254", "255-351", ">351"))
S2B <- aggregate(cbind(unique_sequences, merged_reads) ~ interval, lengths_merged, sum)
S2B$pct_reads <- 100 * S2B$merged_reads / sum(S2B$merged_reads)
write.table(S2B, "results/supplement/S2B_merged_lengths.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(S2B, row.names = FALSE)

# 6. Remove chimeras and save read tracking: Table S2C ----
# Use the table created immediately above. No readRDS or previous run is needed.
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads)
stopifnot(identical(rownames(seqtab_nochim), samples$Library))
tracking <- data.frame(ID = samples$ID, Site = samples$Site, Type = samples$Type, Date = samples$Date, Raw = filtered[, "reads.in"],
  Filtered = filtered[, "reads.out"], Denoised_R1 = denoised_F, Denoised_R2 = denoised_R,
  Merged = rowSums(seqtab), Non_chimera = rowSums(seqtab_nochim), row.names = NULL)
stopifnot(all(with(tracking, Raw >= Filtered & Filtered >= Merged & Merged >= Non_chimera)),
          sum(S2B$merged_reads) == sum(tracking$Merged[environmental]),
          all(tapply(S2A$Reads, S2A$Mate, sum) == sum(tracking$Raw[environmental])))
write.table(tracking[environmental, ], "results/supplement/S2C_environmental_read_tracking.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tracking[!environmental, ], "DADA2_output/control_read_tracking.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(colSums(tracking[environmental, c("Raw", "Filtered", "Merged", "Non_chimera")]))
saveRDS(seqtab_nochim, "DADA2_output/sequence_table_nochimera.rds")  # Optional intermediate, never required as input

# 7. Assign SILVA taxonomy ----
taxa <- assignTaxonomy(seqtab_nochim, training_file, minBoot = 80, multithread = threads, verbose = TRUE)
species <- addSpecies(taxa, species_file, verbose = TRUE)
stopifnot(setequal(rownames(taxa), colnames(seqtab_nochim)))
ranks <- intersect(c("Kingdom", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), colnames(taxa))
if (!("Kingdom" %in% ranks || "Domain" %in% ranks)) stop("No domain/kingdom column in the SILVA output.")

# 8. Write the workbook that starts 02_Sediments.R ----
# Number the columns of the non-chimeric sequence table BEFORE sorting or joining taxonomy,
# exactly as in the original DADA2 script. No ASV ID is read from an earlier run.
sequences <- colnames(seqtab_nochim)
ASV <- paste0("ASV_", seq_along(sequences))
counts <- as.data.frame(t(seqtab_nochim), check.names = FALSE)
colnames(counts) <- samples$ID
counts <- data.frame(OTU_ID = ASV, counts[, c(control_ids, samples$ID[environmental]), drop = FALSE],
                     row.names = NULL, check.names = FALSE)
taxonomy <- data.frame(ASV = ASV, Sequence = sequences, taxa[sequences, ranks, drop = FALSE],
                       row.names = NULL, check.names = FALSE)
if (!"Kingdom" %in% names(taxonomy)) names(taxonomy)[names(taxonomy) == "Domain"] <- "Kingdom"
if ("Species" %in% colnames(species)) taxonomy$Species <- species[sequences, "Species"]
# The original nucleotide-key merge sorted on sequence AFTER labelling. This also
# reproduces the row order used as input to microDecon; do not renumber after sorting.
sequence_order <- order(sequences, method = "radix")
counts <- counts[sequence_order, , drop = FALSE]
taxonomy <- taxonomy[sequence_order, , drop = FALSE]
stopifnot(identical(counts$OTU_ID, taxonomy$ASV), !anyDuplicated(taxonomy$Sequence),
          identical(names(counts)[2:4], control_ids), sum(as.matrix(counts[, -1])) == sum(tracking$Non_chimera))
write.table(taxonomy[, c("ASV", "Sequence")], "DADA2_output/ASV_to_sequence.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(as.vector(rbind(paste0(">", taxonomy$ASV), taxonomy$Sequence)), "DADA2_output/ASV_sequences.fasta")
writeLines(capture.output(sessionInfo()), "DADA2_output/sessionInfo.txt")
write_xlsx(list(counts = counts, taxonomy = taxonomy), path = "taxa.silva.xlsx")
message("Created taxa.silva.xlsx. In a fresh R session, run 02_Sediments.R from this folder.")
