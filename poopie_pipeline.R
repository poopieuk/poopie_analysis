#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Human Stool Microbiome 16S Analysis Pipeline (Nextflow-ready)
# Author: Alisha Ahamed (refactored for CLI/Nextflow)
# Runs once for the entire dataset (all samples in --input)
# No plotting; produces JSON + CSV/RDS artifacts in subfolders
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(R.utils)
})

# -------- Helpers --------
fail <- function(msg, status = 1) {
  message("ERROR: ", msg)
  quit(save = "no", status = status)
}

info <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = " "))
}

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# -------- CLI --------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input directory containing FASTQ files"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory (root)"),
  make_option(c("-T", "--taxonomy_train"), type = "character", help = "SILVA train set fasta.gz (e.g., silva_nr_v138_train_set.fa.gz)"),
  make_option(c("-S", "--taxonomy_species"), type = "character", help = "SILVA species assignment fasta.gz (e.g., silva_species_assignment_v138.fa.gz)"),
  make_option(c("-t", "--threads"), type = "integer", default = max(1L, parallel::detectCores() - 1L), help = "Threads (default: available-1)"),
  make_option(c("--trimLeftF"), type = "integer", default = 10L, help = "Trim left bases for forward reads"),
  make_option(c("--trimLeftR"), type = "integer", default = 10L, help = "Trim left bases for reverse reads"),
  make_option(c("--truncLenF"), type = "integer", default = 260L, help = "Truncate length for forward reads"),
  make_option(c("--truncLenR"), type = "integer", default = 260L, help = "Truncate length for reverse reads"),
  make_option(c("--maxEE_F"), type = "double", default = 2, help = "maxEE for forward reads"),
  make_option(c("--maxEE_R"), type = "double", default = 5, help = "maxEE for reverse reads"),
  make_option(c("--min_abundance"), type = "integer", default = 10L, help = "Minimum total abundance per ASV to retain"),
  make_option(c("--top_genera"), type = "integer", default = 10L, help = "Top N genera per sample for JSON")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$taxonomy_train) || is.null(opt$taxonomy_species)) {
  fail("Missing required arguments. Use --input, --output, --taxonomy_train, --taxonomy_species.")
}

# Normalize paths
in_dir  <- normalizePath(opt$input, mustWork = TRUE)
out_dir <- normalizePath(opt$output, mustWork = FALSE)
safe_dir(out_dir)

# Subfolders
json_dir   <- file.path(out_dir, "json")
tables_dir <- file.path(out_dir, "tables")
rds_dir    <- file.path(out_dir, "rds")
logs_dir   <- file.path(out_dir, "logs")

invisible(lapply(c(json_dir, tables_dir, rds_dir, logs_dir), safe_dir))

# Log session info early
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "sessionInfo.txt"))

# Check taxonomy files
train_fp   <- normalizePath(opt$taxonomy_train, mustWork = TRUE)
species_fp <- normalizePath(opt$taxonomy_species, mustWork = TRUE)

# -------- Locate FASTQs --------
path <- in_dir
fnFs <- sort(list.files(path, pattern = "_R1\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2\\.fastq(\\.gz)?$", full.names = TRUE))

if (length(fnFs) == 0 || length(fnRs) == 0) fail("No FASTQ pairs found matching *_R1_001.fastq.gz and *_R2_001.fastq.gz")
if (length(fnFs) != length(fnRs)) info("WARNING: forward/reverse file counts differ; attempting to proceed by name matching.")

# Ensure pairs by basename prefix (up to _R1/_R2)
baseF <- sub("_R1\\.fastq(\\.gz)?$", "", basename(fnFs))
baseR <- sub("_R2\\.fastq(\\.gz)?$", "", basename(fnRs))
common <- intersect(baseF, baseR)
fnFs <- fnFs[baseF %in% common]
fnRs <- fnRs[baseR %in% common]

if (length(fnFs) == 0) fail("No matching forward/reverse pairs after name harmonization.")

sample.names <- common

# -------- Count reads & discard mismatched pairs by read count --------
getReads <- function(f) R.utils::countLines(f) / 4
readsF <- sapply(fnFs, getReads)
readsR <- sapply(fnRs, getReads)

mismatch <- which(readsF != readsR)
if (length(mismatch) > 0) {
  warn_tbl <- data.frame(sample = sample.names[mismatch], forward_reads = readsF[mismatch], reverse_reads = readsR[mismatch])
  write.csv(warn_tbl, file.path(logs_dir, "mismatched_read_counts.csv"), row.names = FALSE)
  info("WARNING: ", length(mismatch), " samples had unequal read counts; excluding them. See logs/mismatched_read_counts.csv")
  keep <- (readsF == readsR)
  fnFs <- fnFs[keep]; fnRs <- fnRs[keep]; sample.names <- sample.names[keep]
}

if (length(sample.names) == 0) fail("All pairs excluded due to mismatched read counts.")

# -------- Filtering & trimming --------
filt_path <- file.path(out_dir, "filtered")
safe_dir(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

info("Filtering & trimming (", length(sample.names), " samples)…")

out_filt <- dada2::filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  trimLeft = c(opt$trimLeftF, opt$trimLeftR),
  truncLen = c(opt$truncLenF, opt$truncLenR),
  maxN = 0, maxEE = c(opt$maxEE_F, opt$maxEE_R),
  truncQ = 2, rm.phix = TRUE, compress = TRUE,
  multithread = opt$threads
)

write.csv(as.data.frame(out_filt), file.path(tables_dir, "filtering_summary.csv"), row.names = TRUE)

# Remove samples with zero reads after filtering
nonzero <- (out_filt[,2] > 0 & out_filt[,4] > 0)
if (!all(nonzero)) {
  dropped <- sample.names[!nonzero]
  info("WARNING: ", sum(!nonzero), " samples dropped after filtering (zero reads)")
  write.csv(data.frame(sample = dropped), file.path(logs_dir, "dropped_after_filtering.csv"), row.names = FALSE)
  filtFs <- filtFs[nonzero]; filtRs <- filtRs[nonzero]; sample.names <- sample.names[nonzero]
}

if (length(sample.names) == 0) fail("No samples left after filtering.")

# -------- Learn errors --------
info("Learning error models…")
if (file.exists(file.path(rds_dir, "errF.rds"))) {
  errF <- readRDS(file.path(rds_dir, "errF.rds"))
  errR <- readRDS(file.path(rds_dir, "errR.rds"))
} else {
  errF <- dada2::learnErrors(filtFs, multithread = opt$threads, nbases = 1e6)
  errR <- dada2::learnErrors(filtRs, multithread = opt$threads, nbases = 1e6)
  saveRDS(errF, file.path(rds_dir, "errF.rds"))
  saveRDS(errR, file.path(rds_dir, "errR.rds"))
}



# -------- Denoise & merge --------
info("Denoising (dada) & merging pairs…")
dadaFs <- dada2::dada(filtFs, err = errF, multithread = opt$threads)
dadaRs <- dada2::dada(filtRs, err = errR, multithread = opt$threads)
mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- dada2::makeSequenceTable(mergers)
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus", multithread = opt$threads)

# -------- Filter ASVs by total abundance --------
min_abund <- as.integer(opt$min_abundance)
seqtab.filtered <- seqtab.nochim[, colSums(seqtab.nochim) > min_abund, drop = FALSE]
if (ncol(seqtab.filtered) == 0) fail("All ASVs filtered out by abundance threshold.")

# -------- Taxonomy assignment --------
info("Assigning taxonomy (SILVA)…")
tax <- dada2::assignTaxonomy(seqtab.filtered, train_fp, multithread = opt$threads)
tax <- dada2::addSpecies(tax, species_fp)

# -------- Build phyloseq object --------
otu_mat <- t(seqtab.filtered) # taxa rows, samples cols
sampledata <- data.frame(
  Sample = colnames(otu_mat),
  Group  = rep("unknown", ncol(otu_mat)),
  row.names = colnames(otu_mat),
  stringsAsFactors = FALSE
)

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
  phyloseq::tax_table(tax),
  phyloseq::sample_data(sampledata)
)

validObject(ps)

# -------- Export taxonomy table (selected columns) --------
info("Exporting taxonomy table…")
tax_df <- as.data.frame(tax)
tax_df$OTU <- rownames(tax_df)
cols <- intersect(c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(tax_df))
write.csv(tax_df[, cols, drop = FALSE], file.path(tables_dir, "taxonomy_table.csv"), row.names = FALSE)

# -------- Alpha diversity --------
info("Computing alpha diversity…")
alpha_div <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)
write.csv(alpha_div, file.path(tables_dir, "alpha_diversity.csv"), row.names = FALSE)

# -------- Genus-level relative abundance and Top N per sample --------
info("Summarizing genus-level relative abundance…")
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x) * 100)
ps_genus <- phyloseq::tax_glom(ps_rel, taxrank = "Genus")
genus_df <- phyloseq::psmelt(ps_genus) %>%
  select(Sample, Abundance, Kingdom, Phylum, Class, Order, Family, Genus)

# Clean missing taxonomy labels
for (col in c("Kingdom","Phylum","Class","Order","Family","Genus")) {
  genus_df[[col]][is.na(genus_df[[col]]) | genus_df[[col]] == ""] <- paste0("Unclassified_", col)
}

genus_df$Abundance <- round(genus_df$Abundance, 3)

# Top N genera per sample
N <- as.integer(opt$top_genera)
genus_topN <- genus_df %>%
  group_by(Sample) %>%
  arrange(desc(Abundance), .by_group = TRUE) %>%
  slice_head(n = N) %>%
  ungroup()

write.csv(genus_topN, file.path(tables_dir, sprintf("genus_top%d_per_sample.csv", N)), row.names = FALSE)

# -------- Combined JSON for Nextflow artifact --------
info("Writing combined JSON artifact…")
alpha_export <- alpha_div %>% select(Sample, Observed, Shannon, Simpson)
combined_json <- list(
  abundance = genus_topN,
  diversity = alpha_export
)

final_json_path <- file.path(json_dir, "full_microbiome_summary.json")
jsonlite::write_json(combined_json, final_json_path, pretty = TRUE, auto_unbox = TRUE)

# -------- Save important R objects (optional) --------
saveRDS(ps, file.path(rds_dir, "phyloseq_object.rds"))
saveRDS(seqtab.filtered, file.path(rds_dir, "seqtab_filtered.rds"))
saveRDS(tax, file.path(rds_dir, "taxonomy_table.rds"))
saveRDS(ps_rel, file.path(rds_dir, "ps_rel.rds"))

# -------- Done --------
info("Completed successfully. Samples processed: ", length(sample.names))
cat("OUTPUT_JSON\t", final_json_path, "\n", sep = "")  # easy to grep in logs
quit(save = "no", status = 0)
