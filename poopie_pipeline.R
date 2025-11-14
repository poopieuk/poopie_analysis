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
  make_option(c("--top_genera"), type = "integer", default = 10L, help = "Top N genera per sample for JSON"),
  make_option(c("--sample_id"), type = "character", help = "Sample identifier")

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

# -------- Locate FASTQs (generalized pattern, no hardcoding) --------
path <- in_dir

# Automatically detect all R1/R2 pairs

fnFs <- sort(list.files(in_dir, pattern="_R1\\.fastq(\\.gz)?$", full.names=TRUE))
fnRs <- sort(list.files(in_dir, pattern="_R2\\.fastq(\\.gz)?$", full.names=TRUE))

# -------- Optional: Restrict to one sample if --sample_id provided --------
if (!is.null(opt$sample_id)) {
  info("Filtering FASTQ files for sample_id: ", opt$sample_id)
  
  # Keep only files that match the sample_id string
  fnFs <- fnFs[grepl(opt$sample_id, basename(fnFs))]
  fnRs <- fnRs[grepl(opt$sample_id, basename(fnRs))]
  
  if (length(fnFs) == 0 || length(fnRs) == 0) {
    fail(paste("No FASTQ files found matching sample_id:", opt$sample_id))
  }
  
  # Restrict sample names as well
  baseF <- sub("_R1\\.fastq(\\.gz)?$", "", basename(fnFs))
  #sample.names <- gsub("_S\\d+_L\\d+", "", baseF)
  sample.names <- baseF
  
  info("Processing single sample: ", paste(sample.names, collapse = ", "))
}


if (length(fnFs) == 0 || length(fnRs) == 0)
  fail("No FASTQ pairs found matching *_R1.fastq.gz and *_R2.fastq.gz")

if (length(fnFs) != length(fnRs))
  info("WARNING: forward/reverse file counts differ; attempting to proceed by name matching.")

sample_id <- basename(in_dir)
info("Detected sample folder:", sample_id)

# Ensure pairs by basename prefix (up to _R1/_R2)
baseF <- sub("_R1_001\\.fastq(\\.gz)?$", "", basename(fnFs))
baseR <- sub("_R2_001\\.fastq(\\.gz)?$", "", basename(fnRs))
common <- intersect(baseF, baseR)
fnFs <- fnFs[baseF %in% common]
fnRs <- fnRs[baseR %in% common]

if (length(fnFs) == 0)
  fail("No matching forward/reverse pairs after name harmonization.")

# Derive clean sample names directly from FASTQ filenames
sample.names <- gsub("_R[12]\\.fastq(\\.gz)?$", "", basename(fnFs))
#sample.names <- gsub("_S\\d+_L\\d+", "", sample.names)  # remove S##/L##
sample.names <- gsub("_R[12]$", "", sample.names)
sample.names <- basename(sample.names)

info("Detected sample names: ", paste(sample.names, collapse = ", "))

# Assign sample names to files
names(fnFs) <- sample.names
names(fnRs) <- sample.names

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

# Assign sample names so DADA2 keeps them
names(fnFs) <- sample.names
names(fnRs) <- sample.names
names(filtFs) <- sample.names
names(filtRs) <- sample.names

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

# ---- Handle both single and multiple samples ----
if (ncol(out_filt) >= 4) {
  nonzero <- (out_filt[,2] > 0 & out_filt[,4] > 0)
} else if (ncol(out_filt) == 2) {
  nonzero <- (out_filt[,2] > 0)
} else {
  fail("Unexpected structure returned by filterAndTrim().")
}

if (!all(nonzero)) {
  dropped <- sample.names[!nonzero]
  info("WARNING: ", sum(!nonzero), " samples dropped after filtering (zero reads)")
  write.csv(data.frame(sample = dropped), file.path(logs_dir, "dropped_after_filtering.csv"), row.names = FALSE)
  filtFs <- filtFs[nonzero]; filtRs <- filtRs[nonzero]; sample.names <- sample.names[nonzero]
}

if (length(sample.names) == 0) fail("No samples left after filtering.")

# -------- Learn errors --------
info("Learning error models…")
errF <- dada2::learnErrors(filtFs, multithread = opt$threads, nbases = 1e6)
errR <- dada2::learnErrors(filtRs, multithread = opt$threads, nbases = 1e6)
saveRDS(errF, file.path(rds_dir, "errF.rds"))
saveRDS(errR, file.path(rds_dir, "errR.rds"))

# -------- Denoise & merge --------
info("Denoising (dada) & merging pairs…")
dadaFs <- dada2::dada(filtFs, err = errF, multithread = opt$threads)
dadaRs <- dada2::dada(filtRs, err = errR, multithread = opt$threads)
mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- dada2::makeSequenceTable(mergers)
# Maintain sample names in sequence table
# After removeBimeraDenovo
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus", multithread = opt$threads)

# Ensure rows are samples
rownames(seqtab.nochim) <- sample.names

# Filter ASVs by abundance
min_abund <- as.integer(opt$min_abundance)
seqtab.filtered <- seqtab.nochim[, colSums(seqtab.nochim) > min_abund, drop = FALSE]
if (ncol(seqtab.filtered) == 0) fail("All ASVs filtered out by abundance threshold.")



# Assign proper rownames (samples)
rownames(seqtab.nochim) <- names(filtFs)
if (is.null(rownames(seqtab.nochim))) rownames(seqtab.nochim) <- sample.names

# -------- Filter ASVs by total abundance --------
min_abund <- as.integer(opt$min_abundance)
seqtab.filtered <- seqtab.nochim[, colSums(seqtab.nochim) > min_abund, drop = FALSE]
if (ncol(seqtab.filtered) == 0) fail("All ASVs filtered out by abundance threshold.")

# ============================================================
# === TAXONOMY ASSIGNMENT (SILVA) ============================
# ============================================================

info("Assigning taxonomy (SILVA)...")

# use normalized paths (train_fp/species_fp) you already computed
taxa <- assignTaxonomy(seqtab.filtered, train_fp, multithread = opt$threads)
taxa <- tryCatch({
  addSpecies(taxa, species_fp)
}, error = function(e) {
  warning(paste("addSpecies() failed:", e$message))
  taxa
})

info("✅ Taxonomy assignment completed successfully.")

# -------- Build phyloseq object --------
info("Building Phyloseq object...")

# Confirm dimensions and rownames
info("DEBUG: seqtab.filtered dim = ", paste(dim(seqtab.filtered), collapse = " x "))
info("DEBUG: rownames(seqtab.filtered) (samples): ", paste(rownames(seqtab.filtered), collapse = ", "))

# Ensure the rownames are sample names (not sequences)
# If DADA2 lost names, reassign from our detected sample list
if (is.null(rownames(seqtab.filtered)) || all(grepl("^[ACGT]+$", rownames(seqtab.filtered)))) {
  rownames(seqtab.filtered) <- sample.names
  info("✅ Assigned sample names to seqtab.filtered rownames.")
} else {
  info("✅ Sample names already present in seqtab.filtered.")
}

# Create phyloseq object with ASVs as taxa, samples as rownames
ps <- phyloseq(
  otu_table(seqtab.filtered, taxa_are_rows = FALSE),
  tax_table(as.matrix(taxa)),
  sample_data(data.frame(SampleID = rownames(seqtab.filtered), row.names = rownames(seqtab.filtered)))
)


saveRDS(ps, file.path(rds_dir, "phyloseq_object.rds"))
info("✅ Phyloseq object successfully saved at: ", file.path(rds_dir, "ps_rel.rds"))
validObject(ps)


# -------- Export taxonomy table (selected columns) --------
info("Exporting taxonomy table…")
tax_df <- as.data.frame(taxa)
tax_df$OTU <- rownames(tax_df)
cols <- intersect(c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(tax_df))
write.csv(tax_df[, cols, drop = FALSE], file.path(tables_dir, "taxonomy_table.csv"), row.names = FALSE)


# -------- Alpha diversity --------
info("Computing alpha diversity…")
alpha_div <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)

# Normalize sample IDs in alpha_div
alpha_div$Sample <- gsub("\\.fastq.*$", "", alpha_div$Sample)
alpha_div$Sample <- gsub("_R[12]$", "", alpha_div$Sample)
#alpha_div$Sample <- gsub("_S\\d+_L\\d+", "", alpha_div$Sample)
write.csv(alpha_div, file.path(tables_dir, "alpha_diversity.csv"), row.names = FALSE)

# -------- Genus-level relative abundance and Top N per sample --------
info("Summarizing genus-level relative abundance…")
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x) * 100)
ps_genus <- phyloseq::tax_glom(ps_rel, taxrank = "Genus")
genus_df <- phyloseq::psmelt(ps_genus) %>%
  dplyr::select(Sample, Abundance, Kingdom, Phylum, Class, Order, Family, Genus)

# Normalize Sample IDs
genus_df$Sample <- gsub("\\.fastq.*$", "", genus_df$Sample)
genus_df$Sample <- gsub("_R[12]$", "", genus_df$Sample)
#genus_df$Sample <- gsub("_S\\d+_L\\d+", "", genus_df$Sample)

# Replace missing taxonomy with "Unclassified_*"
for (col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
  genus_df[[col]][is.na(genus_df[[col]]) | genus_df[[col]] == ""] <- paste0("Unclassified_", col)
}
genus_df$Abundance <- round(genus_df$Abundance, 3)

# ---- Compute Top N genera safely ----
N <- as.integer(opt$top_genera)
if (!is.finite(N) || N <= 0) N <- 10

if (nrow(genus_df) > 0) {
  genus_topN <- genus_df %>%
    dplyr::group_by(Sample) %>%
    dplyr::arrange(dplyr::desc(Abundance), .by_group = TRUE) %>%
    dplyr::slice_head(n = N) %>%
    dplyr::ungroup()
} else {
  info("WARNING: genus_df is empty; creating empty genus_topN table.")
  genus_topN <- genus_df
}

# Save the Top-N CSV
write.csv(
  genus_topN,
  file.path(tables_dir, sprintf("genus_top%d_per_sample.csv", N)),
  row.names = FALSE
)

# -------- Combined JSON for Nextflow artifact --------
info("Writing combined JSON artifact…")

# Normalize sample names for JSON
genus_topN$Sample <- gsub("\\.fastq.*$", "", genus_topN$Sample)
genus_topN$Sample <- gsub("_R[12]$", "", genus_topN$Sample)
#genus_topN$Sample <- gsub("_S\\d+_L\\d+", "", genus_topN$Sample)
#genus_topN$Sample <- gsub("_L\\d+$", "", genus_topN$Sample)

metadata <- data.frame(Sample = unique(genus_df$Sample))
alpha_export <- alpha_div %>%
  dplyr::select(Sample, Observed, Shannon, Simpson)

combined_json <- list(
  abundance = genus_topN,
  diversity = alpha_export,
  metadata  = metadata
)

final_json_path <- file.path(json_dir, "full_microbiome_summary.json")
jsonlite::write_json(combined_json, final_json_path, pretty = TRUE, auto_unbox = TRUE)

combined_json <- list(
  abundance = genus_topN,
  diversity = alpha_export,
  metadata = metadata
)


final_json_path <- file.path(json_dir, "full_microbiome_summary.json")
# Normalize Sample names in JSON export
genus_topN$Sample <- gsub("\\.fastq.*$", "", genus_topN$Sample)
genus_topN$Sample <- gsub("_R[12]$", "", genus_topN$Sample)
#genus_topN$Sample <- gsub("_S\\d+_L\\d+", "", genus_topN$Sample)
#genus_topN$Sample <- gsub("_L\\d+$", "", genus_topN$Sample)
jsonlite::write_json(combined_json, final_json_path, pretty = TRUE, auto_unbox = TRUE)

# -------- Save important R objects (optional) --------
saveRDS(ps, file.path(rds_dir, "phyloseq_object.rds"))
saveRDS(seqtab.filtered, file.path(rds_dir, "seqtab_filtered.rds"))
saveRDS(taxa, file.path(rds_dir, "taxonomy_table.rds"))
saveRDS(ps_rel, file.path(rds_dir, "ps_rel.rds"))

# -------- Done --------
info("Completed successfully. Samples processed: ", length(sample.names))
cat("OUTPUT_JSON\t", final_json_path, "\n", sep = "")  # easy to grep in logs
quit(save = "no", status = 0)
