#!/usr/bin/env Rscript
# ─────────────────────────────────────────────
# biomarker_single.R
# Biomarker discovery for ONE microbiome sample
# Author: Alisha Ahamed (Nextflow-ready version)
# ─────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
})

# ─────────────────────────────────────────────
# CLI OPTIONS
# ─────────────────────────────────────────────
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", help = "Path to ps_rel.rds file"),
  make_option(c("-s", "--sample"), type = "character", help = "Sample ID to analyze"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory (default: ./results)", default = "results")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rds) || is.null(opt$sample)) {
  stop("❌ Missing required arguments. Use --rds and --sample.")
}

# ─────────────────────────────────────────────
# Load phyloseq object
# ─────────────────────────────────────────────
if (!file.exists(opt$rds)) stop(paste("❌ RDS not found:", opt$rds))
ps_rel <- readRDS(opt$rds)
cat("✅ Loaded phyloseq object from", opt$rds, "\n")

sample_id <- opt$sample
if (!sample_id %in% sample_names(ps_rel)) stop(paste("❌ Sample", sample_id, "not found in phyloseq object."))

# ─────────────────────────────────────────────
# Prepare output dirs
# ─────────────────────────────────────────────
dirs <- file.path(opt$output, c("biomarkers_csv", "biomarkers_plots"))
lapply(dirs, function(d) if(!dir.exists(d)) dir.create(d, recursive=TRUE))

# ─────────────────────────────────────────────
# Define biomarker groups
# ─────────────────────────────────────────────
good_bugs <- c("Bifidobacterium breve", "Eubacterium rectale", "Faecalibacterium prausnitzii",
               "Roseburia", "Subdoligranulum", "Anaerostipes", "Coprococcus", "Lactobacillus",
               "Christensenella", "Butyricimonas", "Oxalobacter", "Phascolarctobacterium")

bad_bugs <- c("Klebsiella", "Escherichia", "Proteus", "Campylobacter", "Fusobacterium",
              "Pseudomonas", "Staphylococcus", "Bilophila", "Desulfovibrio", "Clostridium")

variable_bugs <- c("Bacteroides", "Prevotella", "Alistipes", "Parabacteroides",
                   "Blautia", "Veillonella", "Sutterella", "Dorea")

biomarker_colors <- c(
  setNames(rep("green3", length(good_bugs)), good_bugs),
  setNames(rep("red3", length(bad_bugs)), bad_bugs),
  setNames(rep("orange", length(variable_bugs)), variable_bugs)
)
all_bugs <- c(good_bugs, bad_bugs, variable_bugs)

# ─────────────────────────────────────────────
# Extract abundance for this sample
# ─────────────────────────────────────────────
otu_df <- as.data.frame(otu_table(ps_rel))
if (!taxa_are_rows(ps_rel)) otu_df <- t(otu_df)
tax_df <- as.data.frame(tax_table(ps_rel))

# Combine tax + abundance for one sample
otu_tax <- cbind(tax_df, Abundance = otu_df[, sample_id])
otu_tax <- otu_tax %>%
  filter(!is.na(Genus)) %>%
  filter(sapply(tolower(Genus), function(x) any(grepl(paste(tolower(all_bugs), collapse = "|"), x))))

if (nrow(otu_tax) == 0) {
  cat("⚠️ No biomarkers found for sample:", sample_id, "\n")
  quit(save = "no", status = 0)
}

# Assign biomarker category color
otu_tax$Color <- sapply(otu_tax$Genus, function(g) {
  hit <- names(biomarker_colors)[sapply(names(biomarker_colors), function(y) grepl(tolower(y), tolower(g)))]
  if (length(hit) > 0) biomarker_colors[hit[1]] else "grey50"
})

# ─────────────────────────────────────────────
# Save results
# ─────────────────────────────────────────────
csv_path <- file.path(opt$output, "biomarkers_csv", paste0("biomarkers_", sample_id, ".csv"))
write.csv(otu_tax[, c("Genus", "Abundance", "Color")], csv_path, row.names = FALSE)
cat("✅ Saved biomarker CSV:", csv_path, "\n")

# ─────────────────────────────────────────────
# Plot and save
# ─────────────────────────────────────────────
p <- ggplot(otu_tax, aes(x = reorder(Genus, Abundance), y = Abundance, fill = Color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_identity() +
  theme_bw(base_size = 11) +
  labs(title = paste("Biomarker Abundance –", sample_id),
       x = "", y = "Abundance (relative %)")

plot_path <- file.path(opt$output, "biomarkers_plots", paste0("biomarkers_", sample_id, ".png"))
ggsave(plot_path, p, width = 8, height = 6)
cat("✅ Saved biomarker plot:", plot_path, "\n")

cat("🎯 Biomarker discovery complete for sample:", sample_id, "\n")
quit(save = "no", status = 0)
