#!/usr/bin/env Rscript
# ─────────────────────────────────────────────
# generate_summary_single.R
# Human Stool Microbiome Summary (Single Sample)
# Author: Alisha Ahamed (optimized for Nextflow)
# ─────────────────────────────────────────────

suppressMessages({
  library(phyloseq)
  library(vegan)
  library(scales)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

# ─────────────────────────────────────────────
# Command-line arguments
# ─────────────────────────────────────────────
option_list <- list(
  make_option(c("-p", "--phyloseq_rds"), type="character", help="Input ps_rel.rds file"),
  make_option(c("-s", "--sample_id"), type="character", help="Sample ID to summarize"),
  make_option(c("-o", "--output_dir"), type="character", default="results/reports", help="Output directory (default: results/reports)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$phyloseq_rds) || is.null(opt$sample_id)) {
  stop("❌ Required arguments missing: --phyloseq_rds and --sample_id")
}

ps_rds <- normalizePath(opt$phyloseq_rds, mustWork = TRUE)
sample_id <- opt$sample_id
out_dir <- opt$output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("✅ Using ps_rel.rds:", ps_rds, "\n")
cat("📍 Processing sample:", sample_id, "\n")

# ─────────────────────────────────────────────
# Load ps_rel object
# ─────────────────────────────────────────────
ps_rel <- readRDS(ps_rds)
if (!(sample_id %in% sample_names(ps_rel))) {
  stop(paste0("❌ Sample '", sample_id, "' not found in ps_rel.rds"))
}

# Subset to the single sample
ps_rel <- prune_samples(sample_id, ps_rel)

# ─────────────────────────────────────────────
# Genus-level abundance table
# ─────────────────────────────────────────────
cat("📊 Generating genus-level abundance table...\n")
genus_table <- tax_glom(ps_rel, taxrank = "Genus")
genus_abundances <- as.data.frame(t(otu_table(genus_table)))

taxa_info <- as.data.frame(tax_table(genus_table))
colnames(genus_abundances) <- taxa_info[rownames(tax_table(genus_table)), "Genus"]
colnames(genus_abundances)[is.na(colnames(genus_abundances))] <- "unclassified"
colnames(genus_abundances) <- make.names(colnames(genus_abundances), unique = TRUE)

genus_abundances$SampleID <- sample_id
write.csv(genus_abundances, file.path(out_dir, "genus_abundance.csv"), row.names = FALSE)
cat("✅ Saved genus_abundance.csv\n")

# ─────────────────────────────────────────────
# Calculate diversity
# ─────────────────────────────────────────────
alpha_div <- estimate_richness(ps_rel, measures = "Shannon")
alpha_div$Sample <- rownames(alpha_div)
diversity_index <- alpha_div$Shannon[1]
if (is.na(diversity_index)) diversity_index <- 0

# ─────────────────────────────────────────────
# Define genus categories
# ─────────────────────────────────────────────
good_bugs <- tolower(c(
  "bifidobacterium","eubacterium","butyricicoccus","anaerostipes",
  "coprococcus","subdoligranulum","phascolarctobacterium","barnesiella",
  "christensenellaceae_r_7_group","lactobacillus","streptococcus",
  "butyricimonas","roseburia","faecalibacterium","marvinbryantia"
))

bad_bugs <- tolower(c(
  "klebsiella","escherichia","proteus","campylobacter","fusobacterium",
  "pseudomonas","staphylococcus","bilophila","desulfovibrio",
  "clostridium_sensu_stricto_1"
))

variable_bugs <- tolower(c(
  "bacteroides","alistipes","prevotella","parabacteroides","blautia",
  "veillonella","sutterella","dorea","lachnospira","odoribacter"
))

abund <- genus_abundances
colnames(abund) <- tolower(colnames(abund))

get_sum <- function(taxa) {
  hits <- grep(paste(taxa, collapse="|"), colnames(abund), value=TRUE)
  if (length(hits) == 0) return(0)
  sum(as.numeric(abund[1, hits, drop=FALSE]), na.rm=TRUE)
}

good_total <- get_sum(good_bugs)
bad_total  <- get_sum(bad_bugs)
variable_total <- get_sum(variable_bugs)
total <- good_total + bad_total + variable_total + 1e-6

good_ratio <- good_total / total
bad_ratio <- bad_total / total
variable_ratio <- variable_total / total

# ─────────────────────────────────────────────
# Compute scores (scaled 0–10)
# ─────────────────────────────────────────────
gut_score        <- rescale((good_ratio * diversity_index) / (bad_ratio + 1e-6), to=c(0,10))
mental_score     <- rescale((good_ratio + variable_ratio) * diversity_index, to=c(0,10))
skin_score       <- rescale((good_ratio / (bad_ratio + 1e-6)) * diversity_index, to=c(0,10))
nutrition_score  <- rescale(good_ratio * diversity_index, to=c(0,10))
pathobiont_score <- rescale(1 / (bad_ratio + 1e-6), to=c(0,10))

summary_df <- data.frame(
  SampleID = sample_id,
  GutScore = round(gut_score, 2),
  MentalScore = round(mental_score, 2),
  SkinScore = round(skin_score, 2),
  NutritionScore = round(nutrition_score, 2),
  PathobiontScore = round(pathobiont_score, 2),
  DiversityIndex = round(diversity_index, 3),
  GoodBugLoad = round(good_total, 5),
  BadBugLoad = round(bad_total, 5),
  VariableBugLoad = round(variable_total, 5)
)

write.csv(summary_df, file.path(out_dir, "microbiome_summary.csv"), row.names = FALSE)
cat("✅ Saved microbiome_summary.csv\n")

# ─────────────────────────────────────────────
# Visualizations
# ─────────────────────────────────────────────
p1 <- ggplot(summary_df, aes(x=GoodBugLoad, y=BadBugLoad)) +
  geom_point(color="#4CAF50", size=3) +
  theme_minimal(base_size=12) +
  labs(title="Good vs Bad Bug Load", x="Good Bugs", y="Bad Bugs")

p2 <- ggplot(summary_df, aes(x=DiversityIndex, y=GutScore)) +
  geom_point(color="#009688", size=3) +
  theme_minimal(base_size=12) +
  labs(title="Gut Score vs Diversity", x="Shannon Diversity", y="Gut Health Score")

ggsave(file.path(out_dir, "good_vs_bad_load.png"), p1, width=5, height=4)
ggsave(file.path(out_dir, "gut_vs_diversity.png"), p2, width=5, height=4)
cat("📊 Saved visualizations to", out_dir, "\n")

# ─────────────────────────────────────────────
# Done
# ─────────────────────────────────────────────
cat("✅ Completed summary for sample:", sample_id, "\n")
quit(save = "no", status = 0)
