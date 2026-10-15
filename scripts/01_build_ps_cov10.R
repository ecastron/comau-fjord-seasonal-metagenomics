# 01_build_ps_cov10.R
# Applies the detection filter: a MAG counts in a sample only where at least 10% of its genome is covered by reads (CoverM); otherwise its count is set to zero. Writes data/Comau_ps90_cov10.RDS.
# Run from the repository root (see README.md).

source("utils.R")
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(readr)})
cm <- read_tsv("data/coverm/coverm513_long.tsv.gz", show_col_types = FALSE)
ps <- readRDS("data/Comau_ps90.RDS")
otu <- as(otu_table(ps), "matrix"); tr <- taxa_are_rows(ps); if (!tr) otu <- t(otu)   # MAGs x samples
B <- matrix(NA_real_, nrow(otu), ncol(otu), dimnames = dimnames(otu))
B[cbind(cm$MAG, cm$sample)] <- cm$breadth
stopifnot(!anyNA(B))
keep <- B >= 0.10
otu_f <- otu * keep
cat("occurrences kept:", sum(otu_f > 0), "of", sum(otu > 0), "\n")
rk <- colSums(otu_f) / colSums(otu)
cat("reads kept per sample (lowest 5):\n"); print(round(sort(rk)[1:5], 3))
w <- names(which.min(rk))
top <- sort(otu[, w] * !keep[, w], decreasing = TRUE)[1:5]
cat("Removed MAGs carrying most reads in", w, ":\n")
print(data.frame(reads = top, breadth = round(B[names(top), w], 3), depth = round(cm$mean[cm$sample == w][match(names(top), cm$MAG[cm$sample == w])], 3),
                 as.data.frame(tax_table(ps))[names(top), c("Family", "Genus")]))
otu_table(ps) <- otu_table(if (tr) otu_f else t(otu_f), taxa_are_rows = tr)
stopifnot(all(sample_sums(ps) > 0), all(taxa_sums(ps) > 0))
saveRDS(ps, "data/Comau_ps90_cov10.RDS")
