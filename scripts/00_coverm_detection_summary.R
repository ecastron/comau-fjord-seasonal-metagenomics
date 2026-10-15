# 00_coverm_detection_summary.R
# Summarizes per-sample CoverM coverage of the 513 MAGs and compares detection rules (observed versus expected breadth, minimum covered fraction).
# Run from the repository root (see README.md).

source("utils.R")
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(tidyr); library(purrr); library(readr)})
D <- "data/coverm"
cm <- list.files(file.path(D, "per_sample"), full.names = TRUE) |>
  map_dfr(function(f) {
    x <- read_tsv(f, show_col_types = FALSE)
    names(x) <- c("MAG", "count", "breadth", "mean", "tmean", "relab")
    x |> filter(MAG != "unmapped") |> mutate(sample = sub("\\.tsv$", "", basename(f)))
  }) |>
  mutate(exp_breadth = 1 - exp(-0.883 * mean), ratio = ifelse(exp_breadth > 0, breadth / exp_breadth, NA))
ps <- load_ps()
otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
stopifnot(setequal(unique(cm$sample), colnames(otu)), setequal(unique(cm$MAG), rownames(otu)))
cm$ps_count <- otu[cbind(cm$MAG, cm$sample)]
cat("pairs:", nrow(cm), "| CoverM count > 0:", sum(cm$count > 0), "| ps count > 0:", sum(cm$ps_count > 0), "\n")
cat("Spearman CoverM vs ps counts:", round(cor(cm$count, cm$ps_count, method = "spearman"), 3), "\n")
cat("median per-sample mapped fraction (100 - unmapped%) not shown; breadth quantiles:\n"); print(quantile(cm$breadth, c(.1,.25,.5,.75,.9)))
cat("mean depth quantiles:\n"); print(quantile(cm$mean, c(.1,.25,.5,.75,.9)))
cat("ratio quantiles (count > 0):\n"); print(quantile(cm$ratio[cm$count > 0], c(.05,.1,.25,.5,.75), na.rm = TRUE))
grid <- expand_grid(ratio_min = c(0.3, 0.5, 0.7), breadth_min = c(0, 0.01, 0.05, 0.1))
res <- pmap_dfr(grid, function(ratio_min, breadth_min) {
  det <- cm$count > 0 & !is.na(cm$ratio) & cm$ratio >= ratio_min & cm$breadth >= breadth_min
  per_s <- tapply(det, cm$sample, sum); per_m <- tapply(det, cm$MAG, sum)
  tibble(ratio_min, breadth_min, pairs_detected = sum(det), pct_pairs = round(100 * mean(det), 1),
         MAGs_per_sample_median = median(per_s), MAGs_per_sample_min = min(per_s),
         MAGs_never_detected = sum(per_m == 0), MAGs_in_all_samples = sum(per_m == 90),
         reads_retained_pct = round(100 * sum(cm$count[det]) / sum(cm$count), 1))
})
print(as.data.frame(res))
write_tsv(res, file.path(D, "coverm513_detection_grid.tsv"))
write_tsv(cm, file.path(D, "coverm513_long.tsv.gz"))
