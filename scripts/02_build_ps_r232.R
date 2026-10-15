# 02_build_ps_r232.R
# Replaces the taxonomy with the GTDB-Tk v2.7.2 / GTDB R232 classification of the 513 MAGs. Writes data/Comau_ps90_cov10_r232.RDS, the object used by all analyses.
# Run from the repository root (see README.md).

suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(tidyr); library(readr)})
D <- "data/gtdbtk_r232"
rd <- function(f) read_tsv(file.path(D, f), show_col_types = FALSE, col_types = cols(.default = "c"))
new <- bind_rows(rd("gtdbtk.bac120.summary.tsv"), rd("gtdbtk.ar53.summary.tsv")) |>
  transmute(MAG = sub("[.]fa$", "", user_genome), classification) |>
  separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";") |>
  mutate(across(Domain:Species, ~ sub("^[a-z]__", "", .x)), across(Domain:Species, ~ ifelse(.x == "", NA, .x)))
ps <- readRDS("data/Comau_ps90_cov10.RDS")
tt <- as.matrix(tibble::column_to_rownames(as.data.frame(new), "MAG"))[taxa_names(ps), colnames(tax_table(ps))]
stopifnot(!anyNA(tt[, "Phylum"]), identical(rownames(tt), taxa_names(ps)))
tax_table(ps) <- tax_table(tt)
saveRDS(ps, "data/Comau_ps90_cov10_r232.RDS")
cat("saved; phyla:", paste(sort(unique(tt[, "Phylum"])), collapse = ", "), "\n")
