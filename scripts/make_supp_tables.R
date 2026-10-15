# make_supp_tables.R
# Supplementary Tables S1 (samples) and S2 (MAGs).
# Run from the repository root (see README.md).

source("utils.R")
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(tidyr); library(readr)})
OUT <- "results/supplementary"
ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")
sd <- data.frame(sample_data(ps)); sd$Sample <- rownames(sd)
num <- function(v) suppressWarnings(as.numeric(as.character(v)))
s3 <- sd |> transmute(Sample, `Collection date` = format(as.Date(Collection_Date), "%Y-%m-%d"), Season = as.character(season),
  `Depth (m)` = num(depth_m), Station = line_site, Latitude = num(latitude), Longitude = num(longitude),
  `SRA run` = Run, BioSample, BioProject, `Sequenced bases` = num(Bases),
  `Temperature (°C)` = num(temperature_C), `Salinity (PSU)` = num(salinity_PSU), pH = num(pH),
  `Chlorophyll a (µg/L)` = num(chlorophyll_A_ugL), `Nitrate (µM)` = num(nitrate_uM), `Nitrite (µM)` = num(nitrite_uM),
  `Phosphate (µM)` = num(phosphate_uM), `Silicate (µM)` = num(silicate_uM), `Particulate C (µg)` = num(C_total_ug),
  `Particulate N (µg)` = num(N_total_ug), `C:N ratio` = num(C_N_ratio)) |>
  arrange(`Collection date`, Station, `Depth (m)`)
stopifnot(nrow(s3) == 90, !anyNA(s3$`SRA run`))
write_tsv(s3, file.path(OUT, "TableS1_samples.tsv"))

otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)          # MAGs x samples
rel <- sweep(otu, 2, colSums(otu), "/") * 100
grp <- paste(sd[colnames(rel), "season"], paste0(sd[colnames(rel), "depth_m"], "m"))
ab <- sapply(c("Summer 5m","Summer 20m","Fall 5m","Fall 20m","Winter 5m","Winter 20m"), function(g) rowMeans(rel[, grp == g, drop = FALSE]))
colnames(ab) <- paste0("Mean rel. abundance (%) ", sub("m$", " m", colnames(ab)))
tax <- as.data.frame(tax_table(ps))
acc <- read_tsv("data/coverm/mag677_download.tsv", col_names = c("MAG","GenBank","url"), show_col_types = FALSE)
len <- read_tsv("data/coverm/mag677_len.tsv", col_names = c("MAG","len"), show_col_types = FALSE)
o <- readRDS("results/functional/fig3_shared_objects.RDS"); H <- as.matrix(o$H_use)
pw <- sapply(names(o$pathway_dict), function(p) colSums(H[intersect(o$pathway_dict[[p]], rownames(H)), , drop = FALSE]))
pretty <- c(Respiration_ETC = "Respiration", Carbon_core = "Carbon core", B12_cobalamin = "Cobalamin (B12)", DMSP_Methyl = "DMSP and methylated compounds",
            Nitrogen_cycle = "Nitrogen cycle", Sulfur_cycle = "Sulfur cycle", Methane_redox = "Methane cycling", Phototrophy = "Phototrophy",
            Hydrogenases = "Hydrogenases", Detox_misc = "Detoxification")
colnames(pw) <- paste0("Gene hits: ", pretty[colnames(pw)])
s4 <- data.frame(MAG = rownames(otu), check.names = FALSE) |>
  left_join(acc[, c("MAG","GenBank")], by = "MAG") |> left_join(len, by = "MAG") |>
  mutate(`Genome length (Mb)` = round(len / 1e6, 3)) |> select(-len) |>
  cbind(tax[rownames(otu), c("Domain","Phylum","Class","Order","Family","Genus","Species")],
        `Samples detected (of 90)` = rowSums(otu > 0), round(ab[rownames(otu), ], 4), pw[rownames(otu), ])
stopifnot(nrow(s4) == 513, !anyNA(s4$GenBank), !anyNA(s4$`Genome length (Mb)`))
write_tsv(s4, file.path(OUT, "TableS2_MAGs.tsv"))
cat("S3:", dim(s3), "| S4:", dim(s4), "\n")
