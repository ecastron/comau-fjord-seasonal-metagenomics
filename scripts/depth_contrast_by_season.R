# depth_contrast_by_season.R
# Compares the 5 m versus 20 m community contrast across seasons (paired Bray-Curtis, within-season PERMANOVA).
# Run from the repository root (see README.md).

source("utils.R")
suppressPackageStartupMessages({library(phyloseq); library(vegan); library(dplyr); library(rstatix)})
set.seed(1)
ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")
rel <- transform_sample_counts(ps, function(x) x / sum(x))
X <- as(otu_table(rel), "matrix"); if (taxa_are_rows(rel)) X <- t(X)
bc <- as.matrix(vegdist(X, "bray"))
sd <- data.frame(sample_data(ps)); sd$id <- rownames(sd)
sd$site_date <- sub("D(05|5|20)$", "", sd$id)
pairs <- sd |> group_by(site_date) |> filter(n() == 2, all(c("5","20") %in% as.character(depth_m))) |>
  summarise(season = first(as.character(season)),
            bc = bc[id[as.character(depth_m) == "5"], id[as.character(depth_m) == "20"]], .groups = "drop") |>
  mutate(season = factor(season, levels = c("Summer","Fall","Winter")))
cat("paired 5-20 m samples per season:\n"); print(table(pairs$season))
print(pairs |> group_by(season) |> summarise(n = n(), median_BC = round(median(bc), 3), mean_BC = round(mean(bc), 3)))
print(kruskal.test(bc ~ season, data = pairs))
print(as.data.frame(dunn_test(as.data.frame(pairs), bc ~ season, p.adjust.method = "BH")))   # rstatix; FSA 0.10.1 dunnTest gave an inconsistent Kruskal-Wallis statistic
cat("\nPERMANOVA depth R2 within each season (999 perm):\n")
for (s in c("Summer","Fall","Winter")) {
  k <- sd$season == s
  a <- adonis2(as.dist(bc[sd$id[k], sd$id[k]]) ~ depth_m, data = sd[k, ], permutations = 999)
  cat(sprintf("  %-6s n = %2d  R2 = %.3f  F = %.2f  P = %.3f\n", s, sum(k), a$R2[1], a$F[1], a$`Pr(>F)`[1]))
}
write.csv(pairs, "results/depth_contrast_by_season_pairs.csv", row.names = FALSE)
