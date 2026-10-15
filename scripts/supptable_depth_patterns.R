# supptable_depth_patterns.R
# Supplementary Table S4: depth dependence of each pathway's seasonal pattern (two-way ANOVA).
# Run from the repository root (see README.md).

source("utils.R")
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(tidyr); library(car)})
OUT <- "results/functional"
o  <- readRDS(file.path(OUT, "fig3_shared_objects.RDS"))
A  <- as.matrix(o$A_rel_use); H <- as.matrix(o$H_use)
P  <- A %*% t(H[, colnames(A)])                       # samples x genes, gene potentials
Prel <- P / rowSums(P)                                 # row-normalized, as for the eigengenes
sd <- as(sample_data(load_ps("data/Comau_ps90_cov10_r232.RDS")), "data.frame")[rownames(P), ]
paths <- c(o$main_pathways, "Phototrophy")
val <- sapply(paths, function(pw) {
  g <- intersect(o$pathway_dict[[pw]], colnames(P))
  if (pw == "Phototrophy") {                           # eigengene, as plotted in Fig. 3A
    m <- Prel[, g, drop = FALSE]; m <- m[, apply(m, 2, sd) > 0, drop = FALSE]
    prcomp(m, center = TRUE, scale. = TRUE)$x[, 1]
  } else rowSums(P[, g, drop = FALSE])
})
d <- data.frame(season = sd$season, depth = factor(sd$depth_m), val)
res <- lapply(paths, function(pw) {
  a <- car::Anova(lm(d[[pw]] ~ season * depth, data = d), type = 2)
  mn <- d |> group_by(depth, season) |> summarise(m = mean(.data[[pw]]), .groups = "drop")
  amp <- mn |> group_by(depth) |> summarise(r = max(m) - min(m))
  p_int <- a["season:depth", "Pr(>F)"]; p_dep <- a["depth", "Pr(>F)"]
  tibble(pathway = pw, P_season = a["season", "Pr(>F)"], P_depth = p_dep, P_interaction = p_int,
         seasonal_range_5m = amp$r[amp$depth == "5"], seasonal_range_20m = amp$r[amp$depth == "20"])
})
res <- bind_rows(res) |> mutate(P_interaction_BH = p.adjust(P_interaction, "BH"),
  ratio = seasonal_range_5m / seasonal_range_20m,
  pattern = dplyr::case_when(
    P_interaction_BH < 0.05 & ratio > 1.25  ~ "Seasonal swing larger at 5 m",
    P_interaction_BH < 0.05 & ratio < 0.8   ~ "Seasonal swing larger at 20 m",
    P_interaction_BH < 0.05                 ~ "Seasonal timing differs between depths",
    P_depth < 0.05                          ~ "Constant depth offset",
    TRUE                                    ~ "No depth effect")) |> select(-ratio)
readr::write_tsv(res, file.path(OUT, "SuppTable_pathway_depth_patterns.tsv"))
print(as.data.frame(res |> mutate(across(where(is.numeric), ~ signif(.x, 3)))))
