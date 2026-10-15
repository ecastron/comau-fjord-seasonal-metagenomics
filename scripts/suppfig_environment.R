# suppfig_environment.R
# Supplementary Figures S1 (predictor correlogram) and S3 (environment by season and depth).
# Run from the repository root (see README.md).

source("utils.R")   # sets wd; load_ps(), pal_season
suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
set.seed(100)  # no stochastic steps; set for convention

OUT <- "results"
ps  <- load_ps()
sd  <- as(sample_data(ps), "data.frame")
stopifnot(nrow(sd) == 90)

# Predictors in the order of the dbRDA groups (physicochemical, nutrients, productivity)
vars <- c(temperature_C = "Temperature (°C)", salinity_PSU = "Salinity (PSU)", pH = "pH",
          nitrate_uM = "Nitrate (µM)", nitrite_uM = "Nitrite (µM)",
          phosphate_uM = "Phosphate (µM)", silicate_uM = "Silicate (µM)",
          chlorophyll_A_ugL = "Chlorophyll a (µg/L)", C_total_ug = "Particulate C (µg)",
          N_total_ug = "Particulate N (µg)", C_N_ratio = "C:N ratio")
env <- sd |> mutate(across(all_of(names(vars)), ~ suppressWarnings(as.numeric(as.character(.x)))))

## ---- Supp Fig S3: season x depth distributions ---------------------------------
long <- env |>
  select(season, depth_m, all_of(names(vars))) |>
  pivot_longer(all_of(names(vars)), names_to = "variable", values_to = "value") |>
  drop_na(value) |>
  mutate(variable = factor(vars[variable], levels = vars),
         depth = factor(paste(depth_m, "m"), levels = c("5 m", "20 m")))

# Table of means behind the ranges quoted in the Results
s3_tab <- long |>
  group_by(variable, depth, season) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), .groups = "drop")
write.table(s3_tab, file.path(OUT, "SuppFigS3_env_season_depth_means.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pS3 <- ggplot(long, aes(depth, value, fill = season)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(0.8), linewidth = 0.3) +
  geom_point(aes(color = season), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 100),
             size = 0.8, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pal_season, name = "Season") +
  scale_color_manual(values = pal_season) +
  labs(x = "Depth", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "SuppFigS3_env_season_depth.pdf"), pS3, width = 200, height = 160, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "SuppFigS3_env_season_depth.png"), pS3, width = 200, height = 160,
       units = "mm", dpi = 300, bg = "white")

## ---- Supp Fig S1: Spearman correlogram ----------------------------------------
X <- env[, names(vars)]
X <- X[complete.cases(X), ]
cat("Samples with complete predictors:", nrow(X), "\n")
k <- ncol(X)
pairs <- expand.grid(i = seq_len(k), j = seq_len(k)) |> filter(i > j)
ct <- pairs |>
  rowwise() |>
  mutate(test = list(suppressWarnings(cor.test(X[[i]], X[[j]], method = "spearman"))),
         rho = unname(test$estimate), p = test$p.value) |>
  ungroup() |>
  mutate(p_BH = p.adjust(p, "BH"),
         var1 = factor(vars[names(X)[i]], levels = vars),
         var2 = factor(vars[names(X)[j]], levels = rev(vars))) |>
  select(var1, var2, rho, p, p_BH)
write.table(ct, file.path(OUT, "SuppFigS1_env_spearman.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pS1 <- ggplot(ct, aes(var2, var1, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(p_BH < 0.05, sprintf("%.2f", rho), "")), size = 2.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
                       name = "Spearman ρ") +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT, "SuppFigS1_env_correlogram.pdf"), pS1, width = 160, height = 140, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "SuppFigS1_env_correlogram.png"), pS1, width = 160, height = 140,
       units = "mm", dpi = 300, bg = "white")

# pH partners, for the Results text
print(ct |> filter(var1 == "pH" | var2 == "pH") |> arrange(desc(abs(rho))), n = 20)
