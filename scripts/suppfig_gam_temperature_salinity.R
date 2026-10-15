# suppfig_gam_temperature_salinity.R
# Supplementary Figures S4 and S5: GAM responses of the MAGs most aligned with temperature and salinity.
# Run from the repository root (see README.md).

source("scripts/figure2.R")
OUT <- "results/community"
gam_stats <- list()
## ============================================================
## “Like-2C” PANELS: one faceted plot PER predictor
## (same selection rule and model as Fig. 2C)
## ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggtext); library(purrr); library(mgcv)
})

# --- pretty labels and transforms (same as before) ---
pretty_lab <- c(
  pH="pH",
  temperature_C="Temperature (°C)",
  salinity_PSU="Salinity (PSU)",
  nitrate_uM="Nitrate (µM)", nitrite_uM="Nitrite (µM)",
  phosphate_uM="Phosphate (µM)", silicate_uM="Silicate (µM)",
  chlorophyll_A_ugL="Chlorophyll <i>a</i> (µg L⁻¹)",
  C_total_ug="Total C (µg L⁻¹)", N_total_ug="Total N (µg L⁻¹)",
  C_N_ratio="C:N (molar)"
)
transform_map <- list(
  pH=identity, temperature_C=identity, salinity_PSU=identity,
  nitrate_uM=log1p, nitrite_uM=log1p, phosphate_uM=log1p, silicate_uM=log1p,
  chlorophyll_A_ugL=log1p, C_total_ug=log1p, N_total_ug=log1p, C_N_ratio=identity
)

# --- dbRDA significant predictors (term-wise) ---
anova_terms <- anova(mod_partial, permutations = 999, by = "term")
sig_preds   <- rownames(anova_terms)[anova_terms$`Pr(>F)` < 0.05]

# keep only predictors that have arrows in the biplot
bp_scores <- scores(mod_partial, display = "bp", scaling = 2) |> as.data.frame()
sig_preds <- intersect(sig_preds, rownames(bp_scores))

preds_to_do <- c("temperature_C", "salinity_PSU")
stopifnot(all(preds_to_do %in% sig_preds))

# --- species (MAG) scores for projection ranking ---
sp_scores <- scores(mod_partial, display = "species", scaling = 2) |> as.data.frame()

select_top_MAGs_along <- function(predictor, k = 6, sp = sp_scores, bp = bp_scores) {
  stopifnot(predictor %in% rownames(bp))
  v <- as.numeric(bp[predictor, c("CAP1","CAP2")])
  u <- v / sqrt(sum(v^2))
  proj <- as.matrix(sp[, c("CAP1","CAP2")]) %*% u
  rownames(sp)[order(proj, decreasing = TRUE)][seq_len(min(k, nrow(sp)))]
}

# --- taxonomy → genus-only labels (italics; fallback to higher rank) ---
tax <- as.data.frame(tax_table(ps))
safe <- function(df, nm) if (nm %in% names(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
tax_df <- tibble(MAG = rownames(tax),
                 Genus=safe(tax,"Genus"), Family=safe(tax,"Family"), Order=safe(tax,"Order"))
genus_only_label <- function(genus, family, order) {
  genus  <- ifelse(is.na(genus)  | genus=="", NA, genus)
  family <- ifelse(is.na(family) | family=="", NA, family)
  order  <- ifelse(is.na(order)  | order=="", NA, order)
  if (!is.na(genus)) sprintf("<i>%s</i>", genus) else if (!is.na(family)) family else if (!is.na(order)) order else "Unclassified"
}

# Put this near your "pretty_lab" definitions
pretty_lab_md <- c(
  pH                  = "pH",
  temperature_C       = "Temperature (°C)",
  salinity_PSU        = "Salinity (PSU)",
  nitrate_uM          = "Nitrate (µM)",
  nitrite_uM          = "Nitrite (µM)",
  phosphate_uM        = "Phosphate (µM)",
  silicate_uM         = "Silicate (µM)",
  chlorophyll_A_ugL   = "Chlorophyll *a* (µg L⁻¹)",
  C_total_ug          = "Total C (µg L⁻¹)",
  N_total_ug          = "Total N (µg L⁻¹)",
  C_N_ratio           = "C:N (molar)"
)

# ---- REPLACE your make_panel_for_predictor() with this version ----
make_panel_for_predictor <- function(predictor, k_MAG = 6) {
  # 1) select MAGs by dbRDA projection (same as 2C)
  sel_MAGs <- select_top_MAGs_along(predictor, k = k_MAG)
  
  # 2) long CLR table
  df <- as.data.frame(clr_all[, intersect(colnames(clr_all), sel_MAGs), drop = FALSE]) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(-SampleID, names_to = "MAG", values_to = "CLR_abund") |>
    dplyr::left_join(meta, by = "SampleID")
  
  # 3) genus-only labels (ensure uniqueness)
  lab_tbl <- tax_df |>
    dplyr::filter(MAG %in% sel_MAGs) |>
    dplyr::mutate(MAG_lbl = mapply(genus_only_label, Genus, Family, Order))
  dup <- duplicated(lab_tbl$MAG_lbl) | duplicated(lab_tbl$MAG_lbl, fromLast = TRUE)
  lab_tbl$MAG_lbl[dup] <- paste0(lab_tbl$MAG_lbl[dup], " (", lab_tbl$MAG[dup], ")")
  df <- df |> dplyr::left_join(lab_tbl[, c("MAG","MAG_lbl")], by = "MAG")
  df$MAG_lbl[is.na(df$MAG_lbl) | df$MAG_lbl == ""] <- df$MAG[is.na(df$MAG_lbl) | df$MAG_lbl == ""]
  
  # 4) GAMs like 2C for THIS predictor
  tfun <- transform_map[[predictor]] %||% identity
  df <- df |> dplyr::filter(is.finite(CLR_abund), is.finite(.data[[predictor]]))
  xr <- range(df[[predictor]], na.rm = TRUE); grid <- seq(xr[1], xr[2], length.out = 100)
  ref_depth  <- as.numeric(names(sort(table(df$depth_m),  decreasing = TRUE))[1])
  ref_season <- names(sort(table(df$season),             decreasing = TRUE))[1]
  
  curves <- df |>
    dplyr::group_by(MAG_lbl) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(dat) {
      dat <- dat |> dplyr::mutate(x_pred = tfun(.data[[predictor]]))
      if (sd(dat$x_pred, na.rm = TRUE) == 0) return(tibble::tibble())
      mod <- mgcv::gam(CLR_abund ~ s(x_pred, k = 5) + depth_m + season, data = dat, method = "REML")
      nd  <- tibble::tibble(x = grid, x_pred = tfun(grid),
                            depth_m = ref_depth,
                            season  = factor(ref_season, levels = levels(df$season)))
      pr  <- predict(mod, newdata = nd, se.fit = TRUE)
      st  <- summary(mod)$s.table
      gam_stats[[length(gam_stats) + 1]] <<- tibble::tibble(predictor = predictor, MAG = unique(dat$MAG),
        label = unique(dat$MAG_lbl), edf = st[1, "edf"], F = st[1, "F"], p = st[1, "p-value"])
      tibble::tibble(MAG_lbl = unique(dat$MAG_lbl),
                     x = grid,
                     fit = as.numeric(pr$fit),
                     lwr = as.numeric(pr$fit - 1.96*pr$se.fit),
                     upr = as.numeric(pr$fit + 1.96*pr$se.fit))
    })
  
  xlab_txt <- pretty_lab_md[[predictor]] %||% predictor
  
  ggplot(curves, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey85", alpha = 0.7) +
    geom_line(linewidth = 1, color = "#2C7FB8") +
    facet_wrap(~ MAG_lbl, scales = "free_y", ncol = 3) +
    labs(x = xlab_txt, y = "CLR abundance") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = ggtext::element_markdown(size = 10),
      axis.title.x = ggtext::element_markdown(),  # <-- render italics/µ/superscripts
      axis.title.y = ggtext::element_markdown()
    )
}


for (p in preds_to_do) {
  pn <- make_panel_for_predictor(p, k_MAG = 6)
  ggsave(file.path(OUT, sprintf("SuppFig_GAM_%s.png", p)), pn, width = 180, height = 120, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(OUT, sprintf("SuppFig_GAM_%s.pdf", p)), pn, width = 180, height = 120, units = "mm", device = cairo_pdf)
}
gs <- dplyr::bind_rows(gam_stats)
write.table(gs, file.path(OUT, "SuppFig_GAM_temp_sal_smooth_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(gs))
