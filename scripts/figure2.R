# figure2.R
# Figure 2: partial dbRDA, unique variance fractions, GAM responses to pH, and distance-decay (MRM).
# Run from the repository root (see README.md).

source("utils.R")
OUT <- "results/community"
set.seed(1)

cran_pkgs <- c(
  "vegan", "mgcv", "dplyr", "tidyr", "tibble",
  "ggplot2", "ggrepel", "scales", "ggtext", "patchwork", "cowplot",
  "ecodist", "geosphere", "lubridate"
)
bioc_pkgs <- c("phyloseq")

to_install_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if (length(to_install_cran))
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")

to_install_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[,"Package"]]
if (length(to_install_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(to_install_bioc, ask = FALSE)
}
suppressPackageStartupMessages(invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE)))

## ---- Read data ----
ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")

## ---- Season colors (from utils.R: Wong 2011 palette matching Fig. 1) & depth shapes ----
season_palette <- pal_season          # Summer=#E69F00, Fall=#D55E00, Winter=#56B4E9
depth_shapes   <- c("5 m" = 16, "20 m" = 17)

## ======================================================================
## SHARED: community (relative), metadata, predictors, VIF-filtered X
## ======================================================================
# Community (relative abundances, samples as rows)
Y_counts <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) Y_counts <- t(Y_counts)
Y_counts <- as.data.frame(Y_counts)
rs <- rowSums(Y_counts); rs[rs == 0] <- 1
Y_rel <- sweep(Y_counts, 1, rs, "/")

# Metadata
meta <- as(sample_data(ps), "data.frame") %>% rownames_to_column("SampleID")
meta$season  <- factor(meta$season, levels = c("Fall","Summer","Winter"))
meta$depth_m <- suppressWarnings(as.numeric(as.character(meta$depth_m)))

# Predictor groups
vars_nutrients    <- c("nitrate_uM","nitrite_uM","phosphate_uM","silicate_uM")
vars_physicals    <- c("salinity_PSU","temperature_C","pH")
vars_productivity <- c("chlorophyll_A_ugL","C_total_ug","N_total_ug","C_N_ratio")
vars_all <- unique(c(
  intersect(vars_nutrients, colnames(meta)),
  intersect(vars_physicals, colnames(meta)),
  intersect(vars_productivity, colnames(meta))
))

# Align & scale predictors
numify <- function(df) as.data.frame(lapply(df, function(v) suppressWarnings(as.numeric(v))))
X_raw <- meta[, vars_all, drop = FALSE]
X_num <- numify(X_raw); rownames(X_num) <- meta$SampleID
common <- intersect(rownames(Y_rel), rownames(X_num))
Y_rel  <- Y_rel[common,, drop=FALSE]
X_num  <- X_num[common,, drop=FALSE]
meta   <- meta[match(common, meta$SampleID),, drop=FALSE]
cc <- complete.cases(X_num, meta$depth_m, meta$season)
Y_rel <- Y_rel[cc,, drop=FALSE]; X_num <- X_num[cc,, drop=FALSE]; meta <- meta[cc,, drop=FALSE]
X <- as.data.frame(scale(X_num))

X_vif <- drop_high_vif_vegan(X, Y_rel, thr = 10)  # defined in utils.R

## =========================
## Panel 2A: partial dbRDA
## =========================
env_partial  <- cbind(X_vif, depth_m = meta$depth_m, season = meta$season)
form_partial <- as.formula(
  paste("Y_rel ~", paste(colnames(X_vif), collapse = " + "),
        "+ Condition(depth_m) + Condition(season)")
)
mod_partial <- capscale(form_partial, data = env_partial, distance = "bray")

# Scores and labels
sc_sites <- scores(mod_partial, display="sites", scaling=2) %>% as.data.frame()
sc_bp    <- scores(mod_partial, display="bp",    scaling=2) %>% as.data.frame()
sc_sites$SampleID <- rownames(sc_sites)
plot_df <- sc_sites %>% left_join(meta, by="SampleID")
eig <- mod_partial$CCA$eig
lab_x <- sprintf("dbRDA1 (%.1f%%)", 100*eig[1]/sum(eig))
lab_y <- sprintf("dbRDA2 (%.1f%%)", 100*eig[2]/sum(eig))
arrow_mul <- vegan:::ordiArrowMul(mod_partial, display="bp")
sc_bp$arrow_x <- sc_bp$CAP1 * arrow_mul; sc_bp$arrow_y <- sc_bp$CAP2 * arrow_mul; sc_bp$var <- rownames(sc_bp)
plot_df$depth_std <- factor(ifelse(plot_df$depth_m==5,"5 m","20 m"), levels=c("5 m","20 m"))

p2A <- ggplot(plot_df, aes(x=CAP1, y=CAP2)) +
  geom_point(aes(color=season, shape=depth_std), size=2.5, alpha=0.9) +
  scale_color_manual(values=season_palette, drop=FALSE) +
  scale_shape_manual(values=depth_shapes, drop=FALSE) +
  geom_segment(data=sc_bp, aes(x=0,y=0,xend=arrow_x,yend=arrow_y),
               arrow=arrow(length=unit(0.25,"cm")), linewidth=0.5) +
  ggrepel::geom_text_repel(data=sc_bp, aes(x=arrow_x,y=arrow_y,label=var),
                           size=3, max.overlaps=Inf) +
  labs(x=lab_x, y=lab_y,
       title="Figure 2A. Partial dbRDA of MAG communities",
       subtitle=sprintf("R²=%.3f, adj.R²=%.3f; perm.p=%.3g",
                        RsquareAdj(mod_partial)$r.squared,
                        RsquareAdj(mod_partial)$adj.r.squared,
                        anova(mod_partial, permutations=999)$`Pr(>F)`[1])) +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), legend.title=element_blank(),
        plot.title.position="plot")

## =========================
## Panel 2B: unique adj R²
## =========================
fit_unique <- function(Y, X_main, X_condlist, extra_cond=NULL) {
  dat_list <- c(list(X_main), X_condlist, if (!is.null(extra_cond)) list(extra_cond) else NULL)
  dat <- do.call(cbind, dat_list); colnames(dat) <- make.unique(colnames(dat))
  k <- ncol(X_main)
  main_nms <- colnames(dat)[seq_len(k)]; cond_nms <- colnames(dat)[-(seq_len(k))]
  form_txt <- paste("Y ~", paste(main_nms, collapse=" + "),
                    if (length(cond_nms)) paste("+ Condition(", paste(cond_nms, collapse=" + "), ")", sep="") else "")
  mod <- capscale(as.formula(form_txt), data=dat, distance="bray")
  list(R2=RsquareAdj(mod), anova=anova(mod, permutations=999))
}
extra_cond <- data.frame(depth_m=meta$depth_m, season=meta$season)
X_N <- X_vif[, intersect(colnames(X_vif), vars_nutrients),    drop=FALSE]
X_P <- X_vif[, intersect(colnames(X_vif), vars_physicals),    drop=FALSE]
X_R <- X_vif[, intersect(colnames(X_vif), vars_productivity), drop=FALSE]
uN_c <- fit_unique(Y_rel, X_N, list(X_P, X_R), extra_cond=extra_cond)
uP_c <- fit_unique(Y_rel, X_P, list(X_N, X_R), extra_cond=extra_cond)
uR_c <- fit_unique(Y_rel, X_R, list(X_N, X_P), extra_cond=extra_cond)

summ_tbl <- tibble(
  Group = c("Nutrients","Physicochemical","Productivity"),
  R2_adj = c(uN_c$R2$adj.r.squared, uP_c$R2$adj.r.squared, uR_c$R2$adj.r.squared),
  pval   = c(uN_c$anova$`Pr(>F)`[1], uP_c$anova$`Pr(>F)`[1], uR_c$anova$`Pr(>F)`[1])
) %>%
  mutate(sig = case_when(pval<=0.001~"***", pval<=0.01~"**", pval<=0.05~"*",
                         pval<=0.1~".", TRUE~"ns"))

p2B <- ggplot(summ_tbl, aes(x=Group, y=R2_adj*100, fill=Group)) +
  geom_col(width=0.65) +
  geom_text(aes(label=sprintf("adj.R²=%.3f\np=%.3g %s", R2_adj, pval, sig)),
            vjust=-0.35, size=3.1) +
  scale_y_continuous(labels=percent_format(scale=1, accuracy=0.1),
                     expand=expansion(mult=c(0.02,0.15))) +
  scale_fill_manual(values=c("Nutrients"="#9ecae1","Physicochemical"="#a1d99b","Productivity"="#fdae6b")) +
  labs(x=NULL, y="Unique explained variance (adj.R², %)",
       title="Figure 2B. Unique fractions after conditioning on depth and season") +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), legend.position="none", plot.title.position="plot")

## ==================================================
## Panel 2C: GAMs vs pH (taxonomy + fixed CLR baseline)
## ==================================================
# CLR ONCE over full community (fixed baseline; subcomposition-safe)
clr_transform <- function(x, pseudocount = 1e-6) {
  x <- as.matrix(x)
  x <- sweep(x, 1, rowSums(x), "/")
  log((x + pseudocount) / exp(rowMeans(log(x + pseudocount))))
}
Y_all <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) Y_all <- t(Y_all)
clr_all <- clr_transform(Y_all)

# Select MAGs aligned with pH vector in partial dbRDA (top k)
sp <- scores(mod_partial, display="species", scaling=2) %>% as.data.frame()
bp <- scores(mod_partial, display="bp",      scaling=2) %>% as.data.frame()
stopifnot("pH" %in% rownames(bp))
v <- as.numeric(bp["pH", c("CAP1","CAP2")]); u <- v / sqrt(sum(v^2))
proj <- as.matrix(sp[, c("CAP1","CAP2")]) %*% u
n_MAGs_pH <- 6  # set to 12 if you prefer 12 panels
sel_MAGs <- rownames(sp)[order(proj, decreasing = TRUE)][seq_len(min(n_MAGs_pH, nrow(sp)))]

# --- Genus-only facet labels (italics), higher-rank fallbacks plain ---
tax <- as.data.frame(tax_table(ps))
get_safe <- function(df, nm) if (nm %in% colnames(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
tax_df <- tibble::tibble(
  MAG   = rownames(tax),
  Genus = get_safe(tax, "Genus"),
  Family= get_safe(tax, "Family"),
  Order = get_safe(tax, "Order")
)

genus_only_label <- function(genus, family, order) {
  genus  <- ifelse(is.na(genus)  | genus  %in% c("", "NA"), NA, genus)
  family <- ifelse(is.na(family) | family %in% c("", "NA"), NA, family)
  order  <- ifelse(is.na(order)  | order  %in% c("", "NA"), NA, order)
  if (!is.na(genus)) sprintf("<i>%s</i>", genus)
  else if (!is.na(family)) family
  else if (!is.na(order))  order
  else                     "Unclassified"
}

tax_sel <- tax_df %>%
  dplyr::filter(MAG %in% sel_MAGs) %>%
  dplyr::mutate(MAG_lbl = mapply(genus_only_label, Genus, Family, Order))

# ensure uniqueness of labels; append MAG id if needed
dup <- duplicated(tax_sel$MAG_lbl) | duplicated(tax_sel$MAG_lbl, fromLast = TRUE)
tax_sel$MAG_lbl[dup] <- paste0(tax_sel$MAG_lbl[dup], " (", tax_sel$MAG[dup], ")")

# Build long table using CLR AFTER selection (but computed globally)
get_clr_abund_full <- function(clr_all, taxa, meta) {
  as.data.frame(clr_all[, intersect(colnames(clr_all), taxa), drop = FALSE]) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "MAG", values_to = "CLR_abund") %>%
    left_join(meta, by = "SampleID")
}

# join labels to df_clr (keep the rest of your 2C code as-is)
df_clr <- get_clr_abund_full(clr_all, sel_MAGs, meta) %>%
  dplyr::left_join(tax_sel[, c("MAG","MAG_lbl")], by = "MAG") %>%
  dplyr::filter(!is.na(pH), !is.na(CLR_abund))

# Fit per-MAG GAMs and predict smooths at modal depth/season
rng  <- range(df_clr$pH, na.rm=TRUE); grid <- seq(rng[1], rng[2], length.out=100)
ref_depth  <- as.numeric(names(sort(table(df_clr$depth_m),  decreasing=TRUE))[1])
ref_season <- names(sort(table(df_clr$season), decreasing=TRUE))[1]
mods <- df_clr %>% group_by(MAG_lbl) %>%
  group_map(~ gam(CLR_abund ~ s(pH, k = 5) + depth_m + season, data = .x, method = "REML"))
names(mods) <- dplyr::group_keys(dplyr::group_by(df_clr, MAG_lbl))$MAG_lbl  # match group_map() order
pred_list <- lapply(names(mods), function(lbl) {
  nd <- data.frame(pH = grid, depth_m = ref_depth,
                   season = factor(ref_season, levels = levels(df_clr$season)))
  pr <- predict(mods[[lbl]], newdata = nd, se.fit = TRUE)
  tibble(MAG_lbl = lbl, pH = grid, fit = pr$fit,
         lwr = pr$fit - 1.96*pr$se.fit, upr = pr$fit + 1.96*pr$se.fit)
})
smooths_pH <- bind_rows(pred_list)

p2C <- ggplot(smooths_pH, aes(x = pH, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey85", alpha = 0.7) +
  geom_line(linewidth = 1, color = "#3182BD") +
  facet_wrap(~ MAG_lbl, scales = "free_y", ncol = ifelse(n_MAGs_pH >= 12, 4, 3)) +
  labs(x = "pH", y = "CLR abundance",
       title = "Figure 2C. Taxon-level responses to pH",
       subtitle = "GAM smooths (95% CI); model: CLR ~ s(pH) + depth + season") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(size = 9),
        plot.title.position = "plot")

## ==============================================
## Panel 2D: Distance-decay MRM (5 m & 20 m, DOY)
## ==============================================
# build_station_level / build_distances / run_mrm defined in utils.R
# env_candidates / depth_var / lat_var etc. defined in utils.R
depths_to_show <- c(5, 20); depth_tol <- 0.6; use_time <- TRUE; nperm <- 9999

# Run per depth, tidy, and plot
pred_map <- c(B = "Geographic distance", C = "Environmental distance", D = "Seasonal time distance")
run_one_depth <- function(dv) {
  sta   <- build_station_level(ps, depth_value = dv, depth_tolerance = depth_tol)
  dists <- build_distances(sta$X_sta, sta$md_sta, use_time = use_time)
  res   <- run_mrm(A = dists$A, B = dists$B, C = dists$C, D = dists$D, nperm = nperm)
  cf <- as.data.frame(res$mrm$coef) %>% tibble::rownames_to_column("term")
  val_col <- setdiff(colnames(cf), c("term","pval")); stopifnot(length(val_col) == 1)
  list(
    tidy = cf %>%
      filter(term != "(Intercept)", term %in% res$used_predictors) %>%
      transmute(term, estimate = .data[[val_col]], pval, depth = paste0(dv, " m")),
    R2   = res$mrm$r.squared, adjR2 = res$mrm$adj.r.squared
  )
}
res_list <- lapply(depths_to_show, run_one_depth); names(res_list) <- paste0(depths_to_show, " m")
coef_df <- bind_rows(lapply(res_list, `[[`, "tidy")) %>%
  mutate(predictor = recode(term, !!!pred_map),
         depth = factor(depth, levels = paste0(depths_to_show, " m")),
         sig = case_when(pval<=0.001~"***", pval<=0.01~"**", pval<=0.05~"*", pval<=0.1~".", TRUE~"ns"))

p2D <- ggplot(coef_df, aes(x = predictor, y = estimate, fill = depth)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("p=%.3g %s", pval, sig)),
            position = position_dodge(width = 0.7), vjust = -0.45, size = 3.0) +
  scale_fill_manual(values = c("5 m" = "#3182BD", "20 m" = "#9ecae1")) +
  labs(x = NULL, y = "MRM coefficient (beta; predictors standardized)",
       title = "Figure 2D. Multiple regression on distance matrices (MRM)",
       subtitle = "Response: Bray-Curtis community distance; predictors: geographic, environmental, seasonal time") +
  coord_cartesian(ylim = range(c(0, coef_df$estimate)) * c(1, 1.25)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1),
        plot.title.position = "plot")

## =========================
## Assembly (cowplot, 2×2)
## - legends merged (season + depth)
## - legends on right (closer)
## =========================
suppressPackageStartupMessages({ library(cowplot); library(ggplot2) })

# 1) Remove titles/subtitles
p2A_clean <- p2A + labs(title = NULL, subtitle = NULL)
p2B_clean <- p2B + labs(title = NULL, subtitle = NULL)
p2C_clean <- p2C + labs(title = NULL, subtitle = NULL)
p2D_clean <- p2D + labs(title = NULL, subtitle = NULL)

# 2) Extract unified legend (season color + depth shape from 2A)
leg_unified <- get_legend(
  p2A_clean +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.title = element_blank())
)

# 3) Remove legends inside panels
p2A_noleg <- p2A_clean + theme(legend.position = "none")
p2B_noleg <- p2B_clean + theme(legend.position = "none")
p2C_noleg <- p2C_clean + theme(legend.position = "none")
p2D_noleg <- p2D_clean + theme(legend.position = "none") +
  guides(fill = "none")  # remove fill legend from 2D

# 4) Add per-panel padding
pad_theme <- theme(plot.margin = margin(t = 10, r = 12, b = 10, l = 12))
p2A_pad <- p2A_noleg + pad_theme
p2B_pad <- p2B_noleg + pad_theme
p2C_pad <- p2C_noleg + pad_theme
p2D_pad <- p2D_noleg + pad_theme

# 5) Build rows (A|B, C|D)
row_top    <- plot_grid(p2A_pad, p2B_pad, ncol = 2,
                        labels = c("A","B"), label_size = 12, label_fontface = "bold")
row_bottom <- plot_grid(p2C_pad, p2D_pad, ncol = 2,
                        labels = c("C","D"), label_size = 12, label_fontface = "bold")

main_grid <- plot_grid(row_top, row_bottom, ncol = 1, rel_heights = c(1, 1))

# 6) Outer margin for full grid
main_grid_outer <- ggdraw(main_grid) + theme(plot.margin = margin(8, 8, 8, 8))

# 7) Final assembly: grid + merged legend
fig2_symm <- plot_grid(
  main_grid_outer,
  ggdraw(leg_unified) + theme(plot.margin = margin(2, 2, 2, 2)),
  ncol = 2, rel_widths = c(1, 0.15), align = "h"
)

print(fig2_symm)
ggsave(file.path(OUT, "Figure2.pdf"), fig2_symm, width = 300, height = 180, units = "mm")
ggsave(file.path(OUT, "Figure2.png"), fig2_symm, width = 300, height = 180, units = "mm", dpi = 300, bg = "white")









# 1) What survived VIF?
colnames(X_vif)

# 2) Term-wise significance (partial dbRDA)
anova_terms <- anova(mod_partial, permutations = 999, by = "term")
anova_terms[order(anova_terms$`Pr(>F)`), ]

# 3) Which candidates were dropped at which stage?
candidates <- c("pH","temperature_C","salinity_PSU","nitrate_uM","nitrite_uM",
                "phosphate_uM","silicate_uM","chlorophyll_A_ugL",
                "C_total_ug","N_total_ug","C_N_ratio")

data.frame(
  predictor = candidates,
  in_VIF = candidates %in% colnames(X_vif),
  term_p = anova_terms[candidates, "Pr(>F)"],
  sig_in_dbRDA = candidates %in% rownames(anova_terms)[anova_terms$`Pr(>F)` < 0.05]
)


# ===========================
# Extra exports
# ===========================
set.seed(1)
glob <- anova(mod_partial, permutations = 999)
r2   <- RsquareAdj(mod_partial)
readr::write_tsv(tibble::tibble(R2 = r2$r.squared, adjR2 = r2$adj.r.squared,
                                F = glob$F[1], p = glob$`Pr(>F)`[1],
                                predictors_after_VIF = paste(colnames(X_vif), collapse = ";"),
                                n_samples = nrow(Y_rel)),
                 file.path(OUT, "Figure2A_partial_dbRDA_global.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(anova_terms), "term"),
                 file.path(OUT, "Figure2A_partial_dbRDA_terms.tsv"))
readr::write_tsv(tibble::rownames_to_column(sc_bp, "var_row"), file.path(OUT, "Figure2A_biplot_scores.tsv"))
readr::write_tsv(tibble::tibble(group = c("Nutrients","Physicochemical","Productivity"),
                                vars = c(paste(colnames(X_N), collapse=";"), paste(colnames(X_P), collapse=";"),
                                         paste(colnames(X_R), collapse=";")),
                                adjR2 = c(uN_c$R2$adj.r.squared, uP_c$R2$adj.r.squared, uR_c$R2$adj.r.squared),
                                F = c(uN_c$anova$F[1], uP_c$anova$F[1], uR_c$anova$F[1]),
                                p = c(uN_c$anova$`Pr(>F)`[1], uP_c$anova$`Pr(>F)`[1], uR_c$anova$`Pr(>F)`[1])),
                 file.path(OUT, "Figure2B_unique_fractions.tsv"))
# GAM smooth terms, refit per MAG ID so labels cannot be misaligned
gam_tab <- dplyr::bind_rows(lapply(sel_MAGs, function(m) {
  d  <- dplyr::filter(df_clr, MAG == m)
  st <- summary(gam(CLR_abund ~ s(pH, k = 5) + depth_m + season, data = d, method = "REML"))$s.table
  tibble::tibble(MAG = m, label = unique(d$MAG_lbl), edf = st[1, "edf"], F = st[1, "F"], p = st[1, "p-value"])
}))
readr::write_tsv(gam_tab, file.path(OUT, "Figure2C_GAM_smooth_terms.tsv"))
readr::write_tsv(tibble::tibble(mods_name_order = names(mods),
                                group_map_key_order = dplyr::group_keys(dplyr::group_by(df_clr, MAG_lbl))$MAG_lbl),
                 file.path(OUT, "Figure2C_label_order_check.tsv"))
readr::write_tsv(dplyr::bind_rows(lapply(names(res_list), function(d)
                   dplyr::mutate(res_list[[d]]$tidy, R2 = res_list[[d]]$R2[1], R2_p = res_list[[d]]$R2[2]))),
                 file.path(OUT, "Figure2D_MRM_coefficients.tsv"))
