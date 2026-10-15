# figure5.R
# Figure 5: phylogenetic (SES-MPD), functional (SES-FDis), and trait-specific (SES-CWSD) community structure against null models.
# Run from the repository root (see README.md).

source("utils.R")
OUT <- "results/functional"


## ===============================
## Figure 5A: SES-MPD by Season & Depth (abundance-weighted; parallel)
## Inputs: ps (rooted phylogeny; sample_data contains 'season' and 'depth_m')
## Outputs: two side-by-side panels (richness null; frequency null)
## ===============================
ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")
# Load shared objects from Figure 3 (A_rel_use, H_use, pathway_dict, main_pathways)
if (!exists("A_rel_use")) {
  fig3 <- readRDS(file.path(OUT, "fig3_shared_objects.RDS"))
  list2env(fig3, envir = .GlobalEnv)
}
cran_pkgs <- c(
  "picante", "dplyr", "ape", "tidyr", "tibble",
  "ggplot2", "patchwork"
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

## -------- User-tunable parameters --------------------------------------------
RNGkind("L'Ecuyer-CMRG")           # parallel-safe RNG: mclapply streams are reproducible only with this generator
set.seed(123)                  # reproducible permutations & worker RNG streams
runs              <- 999       # ↓ to 199/499 for iteration; 999 for final
min_prev_samples  <- 2         # drop taxa present in < this many samples
min_max_rel_abund <- 1e-4      # drop taxa whose MAX relative abundance < this
min_sample_rich   <- 3         # (QC) flag samples with richness < this
ncores_requested  <- NA        # NA → auto-detect cores; set integer to force

## -------- Preconditions -------------------------------------------------------
stopifnot(inherits(ps, "phyloseq"))
if (is.null(phy_tree(ps))) stop("ps lacks a phylogeny.")
meta <- as(sample_data(ps), "data.frame")
if (!all(c("season","depth_m") %in% names(meta)))
  stop("sample_data(ps) must have 'season' and 'depth_m'.")

## -------- Community: counts → relative; rare-taxa filter ---------------------
otu <- as(otu_table(ps), "matrix")
comm_counts <- if (taxa_are_rows(ps)) t(otu) else otu     # samples × taxa
# Relative abundance per sample
row_tot  <- rowSums(comm_counts, na.rm = TRUE)
comm_rel <- sweep(comm_counts, 1, ifelse(row_tot == 0, 1, row_tot), "/")
comm_rel[is.na(comm_rel)] <- 0

# Rarity metrics (prevalence & max RA)
tax_prev  <- colSums(comm_rel > 0, na.rm = TRUE)
tax_maxRA <- apply(comm_rel, 2, max, na.rm = TRUE)
keep_taxa <- (tax_prev >= min_prev_samples) | (tax_maxRA >= min_max_rel_abund)
comm_rel  <- comm_rel[, keep_taxa, drop = TRUE]
if (ncol(comm_rel) < 5) warning("After filtering, very few taxa remain; consider relaxing thresholds.")

## -------- Align tree and community once (match.phylo.comm) -------------------
phy <- phy_tree(ps); if (!inherits(phy, "phylo")) phy <- ape::as.phylo(phy)
mp <- picante::match.phylo.comm(phy, comm_rel)   # prunes & reorders both
phy_aln  <- mp$phy
comm_aln <- mp$comm

## -------- Aligned metadata (no duplicate 'sample' columns) -------------------
meta_tbl <- meta %>% tibble::rownames_to_column(".sample_row")
meta_aln <- tibble::tibble(sample = rownames(comm_aln)) %>%
  dplyr::left_join(meta_tbl, by = c("sample" = ".sample_row"))
# QC: per-sample richness
rich_tbl <- tibble::tibble(sample = rownames(comm_aln),
                           richness = rowSums(comm_aln > 0, na.rm = TRUE))
meta_aln <- meta_aln %>% dplyr::left_join(rich_tbl, by = "sample")
# Factor ordering
if (!is.factor(meta_aln$season))
  meta_aln$season <- factor(meta_aln$season, levels = c("Summer","Fall","Winter"))
meta_aln$depth_m <- factor(as.character(meta_aln$depth_m), levels = c("5","20"))

## -------- Distances (compute once) -------------------------------------------
Dphy <- ape::cophenetic.phylo(phy_aln)

## -------- Parallel SES-MPD runner (Windows-safe) -----------------------------
split_indices <- function(n, k) {
  if (k <= 1 || n <= 1) return(list(seq_len(n)))
  sz  <- ceiling(n / k)
  unname(split(seq_len(n), ceiling(seq_along(seq_len(n)) / sz)))
}
run_ses_block <- function(row_idx, null_model, runs, comm, D, aw = TRUE) {
  samp_block <- comm[row_idx, , drop = FALSE]
  res <- picante::ses.mpd(
    samp = samp_block, dis = D,
    null.model = null_model, runs = runs,
    abundance.weighted = aw
  )
  res$sample <- rownames(res); res$null <- null_model
  tibble::as_tibble(res)
}
run_ses_parallel <- function(comm, D, null_model, runs, abundance.weighted = TRUE, ncores = NA_integer_) {
  n <- nrow(comm)
  cores <- if (is.na(ncores)) max(1L, min(parallel::detectCores(TRUE) - 1L, 8L)) else max(1L, as.integer(ncores))
  cores <- min(cores, n)
  blocks <- split_indices(n, cores)
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterSetRNGStream(cl, iseed = .Random.seed[1])
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(picante)))
    parallel::clusterExport(cl, c("comm","D","runs","null_model","run_ses_block"), envir = environment())
    res_list <- parallel::parLapply(cl, blocks, function(idx) run_ses_block(idx, null_model, runs, comm, D, aw = TRUE))
  } else {
    res_list <- parallel::mclapply(blocks, function(idx) run_ses_block(idx, null_model, runs, comm, D, aw = TRUE), mc.cores = cores)
  }
  dplyr::bind_rows(res_list)
}

ncores_use <- if (is.na(ncores_requested)) NA_integer_ else as.integer(ncores_requested)

ses_rich <- run_ses_parallel(comm = comm_aln, D = Dphy, null_model = "richness",
                             runs = runs, abundance.weighted = TRUE, ncores = ncores_use) %>%
  dplyr::left_join(meta_aln, by = "sample")

ses_freq <- run_ses_parallel(comm = comm_aln, D = Dphy, null_model = "frequency",
                             runs = runs, abundance.weighted = TRUE, ncores = ncores_use) %>%
  dplyr::left_join(meta_aln, by = "sample")

## -------- Figure 4B: plotting + depth significance (no ggsignif) -------------

# Season order
ses_levels  <- c("Summer","Fall","Winter")

# Coerce factors exactly once
ses_rich$season  <- factor(ses_rich$season, levels = ses_levels)
ses_freq$season  <- factor(ses_freq$season, levels = ses_levels)
ses_rich$depth_m <- factor(ses_rich$depth_m, levels = c("5","20"), labels = c("5 m","20 m"))
ses_freq$depth_m <- factor(ses_freq$depth_m, levels = c("5","20"), labels = c("5 m","20 m"))

# Shared y-limits (safe)
y_all <- c(ses_rich$mpd.obs.z, ses_freq$mpd.obs.z)
y_lim <- if (length(y_all) == 0 || all(is.na(y_all))) 2 else max(2, ceiling(max(abs(y_all), na.rm = TRUE)))
ylims <- c(-y_lim, y_lim)

# Geoms & scales
pd <- position_dodge(width = 0.6)
zero_line <- geom_hline(yintercept = 0, linetype = "solid", color = "#BFBFBF", linewidth = 0.4)

box_sd <- geom_boxplot(
  aes(fill = season, alpha = depth_m, group = interaction(season, depth_m)),
  position = pd, width = 0.50, linewidth = 0.4,
  outlier.shape = NA, color = "grey35"
)
pt_jd <- geom_point(
  aes(shape = depth_m, color = season),
  position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.6),
  alpha = 0.85, size = 1.6, stroke = 0
)

depth_alpha <- scale_alpha_manual(name = "Depth (m)", values = c("5 m" = 0.90, "20 m" = 0.55))
depth_shape <- scale_shape_manual(name = "Depth (m)", values = c("5 m" = 16, "20 m" = 17))

theme_fig <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "white", colour = "grey70"),
    plot.margin      = margin(4, 6, 4, 6),
    legend.position  = "bottom"
  )

common_scales <- list(
  scale_color_manual(values = pal_season, guide = "none"),
  scale_fill_manual(values = pal_season, guide = "none"),
  depth_alpha, depth_shape
)

# Build the two panels
p_rich <- ggplot(ses_rich, aes(x = season, y = mpd.obs.z, group = interaction(season, depth_m))) +
  zero_line + box_sd + pt_jd + common_scales +
  coord_cartesian(ylim = ylims) +
  labs(x = NULL, y = "SES-MPD (abundance-weighted)", title = "") +
  theme_fig

p_freq <- ggplot(ses_freq, aes(x = season, y = mpd.obs.z, group = interaction(season, depth_m))) +
  zero_line + box_sd + pt_jd + common_scales +
  coord_cartesian(ylim = ylims) +
  labs(x = NULL, y = "SES-MPD (abundance-weighted)", title = "") +
  theme_fig

## -------- Significance: 5 m vs 20 m within each Season (Welch t-test) --------
p_stars <- function(p) ifelse(is.na(p), "n/a",
                              ifelse(p < 0.001, "***",
                                     ifelse(p < 0.01,  "**",
                                            ifelse(p < 0.05,  "*",  "ns"))))

depth_tests_annot <- function(df, ylims, dodge_span = 0.28) {
  rng  <- diff(ylims); bump <- 0.06 * rng
  df %>%
    dplyr::filter(!is.na(mpd.obs.z), !is.na(season), !is.na(depth_m)) %>%
    dplyr::group_by(season) %>%
    dplyr::summarise(
      ok = (length(unique(depth_m)) == 2) && all(table(depth_m) >= 2),
      p  = if (ok) stats::t.test(mpd.obs.z ~ depth_m)$p.value else NA_real_,
      y  = max(mpd.obs.z, na.rm = TRUE) + bump,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(p)) %>%
    dplyr::mutate(
      label = p_stars(p),
      xmin  = as.numeric(season) - dodge_span,
      xmax  = as.numeric(season) + dodge_span
    )
}

ses_rich_plot <- dplyr::filter(ses_rich, !is.na(mpd.obs.z), !is.na(season), !is.na(depth_m))
ses_freq_plot <- dplyr::filter(ses_freq, !is.na(mpd.obs.z), !is.na(season), !is.na(depth_m))

ann_rich <- depth_tests_annot(ses_rich_plot, ylims)
ann_freq <- depth_tests_annot(ses_freq_plot, ylims)

# Ensure brackets fit
y_top_needed <- max(
  c(ylims[2],
    if (nrow(ann_rich)) max(ann_rich$y, na.rm = TRUE) else -Inf,
    if (nrow(ann_freq)) max(ann_freq$y, na.rm = TRUE) else -Inf),
  na.rm = TRUE
)
if (is.finite(y_top_needed) && y_top_needed > ylims[2]) {
  ylims <- c(ylims[1], y_top_needed + 0.03 * diff(ylims))
  p_rich <- p_rich + coord_cartesian(ylim = ylims)
  p_freq <- p_freq + coord_cartesian(ylim = ylims)
}

# Draw brackets (manual; version-proof)
draw_sig <- function(p, ann, ylims, tick_frac = 0.03, text_bump = 0.015) {
  if (is.null(ann) || nrow(ann) == 0) return(p)
  rng  <- diff(ylims)
  tick <- tick_frac * rng
  tbmp <- text_bump * rng
  p +
    geom_segment(
      data = ann,
      aes(x = xmin, xend = xmax, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.4, color = "grey20"
    ) +
    geom_segment(
      data = ann,
      aes(x = xmin, xend = xmin, y = y, yend = y - tick),
      inherit.aes = FALSE, linewidth = 0.4, color = "grey20"
    ) +
    geom_segment(
      data = ann,
      aes(x = xmax, xend = xmax, y = y, yend = y - tick),
      inherit.aes = FALSE, linewidth = 0.4, color = "grey20"
    ) +
    geom_text(
      data = ann,
      aes(x = (xmin + xmax)/2, y = y + tbmp, label = label),
      inherit.aes = FALSE, size = 3.2, color = "black", vjust = 0
    )
}

p_rich <- draw_sig(p_rich, ann_rich, ylims)
p_freq <- draw_sig(p_freq, ann_freq, ylims)

## -------- Assemble -----------------------------------------------------------
p_4B <- p_rich | p_freq
print(p_4B)






## ===============================
## Figure 5B: SES-FDis (abundance-weighted) by Season × Depth
## ===============================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(ggplot2) })

# -- Settings / palettes (reuse Fig. 1D) --
ses_levels <- c("Summer","Fall","Winter")
pal_season <- pal_season[ses_levels]  # enforce order

# -- Metadata (season & depth; ordered as Fig. 1D) --
meta <- as(sample_data(ps), "data.frame") %>%
  tibble::rownames_to_column("sample") %>%
  select(sample, season, depth_m) %>%
  mutate(
    season  = factor(season, levels = ses_levels),
    depth_m = factor(as.character(depth_m), levels = c("5","20"))
  )

# -- Community matrix aligned to metadata (samples × MAGs; rel. abundance) --
comm <- A_rel_use[meta$sample, , drop = FALSE]
stopifnot(all(rownames(comm) == meta$sample))

# -- MAG × pathway trait matrix from H_use (6 pathway axes incl. Phototrophy) --
if (!exists("main_pathways")) {
  main_pathways <- c("Methane_redox","Carbon_core","Sulfur_cycle","B12_cobalamin","Nitrogen_cycle")
}
fig3_paths <- unique(c(main_pathways, "Phototrophy"))
gene_cols_H <- rownames(H_use)
avail <- function(g) intersect(g, gene_cols_H)

mag_ids <- colnames(H_use)
mag_path_hits <- lapply(fig3_paths, function(pw){
  gg <- avail(pathway_dict[[pw]])
  if (length(gg) == 0) return(rep(0, length(mag_ids)))
  colSums(H_use[gg, , drop = FALSE], na.rm = TRUE)
})
mag_path_hits <- do.call(cbind, mag_path_hits)
colnames(mag_path_hits) <- fig3_paths
rownames(mag_path_hits) <- mag_ids

# Keep MAGs present in both
common_mag <- intersect(colnames(comm), rownames(mag_path_hits))
comm       <- comm[, common_mag, drop = FALSE]
traits_raw <- mag_path_hits[common_mag, , drop = FALSE]

# Standardize traits (log1p then column-scale)
traits_use <- scale(log1p(traits_raw))
traits_use[is.na(traits_use)] <- 0
keep_traits <- apply(traits_use, 2, function(v) sd(v) > 0)
traits_use  <- traits_use[, keep_traits, drop = FALSE]
stopifnot(ncol(traits_use) >= 2)

# ---------- Fast observed + SES for FDis ----------
fast_fd_obs <- function(comm, traits_use) {
  stopifnot(nrow(traits_use) == ncol(comm))
  rs <- rowSums(comm); rs[!is.finite(rs) | rs == 0] <- 1
  W <- comm / rs
  T <- as.matrix(traits_use)
  C <- W %*% T
  N <- nrow(W); P <- ncol(T)
  FDis <- numeric(N)
  for (s in seq_len(N)) {
    cs  <- as.numeric(C[s, ])
    dev <- sweep(T, 2, cs, FUN = "-")
    di  <- sqrt(rowSums(dev^2))
    FDis[s] <- sum(W[s, ] * di)
  }
  names(FDis) <- rownames(W)
  list(FDis = FDis)
}

ses_fdis <- function(comm, traits_use, null = c("richness","frequency"), R = 999, ncores = 1) {
  null <- match.arg(null)
  rs <- rowSums(comm); rs[!is.finite(rs) | rs == 0] <- 1
  W <- comm / rs
  T0 <- as.matrix(traits_use)
  N  <- nrow(W); S <- nrow(T0)
  
  # observed
  Fobs <- fast_fd_obs(comm, traits_use)$FDis
  
  # one permutation
  perm_once <- function() {
    if (null == "richness") {
      Tperm <- T0
      for (samp in seq_len(N)) {
        idx <- which(W[samp, ] > 0)
        if (length(idx) >= 2)
          Tperm[idx, ] <- T0[sample(idx, length(idx), replace = FALSE), ]
      }
      T_use <- Tperm
    } else {
      T_use <- T0[sample.int(S, S), , drop = FALSE]
    }
    C <- W %*% T_use
    Fp <- numeric(N)
    for (s in seq_len(N)) {
      cs  <- as.numeric(C[s, ])
      dev <- sweep(T_use, 2, cs, FUN = "-")
      di  <- sqrt(rowSums(dev^2))
      Fp[s] <- sum(W[s, ] * di)
    }
    Fp
  }
  
  permlist <- if (ncores > 1 && Sys.info()[["sysname"]] != "Windows") {
    parallel::mclapply(seq_len(R), function(i) perm_once(), mc.cores = ncores)
  } else replicate(R, perm_once(), simplify = FALSE)
  
  Fmat <- do.call(cbind, permlist)      # N × R
  muF  <- rowMeans(Fmat)
  sdF  <- apply(Fmat, 1, sd)
  SES  <- (Fobs - muF) / ifelse(sdF > 0, sdF, 1)
  names(SES) <- names(Fobs)
  SES
}

# Compute SES-FDis (choose null; here richness)
SES_FDis <- ses_fdis(comm, traits_use, null = "richness", R = 999, ncores = 1)

# Tidy for plotting
dfC <- tibble(sample = names(SES_FDis), SES_FDis = as.numeric(SES_FDis)) %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(season), !is.na(depth_m), is.finite(SES_FDis))

# y-limits (symmetric)
ymax  <- max(2, ceiling(max(abs(dfC$SES_FDis), na.rm = TRUE)))
ylims <- c(-ymax, ymax)

# compact depth contrast brackets (5 m vs 20 m per season)
p_stars <- function(p) ifelse(p < 0.001,"***", ifelse(p < 0.01,"**", ifelse(p < 0.05,"*","ns")))
annC <- dfC %>%
  group_by(season) %>%
  summarise(
    ok = (length(unique(depth_m)) == 2) && all(table(depth_m) >= 2),
    p  = if (ok) stats::t.test(SES_FDis ~ depth_m)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  filter(is.finite(p)) %>%
  mutate(
    label = p_stars(p),
    xmin  = as.numeric(season) - 0.28,
    xmax  = as.numeric(season) + 0.28,
    y     = tapply(dfC$SES_FDis, dfC$season, max, na.rm = TRUE)[as.character(season)] + 0.06*diff(ylims)
  )
# keep limits symmetric but leave headroom so every bracket label stays inside the panel
ylims <- { L <- max(abs(ylims), max(annC$y) + 0.10*diff(ylims)); c(-L, L) }

add_brackets <- function(p, ann, ylims, tick = 0.03*diff(ylims), bump = 0.015*diff(ylims)) {
  if (nrow(ann) == 0) return(p)
  p +
    geom_segment(data = ann, aes(x = xmin, xend = xmax, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4, color = "grey20") +
    geom_segment(data = ann, aes(x = xmin, xend = xmin, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4, color = "grey20") +
    geom_segment(data = ann, aes(x = xmax, xend = xmax, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4, color = "grey20") +
    geom_text(data = ann, aes(x = (xmin+xmax)/2, y = y + bump, label = label),
              inherit.aes = FALSE, size = 3.0, vjust = 0)
}

# Plot
pd <- position_dodge(width = 0.6)
p_4C <- ggplot(dfC, aes(season, SES_FDis, group = interaction(season, depth_m))) +
  geom_boxplot(aes(fill = season, alpha = depth_m), position = pd, width = 0.50,
               linewidth = 0.4, outlier.shape = NA, color = "grey35") +
  geom_point(aes(color = season, shape = depth_m),
             position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.6),
             alpha = 0.85, size = 1.6, stroke = 0) +
  scale_color_manual(values = pal_season, guide = "none", drop = FALSE) +
  scale_fill_manual(values = pal_season, guide = "none", drop = FALSE) +
  scale_alpha_manual(name = "Depth (m)", values = c("5"=0.90,"20"=0.55),
                     breaks = c("5","20"), labels = c("5 m","20 m")) +
  scale_shape_manual(name = "Depth (m)", values = c("5"=16,"20"=17),
                     breaks = c("5","20"), labels = c("5 m","20 m")) +
  coord_cartesian(ylim = ylims) +
  labs(x = NULL, y = "SES-FDis (abundance-weighted)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white", colour = "grey70"))

p_4C <- add_brackets(p_4C, annC, ylims)
print(p_4C)





## ===============================
## Figure 5C: SES-CWSD (per-trait dispersion) by Season × Depth
## ===============================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(patchwork) })

# Reuse 'comm', 'traits_use', 'meta', 'ses_levels', 'pal_season' from 4C above

# Extend SES to also return per-trait dispersion (CWSD)
ses_fdis_trait <- function(comm, traits_use, null = c("richness","frequency"), R = 999, ncores = 1) {
  null <- match.arg(null)
  rs <- rowSums(comm); rs[!is.finite(rs) | rs == 0] <- 1
  W <- comm / rs
  T0 <- as.matrix(traits_use)
  N  <- nrow(W); S <- nrow(T0); P <- ncol(T0)
  
  # observed centroids & CWSD (traits × samples)
  C_obs <- W %*% T0
  CWSD_obs <- matrix(NA_real_, P, N)
  for (s in seq_len(N)) {
    cs  <- as.numeric(C_obs[s, ])
    dev <- sweep(T0, 2, cs, FUN = "-")
    CWSD_obs[, s] <- sqrt(colSums(W[s, ] * (dev^2)))
  }
  rownames(CWSD_obs) <- colnames(traits_use); colnames(CWSD_obs) <- rownames(comm)
  
  # one permutation → CWSDp (P × N)
  perm_once <- function() {
    if (null == "richness") {
      Tperm <- T0
      for (samp in seq_len(N)) {
        idx <- which(W[samp, ] > 0)
        if (length(idx) >= 2)
          Tperm[idx, ] <- T0[sample(idx, length(idx), replace = FALSE), ]
      }
      T_use <- Tperm
    } else {
      T_use <- T0[sample.int(S, S), , drop = FALSE]
    }
    C <- W %*% T_use
    CWSDp <- matrix(NA_real_, P, N)
    for (s in seq_len(N)) {
      cs  <- as.numeric(C[s, ])
      dev <- sweep(T_use, 2, cs, FUN = "-")
      CWSDp[, s] <- sqrt(colSums(W[s, ] * (dev^2)))
    }
    CWSDp
  }
  
  permlist <- if (ncores > 1 && Sys.info()[["sysname"]] != "Windows") {
    parallel::mclapply(seq_len(R), function(i) perm_once(), mc.cores = ncores)
  } else replicate(R, perm_once(), simplify = FALSE)
  
  # streaming mean/sd across permutations (Welford)
  P <- nrow(CWSD_obs); N <- ncol(CWSD_obs)
  m  <- 0; M1 <- matrix(0, P, N); M2 <- matrix(0, P, N)
  for (jj in seq_along(permlist)) {
    X <- permlist[[jj]]
    m  <- m + 1
    d  <- X - M1
    M1 <- M1 + d / m
    M2 <- M2 + d * (X - M1)
  }
  mu <- M1
  sd <- if (m > 1) sqrt(M2 / (m - 1)) else 0
  
  SES_CWSD <- (CWSD_obs - mu) / ifelse(sd > 0, sd, 1)
  rownames(SES_CWSD) <- rownames(CWSD_obs); colnames(SES_CWSD) <- colnames(CWSD_obs)
  SES_CWSD
}

# Compute SES-CWSD (choose null; here richness)
SES_CWSD <- ses_fdis_trait(comm, traits_use, null = "richness", R = 999, ncores = 1)

# Tidy long
panel_traits <- colnames(traits_use)  # will include Phototrophy + 5 main
dfD <- as.data.frame(SES_CWSD) %>%
  tibble::rownames_to_column("trait") %>%
  pivot_longer(-trait, names_to = "sample", values_to = "SES_CWSD") %>%
  left_join(meta, by = "sample") %>%
  filter(trait %in% panel_traits, !is.na(season), !is.na(depth_m), is.finite(SES_CWSD)) %>%
  mutate(
    trait  = factor(trait, levels = panel_traits),
    season = factor(season, levels = ses_levels)
  )

# Brackets (5 m vs 20 m within each season × trait)
p_stars <- function(p) ifelse(p < 0.001,"***", ifelse(p < 0.01,"**", ifelse(p < 0.05,"*","ns")))
make_ann <- function(df) {
  rng <- diff(range(df$SES_CWSD, na.rm = TRUE)); if (!is.finite(rng) || rng == 0) rng <- 1
  df %>%
    group_by(trait, season) %>%
    summarise(
      ok = (length(unique(depth_m)) == 2) && all(table(depth_m) >= 2),
      p  = suppressWarnings(if (ok) stats::t.test(SES_CWSD ~ depth_m)$p.value else NA_real_),
      top = max(SES_CWSD, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(p)) %>%
    mutate(
      label = p_stars(p),
      xmin  = as.numeric(season) - 0.28,
      xmax  = as.numeric(season) + 0.28,
      y     = top + 0.08 * rng
    )
}
add_brackets <- function(p, ann) {
  if (nrow(ann) == 0) return(p)
  p +
    geom_segment(data = ann, aes(x = xmin, xend = xmax, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.35, color = "grey20") +
    geom_text(data = ann, aes(x = (xmin+xmax)/2, y = y, label = label),
              inherit.aes = FALSE, size = 3.0, vjust = 0)
}
annD <- make_ann(dfD)

# y-lims (symmetric)
ymaxD  <- max(2, ceiling(max(abs(dfD$SES_CWSD), na.rm = TRUE)))
ylimsD <- c(-ymaxD, ymaxD)

# Plot
pd <- position_dodge(width = 0.6)
p_4D <- ggplot(dfD, aes(season, SES_CWSD, group = interaction(season, depth_m))) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_boxplot(aes(fill = season, alpha = depth_m), position = pd, width = 0.50,
               linewidth = 0.4, outlier.shape = NA, color = "grey35") +
  geom_point(aes(color = season, shape = depth_m),
             position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.6),
             alpha = 0.85, size = 1.6, stroke = 0) +
  facet_wrap(~ trait, nrow = 2, scales = "fixed", labeller = as_labeller(c(DMSP_Methyl = "DMSP", Carbon_core = "Carbon core", Sulfur_cycle = "Sulfur cycle", B12_cobalamin = "Cobalamin (B12)", Nitrogen_cycle = "Nitrogen cycle", Phototrophy = "Phototrophy", Methane_redox = "Methane cycling"))) +
  scale_color_manual(values = pal_season, guide = "none", drop = FALSE) +
  scale_fill_manual(values = pal_season, guide = "none", drop = FALSE) +
  scale_alpha_manual(name = "Depth (m)", values = c("5"=0.90,"20"=0.55),
                     breaks = c("5","20"), labels = c("5 m","20 m")) +
  scale_shape_manual(name = "Depth (m)", values = c("5"=16,"20"=17),
                     breaks = c("5","20"), labels = c("5 m","20 m")) +
  coord_cartesian(ylim = ylimsD) +
  labs(x = NULL, y = "SES-CWSD") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white", colour = "grey70"))

p_4D <- add_brackets(p_4D, annD)
print(p_4D)









## ===============================
## Figure 5: Assembly (force legend bottom)
## ===============================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

# Strip titles
p5A <- p_freq + labs(title = NULL)
p5B <- p_4C   + labs(title = NULL)
p5C <- p_4D + labs(title = NULL)

# Shared layout
design <- "
AB
CC
"

fig5 <- p5A + p5B + p5C +
  plot_layout(design = design, guides = "collect", heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")   # force legend placement

# Print and save
print(fig5)
ggsave(file.path(OUT, "Figure5.pdf"),  fig5, width = 180, height = 230, units = "mm", device = cairo_pdf)
#ggsave("Figure5_final.tiff", fig5, width = 180, height = 160, units = "mm", dpi = 600, compression = "lzw")

ggsave(file.path(OUT, "Figure5.png"), fig5, width = 180, height = 230, units = "mm", dpi = 300, bg = "white")

# ===========================
# Extra exports
# ===========================
vs_zero <- function(df, value, ...) {
  df |> dplyr::filter(is.finite(.data[[value]])) |>
    dplyr::group_by(...) |>
    dplyr::summarise(n = dplyr::n(), median = median(.data[[value]]), mean = mean(.data[[value]]),
                     t_p = if (n >= 2) stats::t.test(.data[[value]], mu = 0)$p.value else NA_real_,
                     wilcox_p = if (n >= 2) suppressWarnings(stats::wilcox.test(.data[[value]], mu = 0)$p.value) else NA_real_,
                     .groups = "drop") |>
    dplyr::mutate(t_p_BH = p.adjust(t_p, "BH"), wilcox_p_BH = p.adjust(wilcox_p, "BH"))
}
readr::write_tsv(ses_rich, file.path(OUT, "Figure5A_SES_MPD_richness_per_sample.tsv"))
readr::write_tsv(ses_freq, file.path(OUT, "Figure5A_SES_MPD_frequency_per_sample.tsv"))
readr::write_tsv(vs_zero(ses_rich, "mpd.obs.z", season, depth_m), file.path(OUT, "Figure5A_SES_MPD_richness_vs0.tsv"))
readr::write_tsv(vs_zero(ses_freq, "mpd.obs.z", season, depth_m), file.path(OUT, "Figure5A_SES_MPD_frequency_vs0.tsv"))
readr::write_tsv(dplyr::select(ann_rich, season, p), file.path(OUT, "Figure5A_SES_MPD_richness_depth_Welch.tsv"))
readr::write_tsv(dplyr::select(ann_freq, season, p), file.path(OUT, "Figure5A_SES_MPD_frequency_depth_Welch.tsv"))
readr::write_tsv(dfC, file.path(OUT, "Figure5B_SES_FDis_per_sample.tsv"))
readr::write_tsv(vs_zero(dfC, "SES_FDis", season, depth_m), file.path(OUT, "Figure5B_SES_FDis_vs0.tsv"))
readr::write_tsv(dplyr::select(annC, season, p), file.path(OUT, "Figure5B_SES_FDis_depth_Welch.tsv"))
readr::write_tsv(dfD, file.path(OUT, "Figure5C_SES_CWSD_per_sample.tsv"))
readr::write_tsv(vs_zero(dfD, "SES_CWSD", trait, season, depth_m), file.path(OUT, "Figure5C_SES_CWSD_vs0.tsv"))
readr::write_tsv(dplyr::select(annD, dplyr::any_of(c("trait","season","p"))), file.path(OUT, "Figure5C_SES_CWSD_depth_Welch.tsv"))

# ===========================
# Second null for SES-FDis / SES-CWSD (pool-wide trait shuffle, null = "frequency").
# Run after the figure so the original richness-null results and RNG stream are untouched.
# ===========================
set.seed(123)
SES_FDis_freq <- ses_fdis(comm, traits_use, null = "frequency", R = 999, ncores = 1)
dfC_freq <- tibble(sample = names(SES_FDis_freq), SES_FDis = as.numeric(SES_FDis_freq)) %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(season), !is.na(depth_m), is.finite(SES_FDis))
set.seed(123)
SES_CWSD_freq <- ses_fdis_trait(comm, traits_use, null = "frequency", R = 999, ncores = 1)
dfD_freq <- as.data.frame(SES_CWSD_freq) %>%
  tibble::rownames_to_column("trait") %>%
  pivot_longer(-trait, names_to = "sample", values_to = "SES_CWSD") %>%
  left_join(meta, by = "sample") %>%
  filter(trait %in% panel_traits, !is.na(season), !is.na(depth_m), is.finite(SES_CWSD))
readr::write_tsv(dfC_freq, file.path(OUT, "Figure5B_SES_FDis_freqnull_per_sample.tsv"))
readr::write_tsv(vs_zero(dfC_freq, "SES_FDis", season, depth_m), file.path(OUT, "Figure5B_SES_FDis_freqnull_vs0.tsv"))
readr::write_tsv(dfD_freq, file.path(OUT, "Figure5C_SES_CWSD_freqnull_per_sample.tsv"))
readr::write_tsv(vs_zero(dfD_freq, "SES_CWSD", trait, season, depth_m), file.path(OUT, "Figure5C_SES_CWSD_freqnull_vs0.tsv"))
