# utils.R: shared palettes, helpers, and data loader for all Comau figure scripts
# Usage: source("utils.R"); ps <- load_ps()


# ---- Canonical season palette (Wong 2011 colorblind-friendly; Fig. 1 reference) ----
pal_season  <- c("Summer" = "#E69F00", "Fall" = "#D55E00", "Winter" = "#56B4E9")
shape_depth <- c("5" = 16, "20" = 17)

pal_phylum <- c(
  "Actinobacteriota"   = "#E41A1C",
  "Bacteroidota"       = "#377EB8",
  "Chloroflexota"      = "#4DAF4A",
  "Cyanobacteria"      = "#984EA3",
  "Desulfobacterota_B" = "#FF7F00",
  "Desulfobacterota_D" = "#A65628",
  "Gemmatimonadota"    = "#F781BF",
  "Margulisbacteria"   = "#999999",
  "Marinisomatota"     = "#66C2A5",
  "Myxococcota"        = "#FC8D62",
  "Nitrospinota"       = "#8DA0CB",
  "Patescibacteria"    = "#E78AC3",
  "Planctomycetota"    = "#A6D854",
  "Proteobacteria"     = "#FFD92F",
  "SAR324"             = "#E5C494",
  "Thermoplasmatota"   = "#B3B3B3",
  "Verrucomicrobiota"  = "#B15928",
  "Other"              = "#BFBFBF",
  # GTDB R232 names; renamed phyla keep their R202 colors
  "Actinomycetota"     = "#E41A1C",
  "Cyanobacteriota"    = "#984EA3",
  "Desulfobacterota"   = "#FF7F00",
  "Fidelibacterota"    = "#66C2A5",
  "Myxococcota_A"      = "#FC8D62",
  "Patescibacteriota"  = "#E78AC3",
  "Pseudomonadota"     = "#FFD92F"
)

# ---- Season assignment (Southern Hemisphere) ----
get_season_sh <- function(d) {
  m <- as.integer(format(as.Date(d), "%m"))
  ifelse(m %in% c(12, 1, 2), "Summer",
         ifelse(m %in% c(3, 4, 5), "Fall",
                ifelse(m %in% c(6, 7, 8), "Winter", NA_character_)))
}

# ---- Load and annotate phyloseq object (season + depth as ordered factors) ----
load_ps <- function(path = "data/Comau_ps90.RDS") {
  ps <- readRDS(path)
  sd <- as.data.frame(phyloseq::sample_data(ps))
  sd$Collection_Date <- as.Date(sd$Collection_Date)
  sd$season  <- factor(get_season_sh(sd$Collection_Date), levels = c("Summer", "Fall", "Winter"))
  sd$depth_m <- factor(as.character(sd$depth_m), levels = c("5", "20"))
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sd)
  ps
}

# ---- VIF filtering (iterative; drops predictors until all VIF ≤ thr) ----
drop_high_vif_vegan <- function(Xdf, Y, thr = 10) {
  keep <- colnames(Xdf)
  repeat {
    mod <- vegan::rda(as.matrix(Y) ~ ., data = as.data.frame(Xdf[, keep, drop = FALSE]))
    v   <- vegan::vif.cca(mod)
    if (max(v, na.rm = TRUE) > thr && length(keep) > 1)
      keep <- setdiff(keep, names(which.max(v)))
    else break
  }
  Xdf[, keep, drop = FALSE]
}

# ---- MRM: field names and environmental candidates (override locally if needed) ----
depth_var <- "depth_m"; line_var <- "line"; site_var <- "site"
lat_var   <- "latitude"; lon_var  <- "longitude"; date_var <- "Collection_Date"

env_candidates <- c(
  "temperature_C", "salinity_PSU", "pH",
  "nitrate_nitrite_uM", "nitrate_uM", "nitrite_uM",
  "phosphate_uM", "silicate_uM",
  "chlorophyll_A_ugL", "light_lux"
)

# ---- Circular mean of day-of-year (year-agnostic) ----
circular_mean_doy <- function(doy_vec, year_len = 365) {
  x <- doy_vec[!is.na(doy_vec)]
  if (!length(x)) return(NA_real_)
  ang <- 2 * pi * (x / year_len)
  mu  <- atan2(mean(sin(ang)), mean(cos(ang)))
  if (mu < 0) mu <- mu + 2 * pi
  round((mu / (2 * pi)) * year_len)
}

# ---- Build station-level community + metadata at a given depth ----
build_station_level <- function(ps, depth_value, depth_tolerance = 0.6) {
  md_all <- as(phyloseq::sample_data(ps), "data.frame")
  stopifnot(all(c(depth_var, line_var, site_var, lat_var, lon_var) %in% names(md_all)))

  md_all[[depth_var]] <- suppressWarnings(as.numeric(as.character(md_all[[depth_var]])))
  keep_idx <- !is.na(md_all[[depth_var]]) &
    !is.na(md_all[[lat_var]]) & !is.na(md_all[[lon_var]]) &
    abs(md_all[[depth_var]] - depth_value) <= depth_tolerance

  ps_d <- phyloseq::prune_samples(keep_idx, ps)
  ps_d <- phyloseq::prune_samples(phyloseq::sample_sums(ps_d) > 0, ps_d)
  if (phyloseq::nsamples(ps_d) < 2)
    stop(sprintf("Only %d sample(s) near %.1f m; need >=2.", phyloseq::nsamples(ps_d), depth_value))

  ps_d <- phyloseq::transform_sample_counts(ps_d, function(x) if (sum(x) > 0) x / sum(x) else x)
  md   <- as(phyloseq::sample_data(ps_d), "data.frame")
  md$.station_id <- paste0(md[[line_var]], "-", md[[site_var]])
  md$.DOY_sample <- NA_real_
  if (date_var %in% names(md) && any(!is.na(md[[date_var]])))
    md$.DOY_sample <- lubridate::yday(as.Date(as.character(md[[date_var]])))

  X <- as(phyloseq::otu_table(ps_d), "matrix")
  if (phyloseq::taxa_are_rows(ps_d)) X <- t(X)
  g     <- split(seq_len(nrow(X)), md$.station_id)
  X_sta <- t(sapply(g, function(ix) colMeans(X[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(X_sta) <- names(g)

  env_present <- intersect(env_candidates, names(md))
  md_sta <- md |>
    dplyr::select(.station_id, dplyr::all_of(
      c(lat_var, lon_var, line_var, site_var, ".DOY_sample", env_present))) |>
    dplyr::group_by(.station_id) |>
    dplyr::summarise(
      !!lat_var  := mean(.data[[lat_var]],  na.rm = TRUE),
      !!lon_var  := mean(.data[[lon_var]],  na.rm = TRUE),
      !!line_var := dplyr::first(.data[[line_var]]),
      !!site_var := dplyr::first(.data[[site_var]]),
      .DOY = circular_mean_doy(.DOY_sample),
      dplyr::across(dplyr::all_of(env_present), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    as.data.frame()
  rownames(md_sta) <- md_sta$.station_id

  common <- intersect(rownames(X_sta), rownames(md_sta))
  list(
    X_sta   = X_sta[common, , drop = FALSE],
    md_sta  = md_sta[common, , drop = FALSE],
    env_used = intersect(env_present, colnames(md_sta))
  )
}

# ---- Build distance matrices: community (A), geographic (B), env (C), temporal (D) ----
build_distances <- function(X_sta, md_sta, use_time = FALSE) {
  A      <- vegan::vegdist(X_sta, method = "bray")
  coords <- md_sta[, c(lon_var, lat_var)]
  B      <- as.dist(geosphere::distm(coords, fun = geosphere::distGeo) / 1000)

  env_cols <- setdiff(colnames(md_sta),
                      c(".station_id", line_var, site_var, lat_var, lon_var, ".DOY"))
  E <- md_sta[, env_cols, drop = FALSE]
  if (ncol(E) > 0) {
    ok <- vapply(E, function(v) sum(!is.na(v)) > 1 && sd(v, na.rm = TRUE) > 0, logical(1))
    E  <- E[, ok, drop = FALSE]
  }
  C <- if (ncol(E) > 0) dist(scale(E)) else NULL

  D <- NULL
  if (use_time && ".DOY" %in% names(md_sta) && any(!is.na(md_sta$.DOY))) {
    doy <- md_sta$.DOY
    if (length(unique(stats::na.omit(doy))) > 1) {
      Dmat <- outer(doy, doy, function(a, b) { d <- abs(a - b); pmin(d, 365 - d) })
      D    <- as.dist(Dmat)
    }
  }
  list(A = A, B = B, C = C, D = D)
}

# ---- Robust MRM: drops zero-variance and collinear predictors before fitting ----
run_mrm <- function(A, B, C = NULL, D = NULL, nperm = 9999) {
  vec <- function(d) if (is.null(d)) NULL else as.vector(d)
  df  <- data.frame(A = as.vector(A))
  if (!is.null(B)) df$B <- vec(B)
  if (!is.null(C)) df$C <- vec(C)
  if (!is.null(D)) df$D <- vec(D)
  preds <- intersect(c("B", "C", "D"), names(df))
  if (!length(preds)) stop("No predictors available (B/C/D all NULL).")
  df    <- stats::na.omit(df)
  if (!nrow(df)) stop("All pairwise rows are NA after combining predictors.")

  sds   <- sapply(df[, preds, drop = FALSE], sd, na.rm = TRUE)
  preds <- names(sds)[sds > .Machine$double.eps]
  if (!length(preds)) stop("All predictors have zero variance.")

  drop_coll <- character(0)
  if (length(preds) > 1) {
    repeat {
      Xs <- scale(as.matrix(df[, preds, drop = FALSE]))
      if (qr(Xs)$rank == ncol(Xs)) break
      R     <- abs(cor(Xs)); diag(R) <- 0
      worst <- names(which.max(colMeans(R, na.rm = TRUE)))
      drop_coll <- c(drop_coll, worst)
      preds <- setdiff(preds, worst)
      if (length(preds) <= 1) break
    }
  }

  X_final <- as.data.frame(scale(as.matrix(df[, preds, drop = FALSE])))
  df_fit  <- data.frame(A = df$A, X_final)
  form    <- as.formula(paste("A ~", paste(names(X_final), collapse = " + ")))
  fit     <- ecodist::MRM(form, data = df_fit, nperm = nperm)

  list(
    mrm                  = fit,
    used_predictors      = names(X_final),
    dropped_zero_var     = setdiff(names(sds), preds),
    dropped_collinearity = drop_coll
  )
}
