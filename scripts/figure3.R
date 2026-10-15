# figure3.R
# Figure 3 and Supplementary Figure S2: abundance-weighted functional potential by season and depth, sentinel genes, and functional PCoA/PERMANOVA.
# Run from the repository root (see README.md).

source("utils.R")
OUT <- "results/functional"

# ---- Packages ----
cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "purrr", "stringr", "readr",
  "ggplot2", "ggrepel", "scales", "viridisLite",
  "vegan", "ape",
  "cowplot", "patchwork"
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

# ---- Inputs ----
ps      <- load_ps("data/Comau_ps90_cov10_r232.RDS")
hits_fp <- "data/functional/gene_bin_mat.collapsed.csv"  # genes x MAGs (677 MAGs)

# ---- MAG abundance table A (samples x dereplicated MAGs) ----
otu <- as(otu_table(ps), "matrix")
A_counts <- if (taxa_are_rows(ps)) t(otu) else otu   # samples x MAGs
# row-wise relative abundances
row_sums <- rowSums(A_counts)
A_rel <- sweep(A_counts, 1, ifelse(row_sums==0, 1, row_sums), "/")
A_rel[is.na(A_rel)] <- 0
MAG_derep <- colnames(A_rel)                          # 513 MAGs (dereplicated)

# ---- Hits matrix H (genes x MAGs, 677 MAGs in columns) ----
H_tbl <- read_csv(hits_fp, show_col_types = FALSE)
stopifnot("gene" %in% names(H_tbl))
genes <- H_tbl$gene
H_mat <- as.matrix(H_tbl[ , setdiff(names(H_tbl), "gene") ])
rownames(H_mat) <- genes
MAG_hits <- colnames(H_mat)

# ---- MAG ID reconciliation: intersection first ----
MAG_common <- intersect(MAG_derep, MAG_hits)

message(sprintf("Dereplicated MAGs in ps: %d", length(MAG_derep)))
message(sprintf("MAGs with hits (v6):     %d", length(MAG_hits)))
message(sprintf("Intersection MAGs:        %d", length(MAG_common)))

# QC report v1: how many derep MAGs lack hit annotations?
mag_missing_hits <- setdiff(MAG_derep, MAG_common)
message(sprintf("Derep MAGs without hits:  %d", length(mag_missing_hits)))

# QC report v2: per-sample coverage of A contributed by MAGs that have hits
cov_per_sample <- rowSums(A_rel[ , MAG_common, drop = FALSE])
cov_summary <- c(min = min(cov_per_sample), q25 = quantile(cov_per_sample, 0.25),
                 median = median(cov_per_sample), q75 = quantile(cov_per_sample, 0.75),
                 max = max(cov_per_sample))
print(round(cov_summary, 3))

# (Optional) write coverage table for record
cov_tab <- data.frame(sample = rownames(A_rel),
                      season = sample_data(ps)$season,
                      coverage = cov_per_sample,
                      row.names = NULL)
# write_csv(cov_tab, "coverage_per_sample_ps_vs_hits.csv")


##### Sample × gene functional potentials (MAG-weighted) #####

# ---- Restrict to aligned MAGs ----
A_rel_use <- A_rel[, MAG_common, drop = FALSE]          # samples x M
H_use     <- H_mat[, MAG_common, drop = FALSE]          # genes   x M

# ---- Compute potentials: P = A %*% t(H) (samples x genes) ----
P_raw <- A_rel_use %*% t(H_use)                          # samples x genes
P_raw <- as.matrix(P_raw)
colnames(P_raw) <- rownames(H_use)                       # gene names
rownames(P_raw) <- rownames(A_rel_use)                   # sample IDs

# ---- Relative potentials per sample (compositional fingerprints) ----
row_tot <- rowSums(P_raw)
P_rel <- sweep(P_raw, 1, ifelse(row_tot == 0, 1, row_tot), "/")
P_rel[is.na(P_rel)] <- 0

# ---- Season factor for downstream tests/plots ----
season_vec <- sample_data(ps)$season
stopifnot(length(season_vec) == nrow(P_rel))

###### Transformations & pathway scaffolding #####


# =========================
# Step 3: Transformations & pathway scaffolding (tidy)
# =========================

# ---------- inputs from Step 2 ----------
# ps        : phyloseq object with sample_data(ps)$season already set (Summer/Fall/Winter)
# P_raw     : matrix | samples x genes (MAG-weighted hit potentials)
# P_rel     : matrix | samples x genes (row-normalized P_raw)
stopifnot(nrow(P_raw) == nrow(P_rel), all(colnames(P_raw) == colnames(P_rel)))

# ---------- attach sample & season; keep tibbles ----------
sd <- as.data.frame(sample_data(ps))
season_vec <- sd$season
stopifnot(length(season_vec) == nrow(P_raw))

P_raw_tbl <- as_tibble(P_raw, rownames = "sample") %>%
  mutate(season = season_vec, .after = "sample")

P_rel_tbl <- as_tibble(P_rel, rownames = "sample") %>%
  mutate(season = season_vec, .after = "sample")

# Hellinger transform on P_rel (for ordination/PERMANOVA)
# Note: vegan::decostand expects a numeric matrix (no sample/season columns)
gene_cols <- setdiff(names(P_rel_tbl), c("sample", "season"))
P_hell <- vegan::decostand(as.matrix(P_rel_tbl[, gene_cols]), method = "hellinger")
rownames(P_hell) <- P_rel_tbl$sample  # for later multivariate steps

# ---------- pathway dictionary ----------
pathway_dict <- list(
  Respiration_ETC = c("CydA","CyoA","CoxA","CcoN","NuoF","FrdA","SdhA","AtpA","FCC"),
  Carbon_core     = c("RbcL","PEPC","AcnB","CA","CoxL","MDH","Mcr"),   # Mcr = malonyl-CoA reductase (carbon fixation), not methyl-CoM reductase
  B12_cobalamin   = c("BluB","BtuCDF","CobA","CobC","CobU","CbiG","IsiB"),
  DMSP_Methyl     = c("dmdA","dmdB","dmdD","dddL","DddP","mddA","dmoA","Tmm","mtoX","dorA"),
  Nitrogen_cycle  = c("AmoA","AmoB","AmoC","NxrA","NapA","NarG","NirK","NirS","NrfA","NorB","NosZ"),
  Sulfur_cycle    = c("Sqr","SoxB","DsrA"),
  Methane_redox   = c("PmoA"),
  Phototrophy     = c("PsaA","PsbA","RHO"),
  Hydrogenases    = c("NiFe"),
  Detox_misc      = c("ArsC","RdhA","ARO","FdhA")
)

# helper to intersect with available genes safely
avail <- function(g) intersect(g, gene_cols)

# ---------- pathway totals per sample (from P_raw) ----------
pathway_totals <- tibble(sample = P_raw_tbl$sample, season = P_raw_tbl$season)
for (nm in names(pathway_dict)) {
  pathway_totals[[nm]] <- rowSums(P_raw_tbl[, avail(pathway_dict[[nm]]), drop = FALSE])
}

# ---------- pathway eigengenes (PC1 within each pathway from P_rel) ----------
get_pc1 <- function(mat) {
  mat <- mat[, apply(mat, 2, function(v) sd(v) > 0), drop = FALSE]   # cov10func: drop genes with no potential in any sample (e.g., Mcr absent from the 513 MAGs)
  if (ncol(mat) == 0) return(rep(NA_real_, nrow(mat)))      # no genes available in this pathway
  if (ncol(mat) == 1) return(as.numeric(scale(mat)[,1]))    # single gene → z-score
  pc <- prcomp(mat, center = TRUE, scale. = TRUE)
  pc$x[, 1]
}

pathway_eigengenes <- tibble(sample = P_rel_tbl$sample, season = P_rel_tbl$season)
for (nm in names(pathway_dict)) {
  mat <- as.matrix(P_rel_tbl[, avail(pathway_dict[[nm]]), drop = FALSE])
  pathway_eigengenes[[nm]] <- get_pc1(mat)
}

# ---------- QC: per-season means ± SE for pathway totals ----------
pathway_summary <- pathway_totals %>%
  pivot_longer(-c(sample, season), names_to = "pathway", values_to = "potential") %>%
  group_by(pathway, season) %>%
  summarise(
    mean = mean(potential, na.rm = TRUE),
    se   = sd(potential,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ---------- QC: seasonal effect sizes (genes) ----------
compute_effect <- function(df) {
  # df: columns 'val' and 'season'
  fit <- aov(val ~ season, data = df)
  sm  <- summary(fit)[[1]]
  ss  <- sm[,"Sum Sq"]; df_ <- sm[,"Df"]; ms <- sm[,"Mean Sq"]
  ss_season <- ss[1]; df_season <- df_[1]
  ss_resid  <- ss[2]; df_resid  <- df_[2]
  ms_resid  <- ms[2]
  omega2 <- max((ss_season - df_season * ms_resid) / (sum(ss) + ms_resid), 0)
  tibble(F = sm[1,"F value"], p = sm[1,"Pr(>F)"], omega2 = omega2)
}

gene_effects <- P_raw_tbl %>%
  pivot_longer(-c(sample, season), names_to = "gene", values_to = "val") %>%
  group_by(gene) %>%
  group_modify(~ compute_effect(.x)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(desc(omega2))

# ---------- QC: seasonal effect sizes (pathway totals) ----------
pathway_effects <- pathway_totals %>%
  pivot_longer(-c(sample, season), names_to = "pathway", values_to = "val") %>%
  group_by(pathway) %>%
  group_modify(~ compute_effect(.x)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(desc(omega2))

# ---------- ready-to-use tidy objects ----------
# P_raw_tbl  : sample, season, <gene columns>
# P_rel_tbl  : sample, season, <gene columns>
# P_hell     : numeric matrix (samples x genes) for ordination/PERMANOVA (rownames = sample)
# pathway_totals     : sample, season, <pathway columns>   (absolute potentials)
# pathway_eigengenes : sample, season, <pathway columns>   (PC1 scores)
# pathway_summary    : per-season means & SE (totals)
# gene_effects       : per-gene F, p, omega2, BH-adjusted p
# pathway_effects    : per-pathway F, p, omega2, BH-adjusted p


# pathway_dict already defined in Step 3
# gene_effects already computed (columns: gene, F, p, omega2, p_adj)
# P_raw_tbl already available (sample, season, <genes>)

# 1) Make a tidy pathway→gene lookup
pathway_lookup <- enframe(pathway_dict, name = "pathway", value = "gene") %>%
  unnest_longer(gene) %>%
  rename(gene = gene)

# (Optional) keep only genes that actually exist in the effects table
pathway_lookup <- pathway_lookup %>%
  semi_join(gene_effects, by = "gene")

# 2) Join effects and select top 2-3 sentinels per pathway by ω²
sentinels <- gene_effects %>%
  inner_join(pathway_lookup, by = "gene") %>%
  group_by(pathway) %>%
  arrange(desc(omega2), .by_group = TRUE) %>%
  slice_head(n = 3) %>%                     # change to 2 if you prefer two sentinels
  ungroup()








# =========================
# Step 4: Figure 3 Panels (seasonal; Phototrophy by eigengene)
# =========================

# ---------- Tidy lookup (pathway <-> gene) ----------
pathway_lookup <- enframe(pathway_dict, name = "pathway", value = "gene") %>%
  unnest_longer(gene) %>%
  rename(gene = gene)

# ---------- Display names for pathways (figure labels only) ----------
pretty_pw <- c(DMSP_Methyl = "DMSP", Carbon_core = "Carbon core", Sulfur_cycle = "Sulfur cycle", B12_cobalamin = "Cobalamin (B12)", Nitrogen_cycle = "Nitrogen cycle", Phototrophy = "Phototrophy", Methane_redox = "Methane cycling", Respiration_ETC = "Respiration", Hydrogenases = "Hydrogenases", Detox_misc = "Detoxification")
short_pw  <- c(DMSP_Methyl = "DMSP", Carbon_core = "Carbon", Sulfur_cycle = "Sulfur", B12_cobalamin = "B12", Nitrogen_cycle = "Nitrogen", Phototrophy = "Photo.")
# ---------- Which pathways go in the main figure ----------
main_pathways <- c("DMSP_Methyl", "Carbon_core", "Sulfur_cycle",   # DMSP replaces methane cycling (no seasonal effect; PmoA-like hits only)
                   "B12_cobalamin", "Nitrogen_cycle")

# -------------------------
# PANEL A with depth (linetype)
# -------------------------
meta <- as(sample_data(ps), "data.frame") %>%
  tibble::rownames_to_column("sample") %>%
  select(sample, season, depth_m)

# Totals for main pathways
ptot_long <- pathway_totals %>%
  select(sample, season, all_of(main_pathways)) %>%
  pivot_longer(-c(sample, season), names_to = "pathway", values_to = "potential") %>%
  left_join(meta, by = c("sample","season")) %>%
  mutate(source = "Total")

# Phototrophy eigengene
peig_long <- pathway_eigengenes %>%
  select(sample, season, Phototrophy) %>%
  pivot_longer(-c(sample, season), names_to = "pathway", values_to = "potential") %>%
  left_join(meta, by = c("sample","season")) %>%
  mutate(source = "Eigengene")

ptA <- bind_rows(ptot_long, peig_long) %>%
  filter((pathway %in% main_pathways & source == "Total") |
           (pathway == "Phototrophy"   & source == "Eigengene")) %>%
  mutate(
    pathway = factor(pathway, levels = c(main_pathways, "Phototrophy")),
    season  = factor(season,  levels = c("Summer","Fall","Winter")),
    depth_m = factor(as.character(depth_m), levels = c("5","20"))
  )

panelA_summary <- ptA %>%
  group_by(pathway, season, depth_m) %>%
  summarise(mean = mean(potential, na.rm = TRUE),
            se   = sd(potential,   na.rm = TRUE)/sqrt(n()),
            .groups = "drop")

panelA_plot <- panelA_summary %>%
  ggplot(aes(x = season, y = mean, color = pathway, linetype = depth_m, group = interaction(pathway, depth_m))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1, linewidth = 0.5) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 2, labeller = as_labeller(pretty_pw)) +
  scale_color_brewer(palette = "Dark2", labels = pretty_pw, name = "Pathway") +
  scale_linetype_manual(values = c("5" = "solid", "20" = "dashed"), name = "Depth (m)") +
  labs(x = NULL, y = "Functional potential (mean ± SE)",
       title = "Seasonal pathway fingerprints (Phototrophy by eigengene)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

# -------------------------
# PANEL B with depth (linetype)
# -------------------------
# Build sentinel list (top 2-3 per pathway)
sentinels <- gene_effects %>%
  inner_join(pathway_lookup, by = "gene") %>%
  group_by(pathway) %>%
  arrange(desc(omega2), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

# Ensure Phototrophy sentinels appear
phototrophy_sentinels <- intersect(c("PsbA","PsaA","RHO"), colnames(P_raw_tbl))
missing_ps <- setdiff(phototrophy_sentinels, sentinels$gene)
if (length(missing_ps) > 0) {
  sentinels <- bind_rows(
    sentinels,
    tibble(gene = missing_ps,
           F = NA_real_, p = NA_real_, omega2 = NA_real_, p_adj = NA_real_,
           pathway = "Phototrophy")
  )
}

# Sample metadata (for depth)
meta <- as(sample_data(ps), "data.frame") %>%
  tibble::rownames_to_column("sample") %>%
  select(sample, season, depth_m)

# Seasonal × depth summaries for sentinel genes
gene_long <- P_raw_tbl %>%
  pivot_longer(-c(sample, season), names_to = "gene", values_to = "potential") %>%
  left_join(meta, by = c("sample","season"))

sentinel_long <- gene_long %>%
  inner_join(sentinels %>% select(pathway, gene), by = "gene") %>%
  filter(!is.na(pathway)) %>%
  filter(pathway %in% c(main_pathways, "Phototrophy")) %>%
  mutate(depth_m = factor(as.character(depth_m), levels = c("5","20"))) %>%
  group_by(pathway, gene, season, depth_m) %>%
  summarise(mean = mean(potential, na.rm = TRUE),
            se   = sd(potential,   na.rm = TRUE)/sqrt(n()),
            .groups = "drop") %>%
  filter(!is.na(mean)) %>%
  mutate(
    season  = factor(season, levels = c("Summer","Fall","Winter")),
    pathway = factor(pathway, levels = c(main_pathways, "Phototrophy"))
  )

# Label positions: Winter endpoints (one label per gene per facet)
label_points <- sentinel_long %>%
  group_by(pathway, gene) %>%
  filter(as.integer(season) == max(as.integer(season))) %>%
  # Use the *average* of the two depths to place a single clean label
  summarise(y_lab = mean(mean, na.rm = TRUE), season = unique(season)[1], .groups = "drop")

panelB_plot <- ggplot(sentinel_long,
                      aes(x = season, y = mean,
                          group = interaction(gene, depth_m),
                          color = pathway, linetype = depth_m)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.1, linewidth = 0.5) +
  scale_color_brewer(palette = "Dark2", labels = pretty_pw, name = "Pathway") +
  scale_linetype_manual(values = c("5" = "solid", "20" = "dashed"), name = "Depth (m)") +
  facet_wrap(~ pathway, scales = "free_y", nrow = 2, drop = TRUE, labeller = as_labeller(pretty_pw)) +
  ggrepel::geom_text_repel(
    data = label_points,
    inherit.aes = FALSE,                                # <-- key line
    aes(x = season, y = y_lab, label = gene,            # only what we need
        # carry facet mapping so labels land in the right panels:
        # (color by pathway to match facet header color if you like)
    ),
    color = "black", size = 3,
    min.segment.length = 0,
    box.padding = 0.15, point.padding = 0.15, segment.size = 0.2, seed = 42
  ) +
  labs(x = NULL, y = "Sentinel gene potential (mean ± SE)",
       title = "Seasonal dynamics of sentinel genes per pathway (depth shown as linetype)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

# =========================
# PANEL C: Heatmap of sentinel seasonal z-scores (means)
# =========================

panelB_levels <- c(main_pathways, "Phototrophy")

# pool depths (equal-weight mean of the 5 m and 20 m means) so z is taken across
# three seasonal means per gene and each heatmap cell has one tile, as described in Methods.
sentinel_means <- sentinel_long %>%
  group_by(pathway, gene, season) %>%
  summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(season = factor(season, levels = c("Summer","Fall","Winter"))) %>%
  arrange(gene, season)

# z within gene across seasons
sentinel_z <- sentinel_means %>%
  group_by(gene) %>%
  mutate(z = if (sd(mean, na.rm = TRUE) > 0)
    (mean - mean(mean, na.rm = TRUE)) / sd(mean, na.rm = TRUE)
    else 0) %>%
  ungroup()

zlim <- max(abs(sentinel_z$z), na.rm = TRUE)

panelC_plot <- sentinel_z %>%
  mutate(pathway = factor(pathway, levels = panelB_levels)) %>%
  ggplot(aes(x = season, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(pathway ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(short_pw)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"),
                       midpoint = 0, limits = c(-zlim, zlim),
                       name = "z (seasonal mean)") +
  labs(x = NULL, y = NULL, title = "Sentinel gene seasonal z-scores (means)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(face = "bold"))


# =========================
# PANEL D: PCoA with season (color), depth (shape), and SHORT eigengene arrows
# (robust % axis labels)
# =========================

# 1) PCoA on Hellinger distances
dist_h  <- dist(P_hell, method = "euclidean")
pcoa_res <- ape::pcoa(dist_h)

# ---- Helper: robust % variance labels for axes across ape versions ----
get_pcoa_pct <- function(pcoa_res, k = 1, accuracy = 0.1) {
  vals <- pcoa_res$values
  # Some ape versions give a data.frame with named columns
  if (is.data.frame(vals)) {
    if ("Rel_corr_eig" %in% names(vals))   return(scales::percent(vals$Rel_corr_eig[k], accuracy = accuracy))
    if ("Relative_eig" %in% names(vals))   return(scales::percent(vals$Relative_eig[k], accuracy = accuracy))
    # Very old: second column often holds relative eig
    if (ncol(vals) >= 2)                   return(scales::percent(vals[k, 2],          accuracy = accuracy))
  }
  # Fallback: try vectors if present
  if (!is.null(vals$Rel_corr_eig))         return(scales::percent(vals$Rel_corr_eig[k], accuracy = accuracy))
  if (!is.null(vals$Relative_eig))         return(scales::percent(vals$Relative_eig[k], accuracy = accuracy))
  # Last resort
  return("NA")
}

pc1_pct <- get_pcoa_pct(pcoa_res, 1)
pc2_pct <- get_pcoa_pct(pcoa_res, 2)

# 2) Scores + metadata
scores <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(scores)[1:2] <- c("Axis.1","Axis.2")
scores$sample <- rownames(scores)

meta <- as(sample_data(ps), "data.frame") |>
  tibble::rownames_to_column("sample") |>
  dplyr::select(sample, season, depth_m)

scores <- scores |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::mutate(
    season  = factor(season, levels = c("Summer","Fall","Winter")),
    depth_m = factor(as.character(depth_m), levels = c("5","20"))
  )

# 3) Eigengene arrows (correlations → scaled to plot range)
eig_cols <- setdiff(colnames(pathway_eigengenes), c("sample","season"))
eig_mat  <- pathway_eigengenes |>
  dplyr::select(all_of(eig_cols)) |>
  as.matrix() |>
  scale()

ax1 <- pcoa_res$vectors[,1]
ax2 <- pcoa_res$vectors[,2]

arrow_raw <- tibble::tibble(
  pathway = eig_cols,
  x = as.numeric(cor(eig_mat, ax1, use = "pairwise.complete.obs")),
  y = as.numeric(cor(eig_mat, ax2, use = "pairwise.complete.obs"))
)

xr <- diff(range(scores$Axis.1)); yr <- diff(range(scores$Axis.2))
arrow_scale <- 0.35 * min(xr, yr)

arrow_df <- arrow_raw |>
  dplyr::mutate(xend = x * arrow_scale,
                yend = y * arrow_scale)

# 4) Plot (set axis labels ONCE using the robust pc1_pct/pc2_pct)
panelD_plot <- ggplot(scores, aes(x = Axis.1, y = Axis.2, color = season, shape = depth_m)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               inherit.aes = FALSE,
               color = "gray35",
               linewidth = 0.6,
               arrow = arrow(length = unit(0.16, "cm"))) +
  ggrepel::geom_text_repel(data = arrow_df,
                           aes(x = xend, y = yend, label = pretty_pw[as.character(pathway)]),
                           inherit.aes = FALSE,
                           color = "gray20",
                           size = 3,
                           segment.size = 0.2,
                           box.padding = 0.15,
                           point.padding = 0.15,
                           min.segment.length = 0,
                           seed = 42) +
  scale_color_manual(values = pal_season, name = "Season") +
  scale_shape_manual(values = shape_depth, name = "Depth (m)") +
  labs(
    title = "PCoA (Hellinger on MAG-weighted gene potentials) with pathway eigengenes",
    x = paste0("PCoA1 (", pc1_pct, ")"),
    y = paste0("PCoA2 (", pc2_pct, ")")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# =========================
# PERMANOVA (reportable)
# =========================
set.seed(1)
adon <- adonis2(P_hell ~ sample_data(ps)$season, permutations = 999, method = "euclidean")
print(adon)

adon_depth_season_full <- adonis2(
  P_hell ~ season * depth_m,
  data = meta,
  permutations = 999,
  method = "euclidean",
  by = "margin"   # gives marginal (per-term) tests
)
print(adon_depth_season_full)

# with an interaction in the model, by = "margin" tests only season:depth_m.
# Main effects come from the additive model, each tested after the other (marginal).
adon_main_effects <- adonis2(
  P_hell ~ season + depth_m,
  data = meta,
  permutations = 999,
  method = "euclidean",
  by = "margin"
)
print(adon_main_effects)

# Individual panel saves removed: see Figure3_final.pdf / .tiff for the combined figure.











# =========================
# Figure 3: assembly with a single master legend
# =========================

# --- Harmonize panel theming ---
theme_fig <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "white", colour = "grey70"),
    plot.margin      = margin(4, 6, 4, 6)
  )

A <- panelA_plot + theme_fig + labs(title = NULL) + theme(legend.position = "bottom")
B <- panelB_plot + theme_fig + labs(title = NULL) + theme(legend.position = "bottom")
C <- panelC_plot + theme_fig + labs(title = NULL) + theme(legend.position = "right")
D <- panelD_plot + theme_fig + labs(title = NULL) + theme(legend.position = "bottom")

# --- Extract component legends ---

# (i) Pathway (color) + Depth (linetype) for A/B
leg_path_depth <- cowplot::get_legend(
  A + guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    theme(
      legend.box.margin = margin(2, 2, 2, 2),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
)

# (ii) Panel C heatmap legend: HORIZONTAL colorbar
leg_heatmap <- cowplot::get_legend(
  panelC_plot +
    guides(fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(35, "mm"),
      barheight = unit(2.5, "mm")
    )) +
    theme(
      legend.position   = "bottom",
      legend.box.margin = margin(2, 2, 2, 2),
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8)
    )
)

# (iii) Season (color) + Depth (shape) for D
leg_season_depth <- cowplot::get_legend(
  D + guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme(
      legend.box.margin = margin(2, 2, 2, 2),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
)

# Wrap each legend so we can lay them out as one row
legA_plot <- cowplot::ggdraw(leg_path_depth)
legC_plot <- cowplot::ggdraw(leg_heatmap)
legD_plot <- cowplot::ggdraw(leg_season_depth)
spacer    <- cowplot::ggdraw() + theme(plot.margin = margin(0,0,0,0))

# --- Build ONE horizontal master legend row ---
master_legend <- cowplot::plot_grid(   # two rows: pathway/depth legend, then heat map and season/depth legends
  legA_plot,
  cowplot::plot_grid(legC_plot, legD_plot, nrow = 1, rel_widths = c(1, 1.4)),
  ncol = 1
)

# --- Strip legends from panels for the grid ---
A0 <- A + theme(legend.position = "none")
B0 <- B + theme(legend.position = "none")
C0 <- C + theme(legend.position = "none")
D0 <- D + theme(legend.position = "none")

# --- Main figure: pathways (A) on top, heat map (B) and PCoA (C) below; the sentinel-gene panel
#     moves to Supplementary Fig. S2 ---
panel_grid <- A0 / (C0 | D0) +
  patchwork::plot_layout(heights = c(0.85, 1))
suppS2 <- B + theme(legend.position = "bottom", legend.box = "vertical") + labs(title = NULL)
ggsave(file.path(OUT, "SuppFigS2_sentinel_genes.png"), suppS2, width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(file.path(OUT, "SuppFigS2_sentinel_genes.pdf"), suppS2, width = 10, height = 7, units = "in", device = cairo_pdf)
panel_tagged <- panel_grid + patchwork::plot_annotation(tag_levels = "A")

# --- Stack tagged panels + single master legend (no tags on legend row) ---
final_fig3 <- cowplot::plot_grid(
  panel_tagged,
  master_legend,
  ncol = 1,
  rel_heights = c(1, 0.17)
)

# --- Save ---
ggsave(file.path(OUT, "Figure3.pdf"),  final_fig3, width = 11, height = 11, units = "in", device = cairo_pdf)
ggsave(file.path(OUT, "Figure3.png"), final_fig3, width = 11, height = 11, units = "in", dpi = 300, bg = "white")

# --- Persist objects needed by Figure 4 / Figure 5 ---
saveRDS(list(A_rel_use = A_rel_use, H_use = H_use,
             pathway_dict = pathway_dict, main_pathways = main_pathways),
        file.path(OUT, "fig3_shared_objects.RDS"))



# ===========================
# Extra exports
# ===========================
# Alignment checks: season labels and PERMANOVA rows must match P_hell rows
readr::write_tsv(tibble::tibble(
  n_MAG_ps = phyloseq::ntaxa(ps), n_MAG_hits = ncol(H_mat), n_MAG_common = length(MAG_common),
  P_raw_rows_match_sample_data = identical(rownames(P_raw), rownames(sample_data(ps))),
  P_hell_rows_match_meta = identical(rownames(P_hell), meta$sample)),
  file.path(OUT, "Figure3_alignment_checks.tsv"))
readr::write_tsv(gene_effects,    file.path(OUT, "Figure3_gene_season_effects.tsv"))
readr::write_tsv(pathway_effects, file.path(OUT, "Figure3_pathway_totals_season_effects.tsv"))
# Same one-way ANOVA/omega2 applied to the eigengenes plotted for Phototrophy in panel A
eig_effects <- pathway_eigengenes %>%
  pivot_longer(-c(sample, season), names_to = "pathway", values_to = "val") %>%
  filter(!is.na(val)) %>% group_by(pathway) %>% group_modify(~ compute_effect(.x)) %>%
  ungroup() %>% mutate(p_adj = p.adjust(p, method = "BH"))
readr::write_tsv(eig_effects, file.path(OUT, "Figure3_pathway_eigengene_season_effects.tsv"))
readr::write_tsv(sentinels, file.path(OUT, "Figure3_sentinels.tsv"))
readr::write_tsv(panelA_summary_export <- ptA %>% group_by(pathway, season, depth_m) %>%
                   summarise(mean = mean(potential, na.rm = TRUE), se = sd(potential, na.rm = TRUE)/sqrt(n()),
                             n = n(), .groups = "drop"),
                 file.path(OUT, "Figure3A_pathway_means_by_season_depth.tsv"))
readr::write_tsv(sentinel_long, file.path(OUT, "Figure3B_sentinel_means_by_season_depth.tsv"))
readr::write_tsv(sentinel_z,    file.path(OUT, "Figure3C_heatmap_z.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(adon), "term"),
                 file.path(OUT, "Figure3D_PERMANOVA_season.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(adon_depth_season_full), "term"),
                 file.path(OUT, "Figure3D_PERMANOVA_season_x_depth_margin.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(adon_main_effects), "term"),
                 file.path(OUT, "Figure3D_PERMANOVA_main_effects_margin.tsv"))
# PC1 loadings behind each pathway eigengene (same prcomp call as get_pc1), so the sign of the
# eigengene can be interpreted; eigengene-gene correlations are added as a direct check.
eig_loadings <- dplyr::bind_rows(lapply(names(pathway_dict), function(nm) {
  mat <- as.matrix(P_rel_tbl[, avail(pathway_dict[[nm]]), drop = FALSE])
  mat <- mat[, apply(mat, 2, function(v) sd(v) > 0), drop = FALSE]   # cov10func: same zero-gene drop as get_pc1
  if (ncol(mat) < 2) return(NULL)
  pc <- prcomp(mat, center = TRUE, scale. = TRUE)
  tibble::tibble(pathway = nm, gene = rownames(pc$rotation), PC1_loading = pc$rotation[, 1],
                 PC1_var_explained = summary(pc)$importance[2, 1],
                 cor_with_eigengene = as.numeric(cor(mat, pc$x[, 1])))
}))
readr::write_tsv(eig_loadings, file.path(OUT, "Figure3_pathway_eigengene_loadings.tsv"))
