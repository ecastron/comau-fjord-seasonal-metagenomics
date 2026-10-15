# figure1.R
# Figure 1: sampling map, PCoA and PERMANOVA (season, depth, interaction), phylum composition, and abundance-weighted MPD.
# Run from the repository root (see README.md).

source("utils.R")
OUT <- "results/community"

# ---- Packages ----
cran_pkgs <- c(
  "vegan","ggplot2","dplyr","tibble","tidyr","forcats","scales",
  "patchwork","rstatix","multcompView","RColorBrewer","gridExtra","ggpubr",
  "sf","readr","ggrepel","ggspatial","rnaturalearth","shadowtext","ape","picante","ARTool"
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

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")

# ===========================
# PANEL A: Map with locator inset (bottom-left). No titles anywhere.
# Needs: metadata.tsv and CHL_adm_shp/CHL_adm2.shp
# ===========================
# Helpers
make_star_lonlat <- function(lon, lat, r_deg = 0.010, inner_ratio = 0.5, n = 5) {
  lat_rad <- lat * pi/180
  ang <- seq(0, 2*pi, length.out = 2*n + 1)[-(2*n + 1)]
  rad <- rep(c(r_deg, r_deg*inner_ratio), n)
  dx <- (rad * cos(ang)) / max(0.1, cos(lat_rad))
  dy <-  rad * sin(ang)
  data.frame(long = lon + dx, lat = lat + dy, group = "huinay_star")
}
label_lon <- function(x) paste0(formatC(abs(x), format = "f", digits = 2), "°", ifelse(x < 0, "W", "E"))
label_lat <- function(x) paste0(formatC(abs(x), format = "f", digits = 2), "°", ifelse(x < 0, "S", "N"))
sf_polygons_to_df <- function(x) {
  x <- st_cast(x, "POLYGON", warn = FALSE)
  out <- lapply(seq_len(nrow(x)), function(i) {
    cmat <- st_coordinates(x[i, ])
    data.frame(long = cmat[, "X"], lat = cmat[, "Y"], group = paste0("id", i, "_ring", cmat[, "L1"]))
  })
  do.call(rbind, out)
}

# Sites
meta <- readr::read_tsv("data/metadata.tsv", show_col_types = FALSE)
stopifnot(all(c("line","site","longitude","latitude") %in% names(meta)))
pts12 <- meta |>
  transmute(
    line_chr = gsub("^L", "", as.character(line)),
    site_chr = gsub("^S", "", as.character(site)),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  ) |>
  filter(is.finite(longitude), is.finite(latitude)) |>
  group_by(line_chr, site_chr) |>
  summarise(longitude = median(longitude), latitude = median(latitude), .groups = "drop") |>
  mutate(label = paste0("L", line_chr, "S", site_chr))

# Window
pad_x <- 0.05; pad_y <- 0.03
min.long <- min(pts12$longitude) - pad_x
max.long <- max(pts12$longitude) + pad_x
min.lat  <- min(pts12$latitude)  - pad_y
max.lat  <- max(pts12$latitude)  + pad_y

# Polygons and neatlines
adm2       <- sf::st_read("data/CHL_adm_shp/CHL_adm2.shp", quiet = TRUE) |> sf::st_transform(4326)
palena_wgs <- dplyr::filter(adm2, NAME_2 == "Palena")
win_poly_wgs <- sf::st_as_sfc(sf::st_bbox(c(xmin = min.long, xmax = max.long, ymin = min.lat, ymax = max.lat), crs = sf::st_crs(4326)))
palena_m     <- sf::st_transform(palena_wgs, 32718)
win_poly_m   <- sf::st_transform(win_poly_wgs, 32718)
palena_crop  <- suppressWarnings(sf::st_intersection(palena_m, sf::st_buffer(win_poly_m, 200))) |> sf::st_transform(4326)
palena_df    <- sf_polygons_to_df(palena_crop)

n_segments_x <- 4; n_segments_y <- 5
step_x <- (max.long - min.long) / n_segments_x
step_y <- (max.lat  - min.lat ) / n_segments_y
axis.breaks.x <- seq(min.long, max.long, by = step_x)
axis.breaks.y <- seq(min.lat,  max.lat,  by = step_y)
if (tail(axis.breaks.x,1) < max.long - 1e-10) axis.breaks.x <- c(axis.breaks.x, max.long)
if (tail(axis.breaks.y,1) < max.lat  - 1e-10) axis.breaks.y <- c(axis.breaks.y, max.lat)
grid.width <- (max.lat - min.lat) / 90
rects.long <- data.frame(x_start = rep(head(axis.breaks.x, -1), 2))
rects.long$x_end   <- rep(tail(axis.breaks.x, -1), 2)
rects.long$y_start <- c(rep(min.lat, length(axis.breaks.x)-1),
                        rep(max.lat - grid.width, length(axis.breaks.x)-1))
rects.long$y_end   <- c(rep(min.lat + grid.width, length(axis.breaks.x)-1),
                        rep(max.lat, length(axis.breaks.x)-1))
rects.long$color <- rep(c("black","white"), length.out = nrow(rects.long))
rects.lat <- data.frame(y_start = rep(head(axis.breaks.y, -1), 2))
rects.lat$y_end   <- rep(tail(axis.breaks.y, -1), 2)
rects.lat$x_start <- c(rep(min.long, length(axis.breaks.y)-1),
                       rep(max.long - grid.width, length(axis.breaks.y)-1))
rects.lat$x_end   <- c(rep(min.long + grid.width, length(axis.breaks.y)-1),
                       rep(max.long, length(axis.breaks.y)-1))
rects.lat$color <- rep(c("black","white"), length.out = nrow(rects.lat))
rects.black <- rbind(rects.long, rects.lat) |> dplyr::filter(color == "black")
rects.white <- rbind(rects.long, rects.lat) |> dplyr::filter(color == "white")

huinay <- data.frame(longitude = -72.415246, latitude = -42.379839)
huinay_star <- make_star_lonlat(huinay$longitude, huinay$latitude, r_deg = 0.008, inner_ratio = 0.5)
win_poly_sf <- sf::st_as_sf(win_poly_wgs)
chl <- rnaturalearth::ne_countries(country = "Chile", returnclass = "sf") |> sf::st_transform(4326)
bb <- sf::st_bbox(palena_wgs)
palena_box <- sf::st_sfc(sf::st_polygon(list(matrix(
  c(bb["xmin"], bb["ymin"], bb["xmin"], bb["ymax"], bb["xmax"], bb["ymax"],
    bb["xmax"], bb["ymin"], bb["xmin"], bb["ymin"]), ncol = 2, byrow = TRUE))), crs = sf::st_crs(palena_wgs))

main_map <- ggplot() +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = alpha("lightsteelblue1", 0.25), colour = NA),
    panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none",
    plot.margin = unit(c(2,2,2,2), "pt")
  ) +
  geom_sf(data = win_poly_sf, fill = NA, color = NA, inherit.aes = FALSE) +
  geom_polygon(data = palena_df, aes(x = long, y = lat, group = group),
               color = "grey35", fill = "grey95", linewidth = 0.6) +
  geom_polygon(data = huinay_star, aes(x = long, y = lat, group = group),
               fill = "black", color = "black", linewidth = 0.2) +
  shadowtext::geom_shadowtext(data = huinay,
                              aes(x = longitude, y = latitude, label = "Huinay\nStation"),
                              nudge_y = -0.00, nudge_x = 0.030, size = 3.2, colour = "black",
                              bg.colour = "white", bg.r = 0.08, fontface = "bold") +
  geom_point(data = pts12, aes(x = longitude, y = latitude),
             shape = 21, fill = "red", color = "black", size = 3.4, stroke = 0.6) +
  ggrepel::geom_label_repel(data = pts12, aes(x = longitude, y = latitude, label = label),
                            size = 3.2, fill = "white", label.size = 0, box.padding = 0.45, point.padding = 0.45,
                            max.overlaps = Inf, segment.color = NA, seed = 42) +
  coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
  scale_x_continuous(breaks = axis.breaks.x, sec.axis = dup_axis(), labels = label_lon) +
  scale_y_continuous(breaks = axis.breaks.y, sec.axis = dup_axis(), labels = label_lat) +
  geom_rect(data = rects.white, inherit.aes = FALSE,
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
            color = "black", fill = "white") +
  geom_rect(data = rects.black, inherit.aes = FALSE,
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
            color = "black", fill = "black") +
  # North arrow → TOP-LEFT, with a bit of inset from the edges
  annotation_north_arrow(
    location = "tl", which_north = "grid",
    style = north_arrow_minimal(text_size = 6, line_width = 0.5),
    pad_x = unit(0.7, "cm"),   # rightward from the left edge
    pad_y = unit(0.7, "cm")    # downward from the top edge
  ) +
  
  # Scale bar → still BOTTOM-LEFT, but raised “up a bit”
  annotation_scale(
    location = "bl",
    unit_category = "metric",
    style = "bar", width_hint = 0.25, text_cex = 0.8, line_width = 0.3,
    pad_x = unit(0.70, "cm"),
    pad_y = unit(1.00, "cm")   # increase this to raise it further
  )


locator <- ggplot() +
  geom_sf(data = chl, fill = "grey95", color = "grey35", linewidth = 0.3) +
  geom_sf(data = palena_box, fill = NA, color = "red", linewidth = 0.6) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "pt"))

# Build A (map) with locator inset (locator is ONLY used as inset, never as a top-level panel)
p_A_map <- main_map +
  patchwork::inset_element(
    locator + theme(plot.margin = margin(0,0,0,0)),
    left = 0.04, bottom = 0.06, right = 0.42, top = 0.46,
    align_to = "panel"
  )

# ===========================
# PANEL B: PCoA + PERMANOVA inset (no titles)
# ===========================
ps_rel   <- transform_sample_counts(ps, function(x) if (sum(x) > 0) x/sum(x) else x)
dist_mat <- phyloseq::distance(ps_rel, method = "bray")   # computed once; reused by PCoA + PERMANOVA
set.seed(1)
ord  <- ordinate(ps_rel, method = "PCoA", distance = dist_mat)
pdat <- plot_ordination(ps_rel, ord, justDF = TRUE)

p_B <- ggplot(pdat, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = season, shape = depth_m), size = 2.2, alpha = 0.9, na.rm = TRUE) +
  stat_ellipse(data = pdat,
               mapping = aes(x = Axis.1, y = Axis.2, color = season, group = season),
               inherit.aes = FALSE, type = "t", level = 0.95, linewidth = 0.5, alpha = 0.20, na.rm = TRUE) +
  labs(color = "Season", shape = "Depth (m)") +
  scale_color_manual(values = pal_season, drop = FALSE) +
  scale_shape_manual(values = shape_depth, breaks = c("5","20"), drop = FALSE) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_bw(base_size = 11)

# PERMANOVA (stratified by depth): reuses dist_mat from above
meta <- as(sample_data(ps_rel), "data.frame")
labs_dm <- attr(dist_mat, "Labels"); meta <- meta[labs_dm, , drop = FALSE]
meta$season  <- droplevels(meta$season); meta$depth_m <- droplevels(meta$depth_m)
keep <- !is.na(meta$season) & !is.na(meta$depth_m)
meta <- meta[keep, , drop = FALSE]; dmat <- as.matrix(dist_mat)[keep, keep, drop = FALSE]
adon <- vegan::adonis2(dmat ~ season, data = meta, permutations = 999, strata = meta$depth_m)
adon_tab <- as.data.frame(adon)
R2   <- adon_tab$R2[1]; pval <- adon_tab$`Pr(>F)`[1]
lab_perm <- sprintf("PERMANOVA (stratified by depth)\nR² = %.3f\nP = %.4f", R2, pval)

# pairwise by season, stratified by depth
pairwise_permanova_strat <- function(dmat, meta, group_col = "season", strata_col = "depth_m",
                                     permutations = 999, p_adjust = "BH") {
  g  <- droplevels(meta[[group_col]]); st <- droplevels(meta[[strata_col]])
  lv <- levels(g); out <- list(); k <- 1
  for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    pair <- c(lv[i], lv[j]); idx <- g %in% pair
    meta_sub <- droplevels(meta[idx, , drop = FALSE]); st_sub <- droplevels(st[idx])
    dsub_mat <- dmat[idx, idx, drop = FALSE]
    fit <- vegan::adonis2(dsub_mat ~ meta_sub[[group_col]], permutations = permutations, strata = st_sub)
    tab <- as.data.frame(fit)
    out[[k]] <- data.frame(Comparison = paste(pair, collapse = " vs "),
                           R2 = tab$R2[1], p = tab$`Pr(>F)`[1], stringsAsFactors = FALSE)
    k <- k + 1
  }
  df <- do.call(rbind, out); df$`p (BH)` <- p.adjust(df$p, method = p_adjust); df
}
pw_df <- pairwise_permanova_strat(dmat, meta, "season", "depth_m", 999, "BH")
ord_lev <- c("Summer vs Fall","Summer vs Winter","Fall vs Winter")
if (!is.null(pw_df) && nrow(pw_df)>0) {
  pw_df$Comparison <- factor(pw_df$Comparison, levels = ord_lev)
  pw_df <- pw_df[order(pw_df$Comparison), ]
}
pw_compact <- if (nrow(pw_df) == 0) tibble::tibble(Comparison="-", R2="-", `p (BH)`="-") else
  pw_df %>% transmute(Comparison, R2 = sprintf("%.3f", R2),
                      `p (BH)` = scales::pvalue(`p (BH)`, accuracy = 0.0001))

tab_theme <- gridExtra::ttheme_minimal(
  core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 7)),
  colhead = list(fg_params = list(fontface = "bold", fontsize = 7))
)
tab_pw_grob  <- gridExtra::tableGrob(pw_compact, rows = NULL, theme = tab_theme)
tab_pw_plot  <- ggpubr::as_ggplot(tab_pw_grob)

# add global PERMANOVA label (not a title)
p_B <- p_B +
  annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = lab_perm, size = 3, linewidth = 0.25)

# Build B with PERMANOVA table inset (table is ONLY used as inset)
p_B_inset <- p_B +
  patchwork::inset_element(
    tab_pw_plot + theme(plot.margin = margin(0,0,0,0)),
    left = 0.02, bottom = 0.58, right = 0.48, top = 0.98,
    align_to = "panel"
  )

# ===========================
# PANEL C: Depth-resolved phyla across calendar date (dominant ≥1% globally)
# ===========================
mat_counts <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) mat_counts <- t(mat_counts)
phy_vec <- as.vector(tax_table(ps)[, "Phylum"]); phy_vec[is.na(phy_vec) | phy_vec==""] <- "Unclassified"
names(phy_vec) <- rownames(tax_table(ps))
mat_phylum <- rowsum(mat_counts, group = phy_vec, reorder = TRUE)

smeta <- as(sample_data(ps), "data.frame") %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(Collection_Date = as.Date(Collection_Date),
         cal_md   = format(Collection_Date, "%m-%d"),
         cal_date = as.Date(paste0("2000-", cal_md)),
         depth_m  = factor(as.character(depth_m), levels = c("5","20")))
smeta <- smeta[smeta$Sample %in% colnames(mat_phylum), , drop = FALSE]

dfC <- as.data.frame(t(mat_phylum)) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Phylum", values_to = "Count") %>%
  left_join(smeta, by = "Sample") %>%
  filter(!is.na(cal_date), !is.na(depth_m))

dom_phyla <- dfC %>%
  group_by(Phylum) %>% summarise(Total = sum(Count), .groups = "drop") %>%
  mutate(Prop = Total / sum(Total)) %>% filter(Prop >= 0.01) %>% pull(Phylum) %>% as.character()

dfC <- dfC %>%
  mutate(PhylumCollapsed = if_else(Phylum %in% dom_phyla, Phylum, "Other")) %>%
  mutate(PhylumCollapsed = factor(PhylumCollapsed, levels = names(pal_phylum)))

df_day <- dfC %>%
  group_by(cal_date, depth_m, PhylumCollapsed) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

p_C <- ggplot(df_day, aes(x = cal_date, y = Count, fill = PhylumCollapsed)) +
  geom_area(position = "fill", linewidth = 0.1) +
  facet_wrap(~ depth_m, ncol = 1, labeller = labeller(depth_m = c(`5`="5 m", `20`="20 m"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_fill_manual(values = pal_phylum, drop = TRUE) +
  labs(x = "Calendar date", y = "Relative abundance", fill = "Phylum") +
  theme_bw(base_size = 11)

# ===========================
# PANEL D: Abundance-weighted MPD (season × depth) + posthoc letters
# ===========================
stopifnot(!is.null(phy_tree(ps)))
comm <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) comm <- t(comm)
keep_samp <- rowSums(comm) > 0; comm <- comm[keep_samp, , drop = FALSE]
tree0 <- phy_tree(ps); tax_in_both <- intersect(colnames(comm), tree0$tip.label)
if (length(tax_in_both) < 3) stop("Too few overlapping taxa between OTU table and tree.")
comm <- comm[, tax_in_both, drop = FALSE]; tree <- ape::keep.tip(tree0, tax_in_both); tree$tip.label <- tax_in_both
comm_rel <- sweep(comm, 1, rowSums(comm), "/")
D <- cophenetic(tree); D <- D[colnames(comm_rel), colnames(comm_rel)]
aMPD <- picante::mpd(samp = comm_rel, dis = D, abundance.weighted = TRUE)

sd2 <- as(sample_data(ps), "data.frame"); sd2$Sample <- rownames(sd2)
sd2$season  <- droplevels(factor(sd2$season, levels = c("Summer","Fall","Winter")))
sd2$depth_m <- droplevels(factor(as.character(sd2$depth_m), levels = c("5","20")))
adiv_pd <- data.frame(Sample = rownames(comm_rel), aMPD = as.numeric(aMPD)) |>
  dplyr::left_join(sd2[, c("Sample","season","depth_m")], by = "Sample") |>
  dplyr::filter(!is.na(season), !is.na(depth_m))

art_fit <- ARTool::art(aMPD ~ season * depth_m, data = adiv_pd)
art_tab <- anova(art_fit)

dunn_depth <- adiv_pd |>
  dplyr::group_by(depth_m) |>
  rstatix::dunn_test(aMPD ~ season, p.adjust.method = "BH") |>
  dplyr::ungroup()

letters_by_depth <- dunn_depth |>
  dplyr::mutate(comp = paste(group1, group2, sep = "-")) |>
  dplyr::select(depth_m, comp, p.adj) |>
  dplyr::group_split(depth_m) |>
  lapply(function(df) {
    vec <- setNames(df$p.adj, df$comp)
    L   <- multcompView::multcompLetters(vec, threshold = 0.05)$Letters
    tibble::tibble(season = names(L), .label = unname(L))
  })
names(letters_by_depth) <- unique(dunn_depth$depth_m)
letters_df <- dplyr::bind_rows(
  lapply(names(letters_by_depth), function(d) dplyr::mutate(letters_by_depth[[d]], depth_m = d))
)

wilc_season <- adiv_pd |>
  dplyr::group_by(season) |>
  rstatix::wilcox_test(aMPD ~ depth_m, exact = FALSE) |>
  rstatix::adjust_pvalue(method = "BH") |>
  dplyr::ungroup()

adiv_pd$depth_m  <- factor(as.character(adiv_pd$depth_m),  levels = c("5","20"), ordered = TRUE)
letters_df$depth_m <- factor(as.character(letters_df$depth_m), levels = c("5","20"), ordered = TRUE)

y_top_pd <- max(adiv_pd$aMPD, na.rm = TRUE)
p_D <- ggplot(adiv_pd, aes(season, aMPD, fill = season)) +
  geom_boxplot(width = 0.7, alpha = 0.85, outlier.shape = 21) +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.6) +
  geom_text(data = letters_df,
            aes(x = season, y = y_top_pd * 1.05, label = .label),
            inherit.aes = FALSE, vjust = 0, size = 3.4) +
  facet_grid(cols = vars(depth_m), drop = FALSE,
             labeller = labeller(depth_m = c(`5` = "5 m", `20` = "20 m"))) +
  scale_fill_manual(values = pal_season, drop = FALSE) +
  labs(x = NULL, y = "Abundance-weighted MPD") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ==============================================
# FINAL COMPOSITION (one block only)
# - B rendered onto ggplot canvas for alignment
# - Legends for B and C extracted to a right column
# ==============================================
if (!"ggplotify" %in% rownames(installed.packages()))
  install.packages("ggplotify", repos = "https://cloud.r-project.org")
library(ggplotify)

## Extract legends for B and C
leg_B <- cowplot::get_legend(
  p_B + theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.box.margin = margin(2,2,2,2)
  )
)
leg_C <- cowplot::get_legend(
  p_C + theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.box.margin = margin(2,2,2,2)
  )
)

## Remove legends from working panels
p_B_clean <- p_B + theme(legend.position = "none")
p_C_clean <- p_C + theme(legend.position = "none")
p_D_clean <- p_D + theme(legend.position = "none")

## Rebuild B with inset (no legend), then render to ggplot canvas
p_B_inset_nolegend <- p_B_clean +
  patchwork::inset_element(
    tab_pw_plot + theme(plot.margin = margin(0,0,0,0)),
    left = 0.02, bottom = 0.58, right = 0.48, top = 0.98, align_to = "panel"
  )
p_B_canvas <- ggplotify::as.ggplot(function() print(p_B_inset_nolegend))

## Align B with C on the left, then align (C, D) on the same left
aligned_BC <- cowplot::align_plots(p_B_canvas, p_C_clean, align = "h", axis = "l")
p_B_al <- aligned_BC[[1]]; p_C_al <- aligned_BC[[2]]
aligned_CD <- cowplot::align_plots(p_C_al, p_D_clean, align = "h", axis = "l")
p_C_al <- aligned_CD[[1]]; p_D_al <- aligned_CD[[2]]

## Bottom row (C | D) with labels
bottom_right_labeled <- cowplot::plot_grid(
  p_C_al, p_D_al, ncol = 2,
  labels = c("C","D"),
  label_fontface = "bold", label_size = 14,
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1
)

## Right column (B over C|D), then attach legend column at far right
right_col_no_leg <- cowplot::plot_grid(
  p_B_al, bottom_right_labeled, ncol = 1,
  labels = c("B", NULL),
  label_fontface = "bold", label_size = 14,
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1,
  rel_heights = c(1.05, 1.00),
  align = "v", axis = "lr"
)
legend_column <- cowplot::plot_grid(leg_B, leg_C, ncol = 1, align = "v")
right_with_legend <- cowplot::plot_grid(
  right_col_no_leg, legend_column, ncol = 2,
  rel_widths = c(1, 0.20), align = "v"
)

## Final figure: map left, right column (with legend) right
fig1 <- cowplot::plot_grid(
  p_A_map + labs(title = NULL), right_with_legend, ncol = 2,
  labels = c("A", NULL),
  label_fontface = "bold", label_size = 14,
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1,
  rel_widths = c(1.15, 1.60)
) + theme(plot.margin = margin(0,0,0,0), plot.background = element_blank())

## Export
ggsave(file.path(OUT, "Figure1.pdf"), fig1, width = 13, height = 9, units = "in", useDingbats = FALSE)
ggsave(file.path(OUT, "Figure1.png"), fig1, width = 13, height = 9, units = "in", dpi = 300)



# ===========================
# Save stats (optional)
# ===========================
readr::write_tsv(as.data.frame(adon_tab), file.path(OUT, "Figure1_PERMANOVA.tsv"))
readr::write_tsv(as.data.frame(art_tab),  file.path(OUT, "Figure1_aMPD_ART_two_way.tsv"))
readr::write_tsv(dunn_depth,              file.path(OUT, "Figure1_aMPD_Dunn_within_depth.tsv"))
readr::write_tsv(wilc_season,             file.path(OUT, "Figure1_aMPD_Wilcoxon_between_depth.tsv"))

message(sprintf("PERMANOVA (stratified by depth): R^2 = %.3f, P = %.4g", R2, pval))

# Extra exports
readr::write_tsv(as.data.frame(pw_df), file.path(OUT, "Figure1_pairwise_PERMANOVA.tsv"))
readr::write_tsv(adiv_pd, file.path(OUT, "Figure1_aMPD_per_sample.tsv"))
readr::write_tsv(dplyr::summarise(dplyr::group_by(adiv_pd, season, depth_m), n = dplyr::n(),
                                  median = median(aMPD), mean = mean(aMPD), .groups = "drop"),
                 file.path(OUT, "Figure1_aMPD_group_summary.tsv"))

# Season x depth PERMANOVA on Bray-Curtis (mirrors the Fig. 3D design).
# Main effects from the additive model (marginal tests); interaction from the full model.
# meta and dmat are already aligned to attr(dist_mat, "Labels") above.
stopifnot(identical(rownames(meta), rownames(dmat)))
set.seed(1)
adon_main <- vegan::adonis2(dmat ~ season + depth_m, data = meta, permutations = 999, by = "margin")
set.seed(1)
adon_int  <- vegan::adonis2(dmat ~ season * depth_m, data = meta, permutations = 999, by = "margin")
readr::write_tsv(tibble::rownames_to_column(as.data.frame(adon_main), "term"), file.path(OUT, "Figure1_PERMANOVA_season_depth_main.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(adon_int),  "term"), file.path(OUT, "Figure1_PERMANOVA_season_x_depth.tsv"))
