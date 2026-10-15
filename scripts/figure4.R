# figure4.R
# Figure 4: family-level phylogeny with abundance-weighted pathway potential.
# Run from the repository root (see README.md).

source("utils.R")
OUT <- "results/functional"

cran_pkgs <- c(
  "dplyr", "tidyr", "tibble",
  "ggplot2", "ggnewscale", "scales"
)
bioc_pkgs <- c("phyloseq", "ggtree", "treeio", "ggtreeExtra")

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

## ---- 0) Load shared objects from Figure 3 ------------------------------------
ps <- load_ps("data/Comau_ps90_cov10_r232.RDS")
if (!exists("A_rel_use")) {
  fig3 <- readRDS(file.path(OUT, "fig3_shared_objects.RDS"))
  list2env(fig3, envir = .GlobalEnv)
}
stopifnot(exists("A_rel_use"), exists("H_use"), exists("pathway_dict"))
tx <- as.data.frame(tax_table(ps))
stopifnot(all(c("Family","Phylum") %in% colnames(tx)))

## Use Figure 3 pathways (five main + Phototrophy)
fig3_paths <- unique(c(
  if (exists("main_pathways")) main_pathways else c("Methane_redox","Carbon_core","Sulfur_cycle","B12_cobalamin","Nitrogen_cycle"),
  "Phototrophy"
))

## Helper: retain only genes available in H_use
gene_cols_H <- rownames(H_use)
avail <- function(g) intersect(g, gene_cols_H)

## ---- 1) Per-MAG pathway hits from H_use -------------------------------------
mag_ids <- colnames(H_use)
stopifnot(all(mag_ids %in% colnames(A_rel_use)))

mag_path_hits <- lapply(fig3_paths, function(pw){
  genes <- avail(pathway_dict[[pw]])
  if (length(genes) == 0) return(rep(0, length(mag_ids)))
  colSums(H_use[genes, , drop = FALSE], na.rm = TRUE)
})
mag_path_hits <- do.call(cbind, mag_path_hits)
colnames(mag_path_hits) <- fig3_paths
rownames(mag_path_hits) <- mag_ids

## ---- 2) Weight by mean MAG relative abundance -------------------------------
mag_mean_abund <- colMeans(A_rel_use, na.rm = TRUE)
mag_path_potential <- sweep(mag_path_hits, 1, mag_mean_abund[rownames(mag_path_hits)], `*`)

## ---- 3) Aggregate MAG -> Family ----------------------------------------------
mag_tax <- tx %>% rownames_to_column("MAG_id") %>% select(MAG_id, Family, Phylum)
mag_tax <- mag_tax %>% filter(MAG_id %in% rownames(mag_path_potential))

fam_path_potential <- mag_tax %>%
  select(MAG_id, Family) %>%
  left_join(
    mag_path_potential %>% as.data.frame() %>% rownames_to_column("MAG_id"),
    by = "MAG_id"
  ) %>%
  group_by(Family) %>%
  summarise(across(all_of(fig3_paths), ~sum(.x, na.rm = TRUE)), .groups = "drop")

if (nrow(fam_path_potential) == 0L ||
    all(as.matrix(fam_path_potential[, fig3_paths, drop = FALSE]) == 0)) {
  stop("Family x Pathway table is empty or all zeros. Check H_use/pathway_dict coverage.")
}

## ---- 4) Collapse phylogeny to Family & rename tips --------------------------
ps_fam <- tax_glom(ps, taxrank = "Family", NArm = TRUE)
fam_vec <- as.character(tax_table(ps_fam)[, "Family"])
fam_vec[is.na(fam_vec) | fam_vec == ""] <- "Unassigned_Family"
taxa_names(ps_fam) <- make.unique(fam_vec, sep = "__")

td        <- as.treedata(phy_tree(ps_fam))
tip_order <- td@phylo$tip.label
tip_map <- tibble(
  Tip    = taxa_names(ps_fam),
  Family = as.character(tax_table(ps_fam)[, "Family"]),
  Phylum = as.character(tax_table(ps_fam)[, "Phylum"])
) %>%
  mutate(
    Family = if_else(is.na(Family) | Family == "", "Unassigned_Family", Family),
    Phylum = if_else(is.na(Phylum) | Phylum == "", "Other", Phylum)
  )

## ---- 5) Build heat matrix ---------------------------------------------------
tip_pathways <- tip_map %>% left_join(fam_path_potential, by = "Family")
heat_mat <- tip_pathways %>%
  arrange(match(Tip, tip_order)) %>%
  select(all_of(c("Tip", fig3_paths))) %>%
  column_to_rownames("Tip")

## ---- 6) Safe z-score & robust limits ----------------------------------------
safe_z <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- 0
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
heat_scaled <- as.data.frame(apply(heat_mat, 2, safe_z))
# Human-readable column names for display
pretty_pw <- c(DMSP_Methyl = "DMSP", Carbon_core = "Carbon core", Sulfur_cycle = "Sulfur cycle", B12_cobalamin = "Cobalamin (B12)", Nitrogen_cycle = "Nitrogen cycle", Phototrophy = "Phototrophy", Methane_redox = "Methane cycling")
colnames(heat_scaled) <- ifelse(colnames(heat_scaled) %in% names(pretty_pw), pretty_pw[colnames(heat_scaled)], gsub("_", " ", colnames(heat_scaled)))

z_abs    <- abs(unlist(heat_scaled))
z_lim    <- max(1.5, as.numeric(stats::quantile(z_abs, 0.98, na.rm = TRUE)))
z_breaks <- scales::pretty_breaks(n = 5)(c(-z_lim, z_lim))

## ---- 7) Plot ----------------------------------------------------------------
# Reshape heat_scaled to long format for geom_fruit
heat_long <- heat_scaled %>%
  rownames_to_column("Tip") %>%
  pivot_longer(-Tip, names_to = "Pathway", values_to = "z_score") %>%
  mutate(Pathway = factor(Pathway, levels = colnames(heat_scaled)))

p_out <- ggtree(td, layout = "circular", linewidth = 0.3) %<+% tip_map +
  geom_tippoint(aes(color = Phylum), size = 2.5, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = pal_phylum, na.value = "#BFBFBF", name = "Phylum") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_tiplab2(aes(label = Family),
               size = 2.2, offset = 0.07, align = FALSE, linesize = 0.3) +
  new_scale_fill() +
  geom_fruit(
    data    = heat_long,
    geom    = geom_tile,
    mapping = aes(y = Tip, x = Pathway, fill = z_score),
    offset  = 0.15,   # small gap after tips (fraction of tree radius)
    pwidth  = 0.58,   # width of heatmap ring relative to tree radius
    axis.params = list(
      axis       = "x",
      text.angle = 45,
      text.size  = 3,
      hjust      = 0,
      vjust      = 0.5,
      line.size  = 0.4
    ),
    grid.params = list(linetype = "blank")
  ) +
  scale_fill_gradient2(
    low      = "#313695", mid = "white", high = "#A50026",
    midpoint = 0,
    limits   = c(-z_lim, z_lim), breaks = z_breaks,
    oob      = scales::squish,
    name     = "Pathway enrichment\n(family z-score)"
  ) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(5, 6, 5, 6)
  )

## ---- 8) Save ----------------------------------------------------------------
ggsave(file.path(OUT, "Figure4.png"), p_out,
       width = 300, height = 270, units = "mm", dpi = 300, bg = "white")
if (requireNamespace("magick", quietly = TRUE)) {   # crop the empty margin around the circular layout
  magick::image_write(magick::image_border(magick::image_trim(magick::image_read(file.path(OUT, "Figure4.png"))), "white", "40x40"), file.path(OUT, "Figure4.png"))
}
ggsave(file.path(OUT, "Figure4.pdf"), p_out,
       width = 300, height = 270, units = "mm")


# Extra exports
readr::write_tsv(tibble::rownames_to_column(as.data.frame(heat_mat), "Tip"), file.path(OUT, "Figure4_family_pathway_potential_raw.tsv"))
readr::write_tsv(tibble::rownames_to_column(as.data.frame(heat_scaled), "Tip"), file.path(OUT, "Figure4_family_pathway_z.tsv"))
readr::write_tsv(tibble::tibble(z_color_limit = z_lim), file.path(OUT, "Figure4_color_limit.tsv"))
