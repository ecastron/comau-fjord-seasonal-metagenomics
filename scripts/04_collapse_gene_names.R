# 04_collapse_gene_names.R
# Collapses marker-set names to canonical gene symbols and writes the gene-by-MAG matrix used by Figures 3 to 5.
# Run from the repository root (see README.md).

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# ---- Inputs ----
in_long  <- "data/functional/long_table.txt"                      # columns: gene, hit, bin_id
out_wide <- "data/functional/gene_bin_mat.collapsed.csv"
out_long <- "data/functional/gene_bin_mat.collapsed.long.csv"
out_qc1  <- "data/functional/gene_name_collapses.tsv"          # original -> canonical (distinct pairs)
out_qc2  <- "results/functional/unmapped_gene_tokens.tsv"         # canonicals not in expected (should be near-empty)

# ---- Curated set (keeps dhbB, MenF) ----
expected <- c(
  "CydA","CyoA","CoxA","CcoN","FrdA","SdhA","NuoF","AtpA","ArsC","RdhA","ARO","FdhA",
  "PEPC","RbcL","AcnB","CA","Mcr","BtuCDF","BluB","CobA","CobC","CobU","CbiG","IsiB",
  "dmdA","dddL","DddP","dmdD","dmdB","mddA","dmoA","mtoX","Tmm","dorA",
  "AmoA","AmoB","AmoC","NADH","NrfA","NirK","NirS","NxrA","NorB","NapA","NarG","NosZ",
  "PsaA","PsbA","RHO","FCC","Sqr","SoxB","DsrA","PmoA","CoxL","NiFe","MDH",
  "dhbB","MenF"          # <-- kept, as requested
)
expected_ci <- setNames(expected, tolower(expected))
expected_lc <- tolower(expected)

# Generic tails we should avoid if a curated symbol exists
stop_tails <- c("sequence","protein","subunit","component","family","domain",
                "hydrogenase","reductase","oxidase","synthase","dehydrogenase",
                "cytochrome","nitrite","nitrate","aa3","cd1","large","small","type",
                "complex","chain","like","related")

is_numeric_str <- function(x) str_detect(x, "^[0-9]+$")

# ---- Canonicalizer (token scan, rightmost curated match; Tmmm -> Tmm) ----
canonicalize_vec <- function(lbl) {
  s <- str_remove(lbl, "(?i)_(filtered|unfiltered)$")
  s <- str_remove(s, "_[0-9]+$")
  toks_list <- str_split(s, "[-_]")
  
  out <- character(length(toks_list))
  for (i in seq_along(toks_list)) {
    ti <- toks_list[[i]]
    ti <- ti[ti != ""]
    ti <- ti[!is_numeric_str(ti)]
    # normalize close variants at token level
    ti <- ifelse(str_detect(ti, "(?i)^tmmm$"), "Tmm", ti)
    
    ti_lc <- tolower(ti)
    
    # 1) rightmost curated token
    match_idx <- which(ti_lc %in% expected_lc)
    if (length(match_idx) > 0) {
      chosen <- expected_ci[ ti_lc[max(match_idx)] ]
      out[i] <- unname(chosen)
      next
    }
    
    # 2) conservative synonyms on full label
    li <- lbl[i]
    if (str_detect(li, "(?i)rubisco.*rbc[l]") || str_detect(paste(ti, collapse="_"), "(?i)^rubisco$")) {
      out[i] <- "RbcL"; next
    }
    if (str_detect(li, "(?i)microbial.*rhodopsin") || str_detect(paste(ti, collapse="_"), "(?i)^rhodopsin$")) {
      out[i] <- "RHO"; next
    }
    
    # 3) fallback: last non-generic token (avoid tails)
    if (length(ti) == 0) { out[i] <- ""; next }
    idx <- length(ti)
    while (idx >= 1 && tolower(ti[idx]) %in% stop_tails) idx <- idx - 1
    last_tok <- if (idx >= 1) ti[idx] else tail(ti, 1)
    
    # casing normalization
    if (str_detect(last_tok, "(?i)^rbcl$")) last_tok <- "RbcL"
    if (str_detect(last_tok, "(?i)^rbcs$")) last_tok <- "RbcS"
    if (str_detect(last_tok, "(?i)^rho$"))  last_tok <- "RHO"
    if (str_detect(last_tok, "(?i)^nife$")) last_tok <- "NiFe"
    if (str_detect(last_tok, "(?i)^nadh$")) last_tok <- "NADH"
    if (str_detect(last_tok, "(?i)^mdh$"))  last_tok <- "MDH"
    if (str_detect(last_tok, "(?i)^fcc$"))  last_tok <- "FCC"
    
    # map to expected if exact (case-insensitive)
    lc <- tolower(last_tok)
    out[i] <- if (lc %in% names(expected_ci)) expected_ci[[lc]] else last_tok
  }
  out
}

# ---- Read tidy long table ----
df <- read_tsv(
  in_long,
  col_types = cols(gene = col_character(), hit = col_double(), bin_id = col_character())
)

# ---- Apply canonicalization ----
df <- df %>% mutate(gene_canonical = canonicalize_vec(gene))

# ---- Summarize gene × bin ----
gene_bin_long <- df %>%
  group_by(gene = gene_canonical, bin_id) %>%
  summarise(total_hits = sum(hit, na.rm = TRUE), .groups = "drop")

# ---- Wide matrix ----
gene_bin_wide <- gene_bin_long %>%
  pivot_wider(names_from = bin_id, values_from = total_hits, values_fill = 0) %>%
  arrange(gene)

# ---- Outputs ----
write_csv(gene_bin_wide, out_wide)
write_csv(gene_bin_long, out_long)

# ---- QC ----
qc_map <- df %>%
  distinct(gene, gene_canonical) %>%
  arrange(gene_canonical, gene)
write_tsv(qc_map, out_qc1)

unmapped <- qc_map %>%
  filter(!tolower(gene_canonical) %in% expected_lc) %>%
  count(gene_canonical, name = "n_labels") %>%
  arrange(desc(n_labels))
write_tsv(unmapped, out_qc2)


