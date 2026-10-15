#!/bin/bash
# Reruns the analysis chain from the repository root and records R and package versions.
# Usage: bash run_all.sh            (statistics, figures and tables from the processed data)
#        bash run_all.sh --rebuild  (also rebuilds the phyloseq objects and the gene table first)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/community results/functional results/supplementary logs
if [[ "${1:-}" == "--rebuild" ]]; then
  Rscript scripts/00_coverm_detection_summary.R > logs/00.log 2>&1
  Rscript scripts/01_build_ps_cov10.R           > logs/01.log 2>&1
  Rscript scripts/02_build_ps_r232.R            > logs/02.log 2>&1
  python3 scripts/03_make_gene_table.py         > logs/03.log 2>&1
  Rscript scripts/04_collapse_gene_names.R      > logs/04.log 2>&1
  echo "rebuild done"
fi
for s in figure1 figure2 suppfig_environment suppfig_gam_temperature_salinity \
         figure3 figure4 figure5 supptable_depth_patterns depth_contrast_by_season make_supp_tables; do
  start=$(date +%s); Rscript "scripts/$s.R" > "logs/$s.log" 2>&1; echo "$s done ($(( $(date +%s) - start )) s)"
done
python3 scripts/05_supp_table_markers.py > logs/05.log 2>&1
Rscript -e 'pk <- c("phyloseq","vegan","picante","mgcv","ARTool","ecodist","FSA","car","ape","ggplot2","ggtree","ggtreeExtra","patchwork","dplyr","tidyr","geosphere","rstatix","magick")
v <- sapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA))
writeLines(c(R.version.string, paste(names(v), v, sep = " ")), "logs/package_versions_this_run.txt")'
