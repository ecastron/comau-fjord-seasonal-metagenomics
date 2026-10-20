# Marker definitions

Functional potential is estimated from two published resources; no custom sequence database is used.

- `greening_markers.tsv`: the 33 marker sets counted from the Greening lab metabolic marker gene databases
  (Leung PM, Greening C. Compiled Greening lab metabolic marker gene databases. Monash University, 2021.
  <https://doi.org/10.26180/14431208.v1>), with the minimum amino-acid identity used for each (Greening lab thresholds
  for MAGs, as in Hutchinson TF et al. ISME J 2026;20:wrag020; query coverage at least 80% for all) and the number of
  reference sequences per set. `hpc/markerdb.slurm` downloads the databases from figshare and searches all 51 sets.
- `kofam_markers.tsv`: the KEGG KOs searched with KOfam profile HMMs (Aramaki T et al. Bioinformatics 2020;36:2251-2)
  and the marker each KO counts toward. K17224 (soxB) is searched only to cross-check the Greening SoxB calls.
- `kofam_thresholds.tsv`: the KEGG adaptive score thresholds (`ko_list` rows) used for these KOs.

Pathway assignment and per-marker counts are in Supplementary Table S3 (`results/supplementary/Supplementary_Table_S3_markers.xlsx`).
