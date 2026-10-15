# Seasonal structuring of microbial communities and functional potential in a Patagonian fjord
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23401703.svg)](https://doi.org/10.5281/zenodo.23401703)

Code, processed data and the marker protein database for:

> Guajardo-Leiva S, Escobedo VM, Galbán-Malagón C, Berríos-Farías V, Moreno-Pino M, Molina-Trincado V, Castro-Nallar E.
> Seasonal structuring of microbial communities and functional potential in a Patagonian fjord. *In preparation*.

The study follows 90 shotgun metagenomes from the Comau Fjord (northern Patagonia, Chile), sampled at 5 and 20 m
in Summer, Fall and Winter, through 513 metagenome-assembled genomes (MAGs). It tests how season, depth and
position along the fjord structure community composition, phylogenetic structure and functional potential.

## Layout

```
data/
  Comau_ps90.RDS                 phyloseq object: 90 samples x 513 MAGs, read counts, original taxonomy (GTDB R202)
  Comau_ps90_cov10.RDS           same, after the 10% genome-coverage detection filter
  Comau_ps90_cov10_r232.RDS      same, with GTDB R232 taxonomy: the object used by every analysis
  metadata.tsv                   sample metadata and environmental variables
  coverm/                        CoverM per-sample coverage of the 513 MAGs; MAG GenBank accessions and lengths
  functional/                    DIAMOND hits of MAG proteins vs the marker database, gene-by-MAG matrices
  gtdbtk_r232/                   GTDB-Tk v2.7.2 (GTDB R232) classification summaries
database/
  all_metabolic_genes.fasta.gz   marker protein database (88 marker sets, 208,944 sequences)
scripts/                         R and Python analysis scripts (run from the repository root)
hpc/                             Slurm scripts used on the cluster (CoverM, Prodigal + DIAMOND, GTDB-Tk)
results/                         outputs written by the scripts (figures, statistics tables, supplementary tables)
utils.R                          shared helpers (load_ps(), palettes)
run_all.sh                       reruns the full chain
```

## Reproducing the analyses

Requirements: R 4.5 and the packages listed in `package_versions.txt` (Bioconductor: `phyloseq`, `ggtree`,
`ggtreeExtra`; CRAN: the rest, plus `sf`, `readr`, `rstatix`), and Python 3 with `openpyxl` for the marker table.

The map in Figure 1 needs the GADM administrative boundaries of Chile, which are not redistributed here.
Download `CHL_adm.zip` (GADM v2.8, level 2) from <https://biogeo.ucdavis.edu/data/diva/adm/CHL_adm.zip> and unzip it
so that `data/CHL_adm_shp/CHL_adm2.shp` exists.

Then, from the repository root:

```bash
bash run_all.sh            # statistics, figures and tables from the processed data
bash run_all.sh --rebuild  # also rebuilds the filtered phyloseq objects and the gene table first
```

Logs go to `logs/`. All stochastic steps (null models, permutations, jitter) use `set.seed(100)`; the null models
run in parallel with `RNGkind("L'Ecuyer-CMRG")`, so results do not depend on the number of cores.

| Script | Output |
|---|---|
| `00_coverm_detection_summary.R` | per-sample CoverM summary and detection-rule comparison |
| `01_build_ps_cov10.R` | 10% genome-coverage detection filter |
| `02_build_ps_r232.R` | GTDB R232 taxonomy |
| `03_make_gene_table.py`, `04_collapse_gene_names.R` | gene-by-MAG matrix (best hit, at least 50% identity over at least 80% of the protein) |
| `05_supp_table_markers.py` | Supplementary Table S3 (marker sets) |
| `figure1.R` to `figure5.R` | Figures 1 to 5 and their statistics |
| `suppfig_environment.R` | Supplementary Figures S1 and S3 |
| `suppfig_gam_temperature_salinity.R` | Supplementary Figures S4 and S5 |
| `figure3.R` | also Supplementary Figure S2 |
| `supptable_depth_patterns.R` | Supplementary Table S4 |
| `depth_contrast_by_season.R` | 5 m versus 20 m contrast within each season |
| `make_supp_tables.R` | Supplementary Tables S1 and S2 |

## Upstream processing (cluster)

The scripts in `hpc/` document the steps run on the cluster from the public reads and MAGs:

- `coverm_513.slurm`: competitive mapping of each metagenome to the 513 dereplicated MAGs (CoverM), reporting read
  counts and covered fraction per MAG.
- `func677.slurm`: gene prediction (Prodigal) on the 677 MAGs and DIAMOND blastp against the marker database.
- `gtdbtk_513_r232.slurm`: GTDB-Tk v2.7.2 classification against GTDB R232.

Paths in these scripts refer to the cluster and need editing before reuse.

## Data sources

- Metagenome reads: NCBI BioProject [PRJNA729490](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA729490); SRA run
  and BioSample accessions per sample are in Supplementary Table S1 and `data/metadata.tsv`.
- MAGs: NCBI GenBank; accessions for the 677 MAGs are in `data/coverm/mag677_download.tsv` and Supplementary Table S2.

## Marker database

`database/all_metabolic_genes.fasta.gz` combines 51 marker sets from the Greening lab metabolic marker gene databases
(Leung PM, Greening C. Compiled Greening lab metabolic marker gene databases. Monash University, 2021.
<https://doi.org/10.26180/14431208.v1>) with 37 sets retrieved from UniProt (<https://www.uniprot.org>).
Sequence counts per set are in `data/functional/metabolic_db_marker_counts.txt`; source and pathway per set are in
Supplementary Table S3. Please cite both sources when reusing the database.

## License

Code: MIT (`LICENSE`). Data, database and results: CC BY 4.0 (`LICENSE-DATA`).

## Citation

Please cite the article above and this archive (see `CITATION.cff`): https://doi.org/10.5281/zenodo.23401703.
