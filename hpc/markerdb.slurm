#!/bin/bash
#SBATCH --job-name=markerdb
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=logs/markerdb_%j.out
# Marker searches on the 677 MAG proteomes (all677.faa from func677.slurm):
# (1) DIAMOND blastp best hit against all 51 Greening lab metabolic marker gene sets
#     (https://doi.org/10.26180/14431208.v1); identity and coverage thresholds are applied in scripts/03_make_gene_table.py.
# (2) hmmsearch with the KEGG KOfam profiles listed in database/kofam_markers.tsv; KEGG thresholds applied downstream.
# Outputs go to data/functional/ (gzip all677_vs_greening.tsv and all677_vs_kos.tbl).
set -euo pipefail
module load Diamond-2.1.23/diamond-2.1.23 Hmmer/hmmer-3.4
Q="$HOME/comau_mag_coverm/func677/all677.faa"
KOFAM="/path/to/kofam_database"   # KofamScan database: profiles/ and ko_list (ftp://ftp.genome.jp/pub/db/kofam/)
REPO="$HOME/comau-fjord-seasonal-metagenomics"

curl -sSL -o greening.zip https://ndownloader.figshare.com/files/27599390
unzip -q -o greening.zip -d greening
for f in greening/*/*.fasta; do
  m=$(basename "$f" .fasta | tr ' ' '_')
  awk -v m="$m" '/^>/{n++; print ">" m "|" n; next} {print}' "$f"
done > greening_all.faa
diamond makedb --in greening_all.faa -d greening_all --threads "$SLURM_CPUS_PER_TASK"
diamond blastp -d greening_all -q "$Q" -o all677_vs_greening.tsv --threads "$SLURM_CPUS_PER_TASK" --sensitive \
  --max-target-seqs 1 --evalue 1e-5 --id 20 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore

KOS=$(tail -n +2 "$REPO/database/kofam_markers.tsv" | cut -f1)
for k in $KOS; do cat "$KOFAM/profiles/$k.hmm"; done > kos.hmm
hmmsearch --cpu "$SLURM_CPUS_PER_TASK" --noali -E 1e-5 --tblout all677_vs_kos.tbl kos.hmm "$Q" > /dev/null
