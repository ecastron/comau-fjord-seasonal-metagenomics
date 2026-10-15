#!/bin/bash
#SBATCH --job-name=comau513_cov
#SBATCH --partition=basic
#SBATCH --array=0-89%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cov513_%A_%a.out
#SBATCH --error=logs/cov513_%A_%a.err
# Per-sample MAG abundance and detection for the 90 Comau metagenomes (PRJNA729490)
# against the 513 dereplicated MAGs used in Comau_ps90.RDS (MRA 2023, GenBank PRJNA729490).
# Competitive mapping across all 513 genomes, as in the MRA (CoverM, --min-read-aligned-percent 30).
# --min-covered-fraction 0 so no genome is zeroed; detection is decided downstream from
# covered_fraction vs. breadth expected from mean depth.
set -euo pipefail
PROJ="$HOME/comau_mag_coverm"
READS=/home/databases/comau/reads
MAP="$PROJ/sample_reads_map.tsv"   # sample_ID, SRR, library, read-file prefix (90 rows)

read -r SID SRR LIB RD < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MAP")
R1="$READS/${RD}_1.fastq"; R2="$READS/${RD}_2.fastq"
OUT="$PROJ/coverm513/${SID}.tsv"
[ -s "$OUT" ] && { echo "$SID already done"; exit 0; }
echo "[$(date +%F_%T)] $SID ($SRR, $LIB) on $(hostname)"

set +u; module purge; module load CoverM/Coverm; module load SAMTOOLS/samtools.1.20; set -u

coverm genome \
  --genome-fasta-directory "$PROJ/ref513" --genome-fasta-extension fa \
  --coupled "$R1" "$R2" \
  --mapper minimap2-sr \
  --min-read-aligned-percent 30 \
  --min-covered-fraction 0 \
  --methods count covered_fraction mean trimmed_mean relative_abundance \
  --threads "$SLURM_CPUS_PER_TASK" \
  --output-file "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "[$(date +%F_%T)] $SID done -> $OUT"
