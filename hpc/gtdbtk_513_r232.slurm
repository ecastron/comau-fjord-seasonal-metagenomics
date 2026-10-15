#!/bin/bash
#SBATCH --job-name=gtdbtk513
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gtdbtk513_%j.out
#SBATCH --error=logs/gtdbtk513_%j.err
# Reclassify the 513 MAGs in Comau_ps90.RDS with the latest GTDB (GTDB-Tk 2.7.2, release R232); the MRA used R202.
set -eo pipefail   # no -u: the module activation script reads unset JAVA_HOME
module load GTDB-tk/2.7.2
PROJ="$HOME/comau_mag_coverm"
gtdbtk --version
gtdbtk classify_wf --genome_dir "$PROJ/ref513" --extension fa --out_dir "$PROJ/gtdbtk_513_r232" \
  --cpus "$SLURM_CPUS_PER_TASK"
ls -la "$PROJ/gtdbtk_513_r232"/*summary.tsv
