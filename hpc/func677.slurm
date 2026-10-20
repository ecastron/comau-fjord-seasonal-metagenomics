#!/bin/bash
#SBATCH --job-name=func677
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/func677_%j.out
#SBATCH --error=logs/func677_%j.err
# Gene prediction (Prodigal) for the 677 MAGs; the proteins (all677.faa) are searched against the marker
# databases in markerdb.slurm.
set -euo pipefail
module load Prodigal/prodigal-2.6.3
PROJ="$HOME/comau_mag_coverm"; O="$PROJ/func677"
cd "$PROJ/mags677"
ls *.fa | sed "s/\.fa$//" | xargs -P "$SLURM_CPUS_PER_TASK" -I{} sh -c \
  "prodigal -i {}.fa -a $O/faa/{}.faa -p single -q > /dev/null && sed -i \"s/^>/>{}~/\" $O/faa/{}.faa"
echo "faa files: $(ls $O/faa | wc -l)"
cat $O/faa/*.faa > $O/all677.faa
echo "proteins: $(grep -c ">" $O/all677.faa)"
