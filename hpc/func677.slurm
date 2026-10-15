#!/bin/bash
#SBATCH --job-name=func677
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/func677_%j.out
#SBATCH --error=logs/func677_%j.err
# Reproduce long_table.txt: Prodigal proteins for the 677 MRA MAGs, DIAMOND blastp best hit vs the
# metabolic marker DB (all_metabolic_genes.dmnd). Loose search; identity/coverage cutoffs applied downstream.
set -euo pipefail
module load Prodigal/prodigal-2.6.3 Diamond-2.1.23/diamond-2.1.23
PROJ="$HOME/comau_mag_coverm"; O="$PROJ/func677"
DB="$HOME/GMCF3029a_X-504149942/atmos_markers/all_metabolic_genes.dmnd"
cd "$PROJ/mags677"
ls *.fa | sed "s/\.fa$//" | xargs -P "$SLURM_CPUS_PER_TASK" -I{} sh -c \
  "prodigal -i {}.fa -a $O/faa/{}.faa -p single -q > /dev/null && sed -i \"s/^>/>{}~/\" $O/faa/{}.faa"
echo "faa files: $(ls $O/faa | wc -l)"
cat $O/faa/*.faa > $O/all677.faa
echo "proteins: $(grep -c ">" $O/all677.faa)"
diamond blastp -d "$DB" -q $O/all677.faa -o $O/all677_vs_metab.tsv --threads "$SLURM_CPUS_PER_TASK" \
  --sensitive --max-target-seqs 1 --evalue 1e-5 --id 20 \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore
echo "hits: $(wc -l < $O/all677_vs_metab.tsv)"
