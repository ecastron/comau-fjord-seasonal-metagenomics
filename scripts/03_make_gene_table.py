# 03_make_gene_table.py
# Builds the gene-by-MAG hit table used by Figures 3 to 5 from the marker searches in hpc/markerdb.slurm:
#   (1) DIAMOND blastp best hit of each MAG protein against all 51 Greening lab marker sets
#       (https://doi.org/10.26180/14431208.v1), kept at >= 80% query coverage and the marker-specific identity
#       threshold in database/greening_markers.tsv (Greening lab thresholds for MAGs; Hutchinson et al. 2026 ISME J).
#       All 51 sets are searched so that each protein lands on its closest marker; only the sets listed in
#       greening_markers.tsv are counted (the Greening PsbA and ArsC sets also recruit pufM and ArsC-family proteins,
#       so psbA, pufM and arsC come from KOfam).
#   (2) hmmsearch with the KEGG KOfam profiles in database/kofam_markers.tsv, kept above the KEGG adaptive threshold
#       (full-sequence or best-domain score, database/kofam_thresholds.tsv). A protein passing several KOs is assigned
#       to the KO with the highest score/threshold ratio.
# BtuCDF counts complete transporters: per MAG, the minimum of the btuC, btuD and btuF counts.
# Run from the repository root (see README.md).
import csv, collections, gzip

D = "data/functional"
greening = {r["greening_set"]: (r["marker"], float(r["min_identity_pct"]))
            for r in csv.DictReader(open("database/greening_markers.tsv"), delimiter="\t")}
kofam = {r["ko"]: r["marker"] for r in csv.DictReader(open("database/kofam_markers.tsv"), delimiter="\t")
         if r["marker"] != "SoxB_check"}
thr = {r["knum"]: (float(r["threshold"]), r["score_type"])
       for r in csv.DictReader(open("database/kofam_thresholds.tsv"), delimiter="\t")}

calls = {}
with gzip.open(f"{D}/all677_vs_greening.tsv.gz", "rt") as f:
    for q, s, pid, ln, ql, sl, qc, sc, ev, bs in csv.reader(f, delimiter="\t"):
        hit = greening.get(s.split("|")[0])
        if hit and float(qc) >= 80 and float(pid) >= hit[1]:
            calls[q] = hit[0]

best = {}
with gzip.open(f"{D}/all677_vs_kos.tbl.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.split()
        prot, ko, full, dom = x[0], x[2], float(x[5]), float(x[8])
        if ko not in kofam:
            continue
        t, kind = thr[ko]
        ratio = (full if kind == "full" else dom) / t
        if ratio >= 1 and ratio > best.get(prot, (0, None))[0]:
            best[prot] = (ratio, ko)
overlap = sum(p in calls for p in best)
for prot, (_, ko) in best.items():
    calls.setdefault(prot, kofam[ko])

counts = collections.defaultdict(collections.Counter)
for prot, m in calls.items():
    counts[m][prot.split("~")[0]] += 1
mags = sorted(l.split("\t")[0] for l in open("data/coverm/mag677_len.tsv"))
assert len(mags) == 677
for mag in mags:
    n = min(counts["btuC"][mag], counts["btuD"][mag], counts["btuF"][mag])
    if n:
        counts["BtuCDF"][mag] = n
for sub in ("btuC", "btuD", "btuF"):
    counts.pop(sub, None)

with open(f"{D}/gene_bin_mat.csv", "w", newline="") as o:
    wr = csv.writer(o)
    wr.writerow(["gene"] + mags)
    for m in sorted(counts):
        wr.writerow([m] + [counts[m][mag] for mag in mags])
with open(f"{D}/protein_marker_calls.tsv", "w") as o:
    o.write("protein\tmarker\n")
    for prot in sorted(calls):
        o.write(f"{prot}\t{calls[prot]}\n")
print("MAGs:", len(mags), "| markers with hits:", len(counts), "| proteins:", len(calls), "| Greening/KOfam overlaps:", overlap)
for m in sorted(counts):
    print(f"  {m:8s} {sum(counts[m].values()):5d} hits in {sum(v > 0 for v in counts[m].values()):3d} MAGs")
