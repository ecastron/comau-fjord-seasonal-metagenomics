# 03_make_gene_table.py
# Builds the gene-by-MAG hit table from DIAMOND blastp best hits of Prodigal proteins (677 MAGs) against the marker database, keeping hits with at least 50% identity over at least 80% of the protein.
# Run from the repository root (see README.md).

import gzip, csv, collections
c = collections.Counter()
with gzip.open("data/functional/all677_vs_metab.tsv.gz", "rt") as f:
    for q, s, pid, ln, ql, sl, qc, sc, ev, bs in csv.reader(f, delimiter="\t"):
        if float(pid) >= 50 and float(qc) >= 80:
            c[(s, q.split("~")[0])] += 1
with open("data/functional/long_table.txt", "w") as o:
    o.write("gene\thit\tbin_id\n")
    for (g, b), n in sorted(c.items(), key=lambda x: (x[0][1], x[0][0])):
        o.write(f"{g}\t{n}\t{b}\n")
print("rows:", len(c), "| total hits:", sum(c.values()), "| MAGs:", len({b for _, b in c}))
