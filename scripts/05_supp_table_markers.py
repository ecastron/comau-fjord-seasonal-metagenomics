# 05_supp_table_markers.py
# Supplementary Table S3: the 88 marker sets of the reference database, their source (Greening lab databases or UniProt), sequence counts, and pathway assignment.
# Run from the repository root (see README.md).

import re, csv, openpyxl
from openpyxl.styles import Font, Alignment
D = "data/functional/"
counts = [(int(n), m) for n, m in (l.split() for l in open(D + "metabolic_db_marker_counts.txt"))]
# marker set -> gene symbol, from the name-collapsing table of our run (covers every set with at least one hit)
sym = {}
for r in csv.DictReader(open(D + "gene_name_collapses.tsv"), delimiter="\t"):
    sym.setdefault(re.sub(r"_\d+$", "", r["gene"]), r["gene_canonical"])
pathways = {
 "Respiration": ["CydA","CyoA","CoxA","CcoN","NuoF","FrdA","SdhA","AtpA","FCC"],
 "Carbon core": ["RbcL","PEPC","AcnB","CA","CoxL","MDH"],
 "Cobalamin (B12)": ["BluB","BtuCDF","CobA","CobC","CobU","CbiG","IsiB"],
 "DMSP and methylated compounds": ["dmdA","dmdB","dmdD","dddL","DddP","mddA","dmoA","Tmm","mtoX","dorA"],
 "Nitrogen cycle": ["AmoA","AmoB","AmoC","NxrA","NapA","NarG","NirK","NirS","NrfA","NorB","NosZ"],
 "Sulfur cycle": ["Sqr","SoxB","DsrA"],
 "Methane cycling": ["PmoA","Mcr"],
 "Phototrophy": ["PsaA","PsbA","RHO"],
 "Hydrogenases": ["NiFe"],
 "Detoxification and other": ["ArsC","RdhA","ARO","FdhA"],
}
g2p = {g.lower(): p for p, gs in pathways.items() for g in gs}
def is_greening(m):
    return m == "Fe_hydrogenase" or (re.search(r"[A-Za-z]{3,}_[A-Za-z]", m) and m != "Nitrate_reductase_NADH")
rows = []
for n, m in counts:
    s = sym.get(m) or re.split(r"_", m.replace("_sequences", "").replace("_sequence", ""))[-1]
    p = g2p.get(s.lower(), "Not used in pathway analyses")
    src = "Greening lab compiled databases (Leung and Greening 2021)" if is_greening(m) else "UniProt, compiled for this study"
    rows.append((m, s, n, src, p))
assert sum(r[3].startswith("Greening") for r in rows) == 51 and len(rows) == 88
rows.sort(key=lambda r: (r[3], r[0].lower()))
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Supplementary Table 1"
ws.append(["Supplementary Table 1. Marker sets in the functional reference database (208,944 protein sequences, 88 markers). "
           "Malonyl-CoA reductase (Mcr) is a carbon-fixation enzyme and is distinct from methyl-coenzyme M reductase (McrA)."])
ws.append([])
ws.append(["Marker set (database name)", "Gene symbol in analysis", "Reference sequences", "Source", "Pathway in this study"])
for c in ws[3]: c.font = Font(bold=True)
for r in rows: ws.append(list(r))
for col, w in zip("ABCDE", (58, 16, 12, 52, 30)): ws.column_dimensions[col].width = w
ws["A1"].alignment = Alignment(wrap_text=True); ws.merge_cells("A1:E1"); ws.row_dimensions[1].height = 32
wb.save("results/supplementary/Supplementary_Table_S3_markers.xlsx")
print("rows:", len(rows), "| total sequences:", sum(r[2] for r in rows))
print([r[:2] + r[4:] for r in rows if r[4] == "Methane cycling"])
