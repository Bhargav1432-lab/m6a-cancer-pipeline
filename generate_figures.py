import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("figures", exist_ok=True)
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

# -------------------------------------------------------------
# 1. FIG 1: Genomic Region Distribution (3'UTR Validation)
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
regions = ["3'UTR", "CDS", "5'UTR", "Unmapped"]
percentages = [61, 34, 4, 1]
colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#d3d3d3']

bars = ax.bar(regions, percentages, color=colors, edgecolor='black', linewidth=0.8, width=0.55)
ax.set_ylabel('Percentage of Sites (%)', fontsize=12, fontweight='bold')
ax.set_title("m6A Topological Distribution (Top 100 Genes)\nValidation of 3'UTR / Stop-Codon Enrichment", fontsize=13, fontweight='bold', pad=12)
ax.set_ylim(0, 75)

for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{yval}%", ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig("figures/fig1_region_distribution.png")
plt.close()

# -------------------------------------------------------------
# 2. FIG 2: Asymmetric Modification Distribution (Pie/Donut)
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
counts = [1362, 383]
labels = ['Lung Cancer (A549)\nHypermethylated (78.1%)', 'Liver Cancer (HepG2)\nHypermethylated (21.9%)']
colors = ['#d62728', '#2ca02c']

wedges, texts, autotexts = ax.pie(counts, labels=labels, autopct='%1.1f%%', startangle=140,
                                  colors=colors, wedgeprops=dict(width=0.4, edgecolor='black', linewidth=1),
                                  textprops=dict(fontsize=11, fontweight='bold'))
for at in autotexts:
    at.set_fontsize(11)
    at.set_color('white')

ax.set_title("Transcriptome Modification Asymmetry\n(1,745 Cross-Validated Genes)", fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig("figures/fig2_modification_direction.png")
plt.close()

# -------------------------------------------------------------
# 3. FIG 3: Pathway Enrichment Top 10
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)
pathways = [
    "Protein processing in ER",
    "Mitochondrial Gene Expression",
    "Mitochondrial Translation",
    "RNA transport",
    "Translation initiation",
    "Ribosome Biogenesis",
    "Proteasome Ubiquitination",
    "mRNA Processing",
    "Ubiquitin Protein Catabolism",
    "Cell cycle"
]
neg_log_p = [15.04, 11.71, 10.41, 9.18, 8.03, 7.77, 7.67, 7.28, 7.17, 6.99]
y_pos = np.arange(len(pathways))

bars = ax.barh(y_pos, neg_log_p, color='#3b528b', edgecolor='black', linewidth=0.8, height=0.65)
ax.set_yticks(y_pos)
ax.set_yticklabels(pathways, fontsize=10, fontweight='bold')
ax.invert_yaxis()
ax.set_xlabel("-log10(Adjusted P-value)", fontsize=12, fontweight='bold')
ax.set_title("Top 10 Enriched Pathways in Differentially Modified Transcripts", fontsize=13, fontweight='bold', pad=12)

for bar in bars:
    w = bar.get_width()
    ax.text(w + 0.2, bar.get_y() + bar.get_height()/2.0, f"{w:.1f}", ha='left', va='center', fontsize=9, fontweight='bold')

ax.set_xlim(0, 17)
plt.tight_layout()
plt.savefig("figures/fig3_pathway_enrichment.png")
plt.close()

# -------------------------------------------------------------
# 4. FIG 4: TCGA Hazard Ratios (Forest Plot)
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5), dpi=300)
genes = ["CD24", "UBE3C", "AHCY", "RPS6KA3", "FADD"]
hrs = [1.48, 1.42, 1.39, 1.35, 1.37]
ci_lower = [1.11, 1.06, 1.04, 1.01, 1.03]
ci_upper = [1.97, 1.90, 1.86, 1.81, 1.83]

y_pos = np.arange(len(genes))
errors = [
    np.array(hrs) - np.array(ci_lower),
    np.array(ci_upper) - np.array(hrs)
]

ax.errorbar(hrs, y_pos, xerr=errors, fmt='o', color='#b2182b', ecolor='#2166ac',
            elinewidth=2.5, capsize=5, markersize=8, markeredgecolor='black')
ax.axvline(1.0, color='gray', linestyle='--', linewidth=1.2, label='Neutral Risk (HR=1.0)')

ax.set_yticks(y_pos)
ax.set_yticklabels(genes, fontsize=11, fontweight='bold')
ax.invert_yaxis()
ax.set_xlabel("Hazard Ratio (95% CI) for Overall Survival", fontsize=11, fontweight='bold')
ax.set_title("TCGA-LUAD Clinical Validation: Lung Cancer Patient Risk", fontsize=13, fontweight='bold', pad=12)
ax.set_xlim(0.8, 2.2)

for i, hr in enumerate(hrs):
    ax.text(ci_upper[i] + 0.04, i, f"HR={hr:.2f} (p<0.05)", va='center', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig("figures/fig4_tcga_hazard_ratios.png")
plt.close()

print("All 4 publication figures generated successfully in figures/ folder.")
