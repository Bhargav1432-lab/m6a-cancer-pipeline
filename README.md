# Direct RNA Methylation Mapping in Cancer Transcriptomes

[![Python](https://img.shields.io/badge/python-3.9+-blue)](https://python.org)
[![Snakemake](https://img.shields.io/badge/snakemake-7.0+-green)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Team 12 | DNA Sequencing and Technologies**  
> K. Siva Pranav | K. Mohana Venkata Krishna | P.V. Bhargav Reddy

---

## Project Overview

A fully computational pipeline that applies two independent AI tools — **m6Anet** and **xPore** — to Oxford Nanopore Technologies (ONT) direct RNA sequencing data from two cancer cell lines to map differential N⁶-methyladenosine (m6A) RNA modifications at single-molecule resolution.

**Cell lines compared:**
- `A549` — Lung adenocarcinoma (ATCC CCL-185)
- `HepG2` — Hepatocellular carcinoma (ATCC HB-8065)

**Data source:** [SG-NEx Project](https://github.com/GoekeLab/sg-nex-data) (Singapore Nanopore Expression Project, Genome Institute of Singapore)

---

## Key Results

| Metric | Value |
|--------|-------|
| m6A positions evaluated (m6Anet) | 187,022 |
| Differentially modified sites (xPore) | 138,512 |
| Cross-validated genes (both tools) | **1,745** |
| Expected by random chance | ~944 |
| Enrichment over chance | **1.85× (p < 0.001)** |
| Genes hypermethylated in lung cancer | 1,362 (78%) |
| 3'UTR site localisation | **61%** (consistent with known m6A biology) |
| Strongest pathway | Protein processing in ER (padj = 9.2 × 10⁻¹⁶) |

---

## Repository Structure

```
m6a_cancer_pipeline/
├── data/
│   ├── A549_rep6/dataprep/          # m6Anet dataprep output for A549
│   └── HepG2_rep6/dataprep/         # m6Anet dataprep output for HepG2
├── results/
│   ├── m6anet/                      # m6Anet site probabilities
│   ├── xpore/                       # xPore differential modification output
│   ├── cross_validation/            # Cross-validated gene lists
│   ├── pathway_enrichment/          # GO and KEGG enrichment results
│   └── tcga_validation/             # TCGA survival correlation data
├── pipeline/
│   ├── Snakefile                    # Main Snakemake workflow
│   ├── config.yaml                  # Pipeline configuration
│   └── envs/                        # Conda environment files
├── notebooks/
│   ├── 01_m6anet_analysis.ipynb     # m6Anet results exploration
│   ├── 02_xpore_analysis.ipynb      # xPore differential analysis
│   ├── 03_cross_validation.ipynb    # Cross-tool validation
│   ├── 04_pathway_enrichment.ipynb  # GO/KEGG enrichment
│   └── 05_tcga_validation.ipynb     # TCGA survival analysis
├── docs/
│   └── methods.md                   # Detailed methods documentation
├── tcga_query.py                    # TCGA cBioPortal API script
├── methods_comparison.py            # m6Anet vs xPore comparison
├── cross_validated_genes.csv        # Final 1,745 cross-validated genes
├── cross_validated_genes_annotated.csv
├── final_ranked_results.csv         # All genes ranked by significance
├── top100_with_regions.csv          # Top 100 genes with genomic regions
├── pathway_enrichment_results.csv   # Full pathway enrichment results
├── top20_pathway_enrichment.csv     # Top 20 enriched pathways
├── region_validation_summary.txt    # 3'UTR/CDS/5'UTR distribution
└── README.md
```

---

## Pipeline Overview

```
Raw Nanopore Signal (fast5/pod5)
        │
        ▼
Basecalling (Dorado GPU)
        │
        ▼
Alignment (minimap2 → reference transcriptome)
        │
        ▼
Signal Realignment (f5c eventalign)
        │
        ├─────────────────────────────────┐
        ▼                                 ▼
  m6Anet dataprep                   xPore dataprep
  m6Anet inference                  xPore diffmod
  (site probabilities,              (differential mod rates,
   stoichiometry)                    p-values)
        │                                 │
        └──────────┬──────────────────────┘
                   ▼
         Cross-Validation
         (hypergeometric test)
                   │
                   ▼
         Annotation & Enrichment
         (MyGene.info, Enrichr API)
                   │
                   ▼
         TCGA Validation
         (cBioPortal API, survival correlation)
```

---

## Quick Start

### Prerequisites

```bash
# Install conda/mamba
conda env create -f pipeline/envs/m6a_env.yaml
conda activate m6a_pipeline

# Or install manually
pip install m6anet xpore pandas numpy scipy requests
```

### Run the Full Pipeline

```bash
# 1. Prepare data (already done if using SG-NEx pre-processed files)
snakemake --snakefile pipeline/Snakefile --cores 8

# 2. Run cross-validation and annotation
python methods_comparison.py

# 3. Run TCGA validation
python tcga_query.py
```

### Using Pre-Processed SG-NEx Data

The SG-NEx project provides pre-processed eventalign files that skip basecalling and alignment:

```bash
# Download A549 data (replace with actual SG-NEx AWS paths)
aws s3 cp s3://sg-nex-data/data/processed_data/m6Anet/SGNex_A549_directRNA_replicate6_run1/ \
    data/A549_rep6/ --recursive --no-sign-request

# Download HepG2 data
aws s3 cp s3://sg-nex-data/data/processed_data/m6Anet/SGNex_HepG2_directRNA_replicate6_run1/ \
    data/HepG2_rep6/ --recursive --no-sign-request
```

---

## Results Summary

### Top Differentially Methylated Genes

| Gene | Direction | Differential Rate | Cancer Relevance |
|------|-----------|-------------------|-----------------|
| UBE3C | Higher in Lung | +0.9999 | Degrades tumour suppressor p53 |
| CD24 | Higher in Lung | +0.9999 | Immune evasion checkpoint |
| RPS6KA3 | Higher in Lung | +0.9999 | Accelerates protein synthesis |
| AHCY | Higher in Lung | +0.913 | Controls methylation fuel (SAM cycle) |
| CTCF | Higher in Lung | — | 3D genome organizer, cancer mutations |
| FADD | Cross-validated | — | Apoptosis signalling |
| METTL13 | Higher in Liver | -0.9999 | Methyltransferase enzyme |

### Top Enriched Pathways

| Pathway | Overlap | Adjusted p-value |
|---------|---------|-----------------|
| Protein processing in ER | 55/171 | 9.2 × 10⁻¹⁶ |
| Mitochondrial Gene Expression | 39/103 | 1.95 × 10⁻¹² |
| RNA transport | 48/186 | 6.6 × 10⁻¹⁰ |
| Ribosome Biogenesis | 42/155 | 1.7 × 10⁻⁸ |
| Ubiquitin-mediated Proteolysis | 66/319 | 2.1 × 10⁻⁸ |

---

## Methods

### m6Anet
- **Publication:** Hendra et al., *Nature Methods* 2022
- **Approach:** Multiple Instance Learning (MIL) neural network
- **Input:** Nanopore direct RNA signal features per DRACH motif position
- **Output:** Site modification probability + per-molecule stoichiometry

### xPore
- **Publication:** Pratanwanich et al., *Nature Biotechnology* 2021
- **Approach:** Gaussian Mixture Model (GMM) statistical framework
- **Input:** Aligned nanopore current signals from two conditions
- **Output:** Differential modification rate (−1 to +1) + adjusted p-value

### Cross-Validation
Both tools were run independently and genes identified by both were retained. The observed overlap (1,745 genes) significantly exceeds random expectation (~944) by hypergeometric test (p < 0.001).

---

## References

1. Hendra et al. (2022). Detection of m6A from direct RNA sequencing using a multiple instance learning framework. *Nature Methods*.
2. Pratanwanich et al. (2021). Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. *Nature Biotechnology*.
3. Maestri et al. (2024). Benchmark of m6A detection tools on direct RNA sequencing data. *Briefings in Bioinformatics*.
4. Chen et al. (2025). A systematic benchmark of Nanopore long-read RNA sequencing for transcript-level analysis in human cell lines. *Nature Methods*.
5. SG-NEx Consortium GitHub: [GoekeLab/sg-nex-data](https://github.com/GoekeLab/sg-nex-data)

---

## License

MIT License — see [LICENSE](LICENSE)

---

## Authors

**Team 12 — DNA Sequencing and Technologies**
- K. Siva Pranav (CB.SC.U4AIE23309)
- K. Mohana Venkata Krishna (CB.SC.U4AIE23313)
- P.V. Bhargav Reddy (CB.SC.U4AIE23329)
