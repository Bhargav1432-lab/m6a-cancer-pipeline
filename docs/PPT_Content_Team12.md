# Final Presentation Content (Slide-by-Slide)
## Project: Direct RNA Methylation Mapping in Cancer Transcriptomes (Team 12)
### Comparative Analysis of A549 (Lung) vs HepG2 (Liver) with TCGA Clinical Validation

---

## SLIDE 1: Title Slide
**Direct RNA Methylation Mapping in Cancer Transcriptomes**
*A Dual-Tool Computational Pipeline Using m6Anet and xPore with Clinical Survival Validation*

- **Course:** DNA Sequencing and Technologies
- **Team 12:** K. Siva Pranav (CB.SC.U4AIE23309), K. Mohana Venkata Krishna (CB.SC.U4AIE23313), P.V. Bhargav Reddy (CB.SC.U4AIE23329)
- **Institution:** Amrita Vishwa Vidyapeetham | October 2026

---

## SLIDE 2: Biological Background & Motivation
- **The Epitranscriptome:** Chemical modifications on RNA dynamically control mRNA stability, translation, and decay.
- **N⁶-methyladenosine (m6A):** The most prevalent and influential internal modification on eukaryotic mRNA.
- **Role in Oncology:** Cancer cells reprogram m6A "writers" (METTL3), "erasers" (ALKBH5), and "readers" (YTHDF1/2) to evade immune detection, accelerate growth, and drive drug resistance.
- **The Challenge:** How can we accurately map and quantify m6A changes between different cancers without expensive, destructive wet-lab assays?

---

## SLIDE 3: The Problem & Research Gaps
1. **The Common-Cancer Bias:** Most published m6A direct RNA maps focus on single cell lines. Systematic comparative analyses across lineages (e.g., Lung vs Liver) remain rare.
2. **The Stoichiometry Gap:** Legacy methods (MeRIP-seq) provide only binary ("methylated or not") signals. They cannot quantify the percentage of modified transcripts.
3. **The Tool Complementarity Gap:** Machine learning callers (m6Anet) and statistical mixture models (xPore) are rarely benchmarked together on the same cancer dataset.
4. **Our Goal:** Build a purely computational, reproducible pipeline using Oxford Nanopore direct RNA sequencing to map differential m6A and validate clinical significance in TCGA.

---

## SLIDE 4: Why Oxford Nanopore Direct RNA Sequencing?
| Feature | Traditional Methods (MeRIP-seq) | Direct RNA Sequencing (ONT) |
|:---|:---|:---|
| **Sample Handling** | Harsh chemical fragmentation & antibodies | Native, intact, unfragmented RNA |
| **Modification Detection** | Indirect (bulk enrichment peak) | Direct (ionic current disruption) |
| **Resolution** | Low (~100–200 nt window) | Single-base & single-molecule |
| **Stoichiometry** | ❌ Cannot measure | ✅ Quantitative (per-molecule) |
| **Wet-Lab Cost** | High experimental overhead | Purely computational on public data |

---

## SLIDE 5: Data Sources & Architecture
- **Source:** Singapore Nanopore Expression Project (SG-NEx; Genome Institute of Singapore).
- **Public Repository:** AWS Open Data Registry (`registry.opendata.aws/sgnex`).
- **Cell Lines Selected:**
  - **A549:** Human lung adenocarcinoma (ATCC CCL-185; Replicate 6, Run 1)
  - **HepG2:** Human hepatocellular carcinoma (ATCC HB-8065; Replicate 6)
- **Clinical Validation Cohorts:** The Cancer Genome Atlas (TCGA-LUAD, n=513; TCGA-LIHC, n=371).

---

## SLIDE 6: Dual-Tool Analytical Pipeline
```
[Raw Nanopore Signal: fast5/pod5]
            │
            ▼
[Dorado Basecalling & minimap2 Alignment]
            │
            ▼
[Signal Realignment: f5c / nanopolish eventalign]
            │
      ┌─────┴─────────────────────────┐
      ▼                               ▼
[m6Anet (Nature Methods 2022)]   [xPore (Nature Biotech 2021)]
• Multiple Instance Learning NN   • Gaussian Mixture Model (GMM)
• 187,022 positions evaluated     • 7,200,000+ positions tested
• Single-molecule stoichiometry   • 138,512 differential sites (padj<0.05)
      └─────┬─────────────────────────┘
            ▼
[Cross-Validation: Hypergeometric Test]
            ▼
[1,745 High-Confidence Shared Genes (1.85× Expected Overlap)]
            ▼
[Genomic Region Mapping (3'UTR) & Functional Enrichment (KEGG/GO)]
            ▼
[Clinical TCGA Kaplan-Meier Survival Validation]
```

---

## SLIDE 7: Methodological Comparison (m6Anet vs xPore)
| Evaluation Metric | m6Anet | xPore |
|:---|:---|:---|
| **Underlying Engine** | Deep Neural Network (MIL) | Gaussian Mixture Model (GMM) |
| **Operational Input** | Single condition (de novo) | Two-condition comparative |
| **Motif Constraints** | DRACH motif restricted | All transcriptomic positions |
| **Output Type** | Posterior probability & stoichiometry | Differential rate (−1.0 to +1.0) |
| **Multiple Testing** | Cutoff threshold (prob ≥ 0.9) | Benjamini–Hochberg FDR (padj < 0.05) |
| **Complementary Value**| Quantifies modified fraction | Identifies lineage-specific shifts |

---

## SLIDE 8: Result 1 — Robust Cross-Validation
- **Observed Overlap:** **1,745 genes** identified by BOTH independent tools.
- **Expected by Random Chance:** **~944 genes**.
- **Enrichment Ratio:** **1.85× over random expectation**.
- **Statistical Significance:** **p < 0.001** (Hypergeometric test).
- **Takeaway:** Cross-validation mathematically confirms that signals reflect authentic RNA biology rather than algorithmic noise.

---

## SLIDE 9: Result 2 — Asymmetric Lung Hypermethylation
- Total Cross-Validated Genes: **1,745**
- **Higher in Lung (A549):** **1,362 genes (78.1%)**
- **Higher in Liver (HepG2):** **383 genes (21.9%)**
- **Asymmetry Ratio:** **3.6-fold bias toward lung cancer hypermethylation**.
- **Ultra-Significant Targets:** 203 genes with padj = 0.0.
- **Biological Implication:** Lung adenocarcinoma exhibits an extensively hypermethylated epitranscriptome to sustain aggressive proliferation and translation.

---

## SLIDE 10: Result 3 — Genomic Topology Validation
- Distribution of top 100 differential sites across transcript regions:
  - **3'UTR:** **61%**
  - **CDS (Coding Sequence):** **34%**
  - **5'UTR:** **4%**
  - **Unmapped:** **1%**
- **Validation Milestone:** Established m6A biology shows characteristic clustering near stop codons and 3'UTRs. Our 61% 3'UTR distribution confirms the algorithm detected genuine m6A marks.

---

## SLIDE 11: Result 4 — Functional Pathway Enrichment
Differentially methylated genes strongly cluster into 4 core oncogenic hallmarks:
1. **Endoplasmic Reticulum & Proteostasis:**
   - *Protein processing in ER:* 55 genes (padj = 9.2 × 10⁻¹⁶)
   - *Ubiquitin-dependent catabolism:* 66 genes (padj = 2.1 × 10⁻⁸)
2. **Translation & Ribosomes:**
   - *Mitochondrial gene expression:* 39 genes (padj = 1.95 × 10⁻¹²)
   - *Ribosome biogenesis:* 42 genes (padj = 1.7 × 10⁻⁸)
   - *Translation initiation:* 55 genes (padj = 9.3 × 10⁻⁹)
3. **RNA Processing & Transport:**
   - *RNA transport:* 48 genes (padj = 6.6 × 10⁻¹⁰)
   - *Spliceosome assembly:* 47 genes (padj = 5.7 × 10⁻⁷)
4. **Proliferation:**
   - *Cell cycle regulation:* 34 genes (padj = 1.0 × 10⁻⁷)

---

## SLIDE 12: Result 5 — Clinical Survival Validation in TCGA
We linked dry-lab m6A candidates to clinical outcomes across **884 cancer patients**:
- **TCGA-LUAD (Lung Adenocarcinoma, n=513)**
- **TCGA-LIHC (Liver Hepatocellular Carcinoma, n=371)**

| Gene | m6A Shift in Lung | LUAD Hazard Ratio (p-value) | LIHC Hazard Ratio (p-value) | Clinical Outcome |
|:---|:---|:---|:---|:---|
| **CD24** | Hypermethylated (+0.9999) | **HR = 1.48 (p = 0.0076)** | HR = 1.34 (p = 0.098) | Adverse Prognosis in Lung |
| **UBE3C** | Hypermethylated (+0.9999) | **HR = 1.42 (p = 0.018)** | **HR = 1.55 (p = 0.014)** | Adverse Prognosis (Both) |
| **AHCY** | Hypermethylated (+0.9131) | **HR = 1.39 (p = 0.026)** | **HR = 1.62 (p = 0.0068)** | Adverse Prognosis (Both) |
| **RPS6KA3**| Hypermethylated (+0.9999) | **HR = 1.35 (p = 0.043)** | HR = 1.21 (p = 0.28) | Selective Lung Risk Factor |
| **FADD** | Cross-validated Target | **HR = 1.37 (p = 0.031)** | **HR = 1.45 (p = 0.039)** | Adverse Prognosis (Both) |

---

## SLIDE 13: Spotlight: CD24 (Immune Evasion)
- **Normal Function:** Cell surface glycoprotein engaging Siglec-10 on macrophages to deliver an inhibitory "don't eat me" signal.
- **In Lung Cancer:** High CD24 blocks phagocytosis, blunting anti-PD-1 immunotherapies.
- **Our Direct RNA Finding:** CD24 mRNA is near-completely m6A-methylated in A549 (+0.9999) compared to HepG2.
- **TCGA Validation:** Patients with high CD24 show significantly worse survival (**HR = 1.48, p = 0.0076**).
- **Conclusion:** m6A hypermethylation stabilizes CD24 mRNA, driving cancer immune evasion.

---

## SLIDE 14: Spotlight: AHCY (Self-Regulating Methylation Loop)
- **The Role of AHCY:** Hydrolyzes S-adenosylhomocysteine (SAH), which is a potent feedback inhibitor of methyltransferases.
- **The Finding:** AHCY itself is hypermethylated in lung cancer (+0.9131).
- **Novel Hypothesis:** Lung cancer uses m6A to regulate the enzyme that controls cellular methylation—forming a self-amplifying transmethylation feedback loop.
- **TCGA Validation:** Significant adverse survival factor in both LUAD (**HR = 1.39, p = 0.026**) and LIHC (**HR = 1.62, p = 0.0068**).

---

## SLIDE 15: Spotlight: UBE3C & RPS6KA3
- **UBE3C (E3 Ubiquitin Ligase):**
  - Targets p53/AHNAK tumor-suppressor complexes for proteasomal degradation.
  - Hypermethylated in A549 (+0.9999); adverse survival in TCGA (**HR = 1.42, p = 0.018**).
- **RPS6KA3 (RSK2 Kinase):**
  - Acts downstream of Ras-MAPK to phosphorylate RPS6 and eIF4B, boosting protein translation.
  - Hypermethylated in A549 (+0.9999); adverse survival in TCGA (**HR = 1.35, p = 0.043**).

---

## SLIDE 16: Pipeline Engineering & Reproducibility
- **End-to-End Orchestration:** Built with **Snakemake** (`pipeline/Snakefile`).
- **Configuration Management:** Centralized YAML configuration (`pipeline/config.yaml`).
- **Dependency Isolation:** Conda environment specification (`pipeline/envs/m6a_env.yaml`).
- **Open Science Artifacts:** Fully structured GitHub-ready repository with modular scripts, notebooks, and result tables (`m6a_cancer_pipeline/`).

---

## SLIDE 17: Limitations & Methodological Reflections
- **In Vitro Models:** A549 and HepG2 are established cell line models; clinical in vivo microenvironments exhibit additional multicellular complexity.
- **Experimental Validation:** All findings are computational; physical assays (e.g., SELECT or m6A-IP-qPCR) are ideal future validation steps.
- **Replicate Scope:** Analyzed high-depth single biological replicates; future runs will scale across multi-replicate series.
- **Exploratory Horizons:** Pipeline is fully architected to ingest multiple myeloma patient samples from SG-NEx as an immediate extension.

---

## SLIDE 18: Summary & Final Conclusions
1. **Pipeline Validated:** Successfully engineered a dual-tool workflow combining m6Anet and xPore for direct RNA sequencing.
2. **Biological Authenticity:** 1,745 cross-validated genes (1.85× expected by chance, p < 0.001) with 61% 3'UTR enrichment.
3. **Phenotypic Discovery:** Lung adenocarcinoma exhibits an asymmetric 3.6-fold hypermethylation program focused on proteostasis, ribosome biogenesis, and immune escape.
4. **Clinical Translation:** Top candidates (CD24, UBE3C, AHCY, RPS6KA3) significantly stratify patient survival in TCGA cohorts.

---

## SLIDE 19: References
1. Hendra C, et al. (2022). m6Anet. *Nature Methods*.
2. Pratanwanich PN, et al. (2021). xPore. *Nature Biotechnology*.
3. Chen Y, et al. (2025). SG-NEx benchmark. *Nature Methods*.
4. Maestri S, et al. (2024). m6A tool benchmark. *Briefings in Bioinformatics*.
5. Meyer KD, et al. (2012). Topology of mRNA methylation. *Cell*.
6. Müller S, et al. (2019). CD24 immune checkpoint. *Cancer Research*.
7. He L, et al. (2021). UBE3C in lung cancer. *Oncogene*.
8. Cerami E, et al. (2012). cBioPortal. *Cancer Discovery*.

---
**Thank You — Questions & Discussion**
