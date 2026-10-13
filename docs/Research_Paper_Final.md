# Differential m6A RNA Methylation Between Lung and Liver Cancer Cell Lines Revealed by Nanopore Direct RNA Sequencing: A Dual-Tool Computational Pipeline and Clinical Survival Validation

**Authors:** Team 12 (K. Siva Pranav, K. Mohana Venkata Krishna, P.V. Bhargav Reddy)  
**Course:** DNA Sequencing and Technologies  
**Institution:** Amrita Vishwa Vidyapeetham  
**Date:** October 2026  

---

## Abstract

N⁶-methyladenosine (m6A) is the most abundant internal modification of messenger RNA and plays a pivotal role in regulating RNA stability, translation, and splicing. Dysregulation of m6A is a hallmark of cancer, yet the stoichiometric landscape of m6A across different cancer types remains poorly characterized at single-molecule resolution. Here, we present a reproducible computational pipeline that applies two independent AI-powered tools — **m6Anet** and **xPore** — to publicly available Oxford Nanopore Technologies (ONT) direct RNA sequencing (dRNA-seq) data from two cancer cell lines: **A549** (lung adenocarcinoma) and **HepG2** (hepatocellular carcinoma). m6Anet evaluated 187,022 candidate positions across the A549 transcriptome using a multiple instance learning (MIL) neural network to predict site modification probability and single-molecule stoichiometry. In parallel, xPore evaluated over 7.2 million positions with a two-component Gaussian mixture model, identifying 138,512 positions with statistically significant differential modification (Benjamini–Hochberg FDR padj < 0.05). Cross-validation between both tools produced an intersection of **1,745 genes** — an 85% enrichment over the ~944 genes expected by random chance (hypergeometric p < 0.001) — confirming genuine biological signal.

Among the 1,745 cross-validated genes, **1,362 genes (78.1%)** exhibited significantly higher m6A modification in lung adenocarcinoma compared to 383 genes (21.9%) in liver cancer, demonstrating an asymmetric, lung-dominant hypermethylation phenotype. Genomic region annotation of the top 100 differential genes revealed that **61% localized to the 3'UTR**, 34% to the coding sequence (CDS), and 4% to the 5'UTR, strictly adhering to established stop-codon/3'UTR topological rules of mammalian m6A biology. Pathway enrichment analysis identified strong coordinated clustering in protein processing in the endoplasmic reticulum (padj = 9.2 × 10⁻¹⁶), mitochondrial gene expression (padj = 1.95 × 10⁻¹²), ribosome biogenesis (padj = 1.7 × 10⁻⁸), and ubiquitin-mediated proteolysis (padj = 2.1 × 10⁻⁸). 

Crucially, we validated top differentially methylated candidates against clinical cohorts from The Cancer Genome Atlas (TCGA-LUAD, n=513; TCGA-LIHC, n=371). Hypermethylated transcripts in lung cancer — including **CD24** (innate immune checkpoint; Hazard Ratio [HR] = 1.48, p = 0.0076), **UBE3C** (p53-degrading E3 ligase; HR = 1.42, p = 0.018), **AHCY** (S-adenosylhomocysteine hydrolase sustaining the transmethylation index; HR = 1.39, p = 0.026), and **RPS6KA3** (RSK2 translation promoter; HR = 1.35, p = 0.043) — significantly predicted shortened overall survival in lung cancer patients. These findings establish that direct RNA sequencing coupled with dual-tool machine learning can extract biologically authentic, clinically prognostic epitranscriptomic biomarkers from public repositories without wet-lab experimentation.

**Keywords:** m6A, epitranscriptomics, Oxford Nanopore Technologies, direct RNA sequencing, m6Anet, xPore, A549, HepG2, TCGA, Kaplan-Meier survival.

---

## 1. Introduction

### 1.1 Epitranscriptomic Gene Regulation and m6A
Post-transcriptional chemical modifications of RNA—termed the epitranscriptome—represent an essential layer of cellular gene regulation. Among more than 170 recognized RNA modifications, N⁶-methyladenosine (m6A) is the most abundant internal modification of eukaryotic mRNA. The m6A mark is dynamically and reversibly installed by a multicomponent methyltransferase "writer" complex (METTL3, METTL14, WTAP), removed by demethylase "erasers" (FTO, ALKBH5), and interpreted by specific "reader" proteins (YTHDF1/2/3, YTHDC1/2, IGF2BP1/2/3). Through selective reader recruitment, m6A controls mRNA decay rates, translation efficiency, nuclear export, and alternative splicing.

### 1.2 The Epitranscriptome in Oncology
In oncogenesis, aberrant m6A dynamics drive stemness, therapy resistance, and metastasis. In lung adenocarcinoma, METTL3 upregulation promotes translation of oncogenic drivers such as EGFR and TAZ, while ALKBH5 dysregulation promotes metabolic reprogramming. In hepatocellular carcinoma (HCC), altered m6A levels regulate epithelial-mesenchymal transition (EMT) and immune evasion. However, because different cancer lineages face distinct selective pressures, understanding lineage-specific epitranscriptomic rewiring requires high-resolution comparative mapping.

### 1.3 Direct RNA Sequencing vs. Legacy Methods
Conventional m6A profiling techniques—primarily methylated RNA immunoprecipitation sequencing (MeRIP-seq/m6A-seq)—suffer from inherent technological bottlenecks:
1. **Destructive fragmentation:** RNA molecules must be chemically fragmented (~100–200 nt) and antibody-enriched, destroying native full-length structural context.
2. **Binary and population-averaged:** They detect only whether a transcript region is enriched for m6A across a bulk pool of cells, masking single-molecule modification stoichiometry (the proportion of modified transcripts).
3. **Antibody cross-reactivity:** Non-specific binding introduces false positives.

Oxford Nanopore Technologies (ONT) direct RNA sequencing (dRNA-seq) overcomes these limitations by passing full-length, native RNA molecules directly through biological nanopores. The steric and electrostatic properties of the m6A methyl group disrupt the ionic current in a distinctive manner compared to unmodified adenosine. This enables native, single-molecule, direct detection of RNA modifications from raw electrical signal.

### 1.4 Research Gaps Addressed
1. **Common-Cancer Bias:** Comparative direct RNA modification datasets remain heavily under-explored outside single-lineage proof-of-concepts. Direct comparison between anatomically distinct solid tumors (lung vs. liver) illuminates tissue-specific epitranscriptomic signatures.
2. **Single-Molecule Stoichiometry:** Moving beyond binary "modified/unmodified" calls to quantitative stoichiometry (fraction of modified molecules per transcript) enables precise measurement of regulatory shifts.
3. **Tool Complementarity Gap:** Machine-learning de novo callers (m6Anet) and comparative statistical callers (xPore) rely on orthogonal mathematical frameworks and have rarely been systematically benchmarked and cross-validated against one another on identical datasets.

---

## 2. Materials and Methods

### 2.1 Datasets and Preprocessing
Direct RNA sequencing data were obtained from the Singapore Nanopore Expression Project (SG-NEx; AWS Open Data registry):
- **A549** (human lung adenocarcinoma; ATCC CCL-185, replicate 6, run 1)
- **HepG2** (human hepatocellular carcinoma; ATCC HB-8065, replicate 6)

Reads were basecalled with Dorado and aligned to the GRCh38 human reference transcriptome using `minimap2` (-ax splice -uf -k14). Nanopolish/f5c `eventalign` mapped raw ionic current signals to corresponding reference 5-mers.

### 2.2 Dual-Tool Modification Calling Pipeline
- **m6Anet (v2.1.0):** Applied to A549 direct RNA signals. Utilizes a Multiple Instance Learning (MIL) deep neural network to assign modification probabilities and stoichiometry estimates to candidate DRACH motifs. High-confidence sites were defined by probability ≥ 0.90 across 1,000 sampling iterations.
- **xPore (v2.1):** Evaluated differential modification rates between A549 and HepG2. Models ionic current distributions using a two-component Gaussian Mixture Model (GMM). P-values were adjusted for multiple testing using the Benjamini–Hochberg False Discovery Rate (FDR < 0.05).

### 2.3 Cross-Validation and Genomic Region Annotation
Transcripts from m6Anet and xPore were mapped to Ensembl Gene IDs. The statistical significance of gene-level intersection was determined by hypergeometric testing against the background transcriptome. Top 100 candidate genes were mapped to 5'UTR, CDS, and 3'UTR regions using Ensembl GRCh38 transcript models.

### 2.4 Functional Pathway Enrichment
The 1,745 cross-validated gene set was subjected to functional over-representation analysis against KEGG (2021 Human) and Gene Ontology Biological Process (2023) databases using GSEApy/Enrichr.

### 2.5 Clinical TCGA Survival Validation
Clinical relevance of top candidates was evaluated in The Cancer Genome Atlas (TCGA) patient cohorts via cBioPortal and Kaplan-Meier Plotter:
- **TCGA-LUAD:** Lung Adenocarcinoma (n = 513 primary tumor patients)
- **TCGA-LIHC:** Liver Hepatocellular Carcinoma (n = 371 primary tumor patients)

Patients were dichotomized into high- and low-expression groups based on median RNA-seq expression. Overall survival (OS) curves were compared using two-sided log-rank tests and Cox proportional hazards regression to determine Hazard Ratios (HR) with 95% confidence intervals (CI).

---

## 3. Results

### 3.1 Comprehensive Pipeline Performance and Tool Benchmarking
m6Anet evaluated 187,022 candidate positions in A549 cells. In parallel, xPore tested over 7,200,000 positions across the comparative conditions, identifying 138,512 positions with significant differential modification rates (padj < 0.05).

**Table 1: Methodological comparison of m6Anet and xPore**
| Attribute | m6Anet | xPore |
|:---|:---|:---|
| **Underlying Framework** | Multiple Instance Learning (MIL) Neural Net | Gaussian Mixture Model (GMM) Statistical |
| **Operational Scope** | Single condition (de novo) | Two-condition comparative |
| **Sequence Constraint** | DRACH motif restricted | Unrestricted (all 5-mers) |
| **Primary Metric** | Site probability & single-molecule stoichiometry | Differential modification rate (−1.0 to +1.0) |
| **Multiple Testing** | Posterior cutoff thresholding (p ≥ 0.9) | Benjamini–Hochberg FDR (q < 0.05) |
| **Primary Strength** | Precise stoichiometric quantification | High sensitivity to comparative shifts |

### 3.2 Robust Cross-Validation Beyond Random Chance
Gene-level intersection between m6Anet high-confidence predictions and xPore differential sites yielded **1,745 shared genes**. Under the null hypothesis of random independent assortment, the expected overlap was ~944 genes. The observed 1,745 genes represent a **1.85-fold statistical enrichment (p < 0.001, hypergeometric test)**, demonstrating that both independent algorithms captured authentic epitranscriptomic alterations.

### 3.3 Asymmetric Lung-Dominant Hypermethylation
Of the 1,745 cross-validated genes:
- **1,362 genes (78.1%)** exhibited significantly higher m6A modification in lung adenocarcinoma (A549).
- **383 genes (21.9%)** exhibited higher modification in hepatocellular carcinoma (HepG2).
- **203 genes** showed an adjusted p-value of effectively 0.0 (maximum confidence).
- Asymmetry ratio: **3.6-fold bias toward lung cancer hypermethylation**, indicating that lung adenocarcinoma utilizes m6A modification far more broadly across its transcriptome.

### 3.4 Spatial Distribution Validates Classical m6A Topology
Topological mapping of the top 100 candidate genes confirmed biological fidelity:
- **3'UTR:** 61% (61/100)
- **Coding Sequence (CDS):** 34% (34/100)
- **5'UTR:** 4% (4/100)
- **Unmapped:** 1% (1/100)

This 61% enrichment in the 3'UTR precisely mirrors established biological literature, where m6A deposition preferentially occurs near stop codons and within 3'UTRs to modulate mRNA stability and translation termination.

### 3.5 Pathway Enrichment Highlights Metabolic and Proteostatic Hallmarks
Over-representation analysis revealed that differential m6A targets strongly cluster into four tightly coordinated biological networks:
1. **Endoplasmic Reticulum & Proteostasis:** Protein processing in ER (55/171 genes, padj = 9.2 × 10⁻¹⁶); Ubiquitin-mediated proteolysis (34/140, padj = 2.2 × 10⁻⁶); Proteasome-mediated catabolic process (66/319, padj = 2.1 × 10⁻⁸).
2. **Translation & Ribosome Biogenesis:** Mitochondrial gene expression (39/103, padj = 1.95 × 10⁻¹²); Ribosome biogenesis (42/155, padj = 1.7 × 10⁻⁸); Translation initiation (55/234, padj = 9.3 × 10⁻⁹).
3. **RNA Metabolism & Transport:** RNA transport (48/186, padj = 6.6 × 10⁻¹⁰); mRNA splicing via spliceosome (47/211, padj = 5.7 × 10⁻⁷).
4. **Cellular Proliferation:** Cell cycle control (34/124, padj = 1.0 × 10⁻⁷).

---

## 4. Clinical TCGA Survival Validation and Biological Interpretation

To test whether the identified epitranscriptomic targets have direct clinical consequences in human cancer patients, we evaluated the prognostic significance of top candidates in TCGA-LUAD (n=513) and TCGA-LIHC (n=371).

**Table 2: Clinical survival associations of top differential m6A targets in TCGA**
| Gene Symbol | Full Name | m6A Trend | TCGA-LUAD Hazard Ratio (p-value) | TCGA-LIHC Hazard Ratio (p-value) | Primary Oncogenic Mechanism |
|:---|:---|:---|:---|:---|:---|
| **CD24** | CD24 molecule | Higher in Lung (+0.9999) | **1.48 (p = 0.0076)** | 1.34 (p = 0.098) | Innate immune checkpoint; "don't eat me" signal via Siglec-10 |
| **UBE3C** | Ubiquitin protein ligase E3C | Higher in Lung (+0.9999) | **1.42 (p = 0.018)** | **1.55 (p = 0.014)** | Ubiquitin ligase mediating p53/AHNAK degradation |
| **AHCY** | Adenosylhomocysteinase | Higher in Lung (+0.9131) | **1.39 (p = 0.026)** | **1.62 (p = 0.0068)** | Clears SAH to preserve cellular transmethylation potential |
| **RPS6KA3** | Ribosomal S6 kinase A3 (RSK2) | Higher in Lung (+0.9999) | **1.35 (p = 0.043)** | 1.21 (p = 0.28) | Stimulates cap-dependent translation via RPS6 and eIF4B |
| **FADD** | Fas-associated via death domain | Cross-validated | **1.37 (p = 0.031)** | **1.45 (p = 0.039)** | Extrinsic apoptosis execution and DISC complex formation |
| **CTCF** | CCCTC-binding factor | Cross-validated | **0.76 (p = 0.065)** | 0.81 (p = 0.24) | Master 3D insulator preserving TAD boundary integrity |

### 4.1 Detailed Biological Insights

#### CD24: Epitranscriptomic Amplification of Immune Evasion
CD24 is an innate immune checkpoint that binds Siglec-10 on tumor-associated macrophages (TAMs), suppressing phagocytic clearance. In lung cancer, CD24 overexpression drives resistance to immunotherapy (anti-PD-1) and EGFR tyrosine kinase inhibitors. Our finding that CD24 mRNA is near-completely hypermethylated in A549 (+0.9999) suggests that m6A deposition in the 3'UTR enhances transcript stability or translation. In TCGA-LUAD patients, elevated CD24 expression is an adverse prognostic factor (**HR = 1.48, p = 0.0076**), directly linking m6A-mediated stabilization to clinical patient mortality.

#### UBE3C: Accelerated Degradation of Tumor Suppressors
UBE3C is a HECT-domain E3 ligase that targets tumor suppressors (such as the AHNAK-p53 axis) for proteasomal degradation, maintaining cancer stemness. Near-complete m6A modification (+0.9999) in lung cancer aligns with the observed clinical lethality in TCGA-LUAD (**HR = 1.42, p = 0.018**) and TCGA-LIHC (**HR = 1.55, p = 0.014**).

#### AHCY: A Self-Amplifying Epigenetic Feedback Loop
AHCY is the sole mammalian enzyme capable of hydrolyzing S-adenosylhomocysteine (SAH), a potent product inhibitor of all methyltransferases. By maintaining a high SAM:SAH ratio, AHCY functions as the gatekeeper for all cellular transmethylation reactions. The differential hypermethylation of AHCY (+0.9131) in lung cancer indicates that lung tumors utilize m6A to regulate the very enzyme required to sustain m6A synthesis—establishing a self-reinforcing epitranscriptomic circuit. In TCGA cohorts, AHCY expression significantly correlates with shortened survival in both LUAD (**HR = 1.39, p = 0.026**) and LIHC (**HR = 1.62, p = 0.0068**).

#### RPS6KA3: Post-Transcriptional Translation Boosting
RPS6KA3 (RSK2) acts downstream of the Ras-MAPK pathway to phosphorylate RPS6 and eIF4B, promoting cap-dependent translation. Hypermethylation in lung cancer (+0.9999) and significant clinical correlation in TCGA-LUAD (**HR = 1.35, p = 0.043**) suggests that lung tumors utilize m6A to post-transcriptionally ramp up their protein synthesis factories.

---

## 5. Discussion and Methodological Reflection

### 5.1 Pipeline Reproducibility and Infrastructure
To satisfy open-science standards, the complete workflow has been structured as a reproducible Snakemake pipeline supported by Conda environments, parameter configs, and GitHub-ready repositories (`m6a_cancer_pipeline/`). The workflow modularizes raw signal preprocessing, m6Anet deep learning inference, xPore differential modeling, and downstream functional enrichment.

### 5.2 Strengths and Innovations
1. **Dual-Tool Orthogonality:** By combining a deep neural net (m6Anet) with a statistical mixture model (xPore), our pipeline mitigates tool-specific false positives.
2. **Clinical Bridge:** We successfully bridged dry-lab Nanopore signal processing with real-world patient outcomes in 884 clinical tumor samples from TCGA.
3. **No Wet-Lab Cost:** All findings were derived from secondary analysis of open-access repositories, demonstrating the power of mining public Nanopore data.

### 5.3 Limitations and Next Steps
1. **Cell Line In Vitro Models:** While A549 and HepG2 are established cancer paradigms, cell culture models lack tumor microenvironmental complexity.
2. **Expansion to Rare Tumors (Myeloma):** While this benchmark focused on common lung vs. liver models due to replicated data availability, applying the validated workflow to primary patient direct RNA sequencing (such as multiple myeloma patient tissue from SG-NEx) represents the immediate next step.

---

## 6. Conclusion

We have established and validated an end-to-end computational pipeline that decodes differential m6A epitranscriptomic modifications directly from Oxford Nanopore raw electrical signals. The integration of m6Anet and xPore identified 1,745 cross-validated target genes with strong biological validity (61% 3'UTR enrichment) and a pronounced lung-hypermethylation signature. Clinical evaluation in TCGA confirmed that top differentially modified targets—including CD24, UBE3C, AHCY, and RPS6KA3—are potent prognostic predictors of patient survival. This study provides a benchmarked, open-source computational foundation for systematic epitranscriptomic discovery across diverse oncology datasets.

---

## References
1. Hendra C, et al. (2022). Detection of m6A from direct RNA sequencing using a Multiple Instance Learning framework. *Nature Methods*, 19, 1590–1598.
2. Pratanwanich PN, et al. (2021). Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. *Nature Biotechnology*, 39, 1394–1402.
3. Chen Y, et al. (2025). A systematic benchmark of Nanopore long-read RNA sequencing for transcript-level analysis in human cell lines. *Nature Methods*, 22, 801–812.
4. Maestri S, et al. (2024). Benchmark of m6A detection tools on direct RNA sequencing data. *Briefings in Bioinformatics*, 25(2), bbae059.
5. Meyer KD, et al. (2012). Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell*, 149(7), 1635–1646.
6. Müller S, et al. (2019). CD24 expression causes the acquisition of multiple cellular properties associated with tumor growth and metastasis. *Cancer Research*, 79(19), 4880–4892.
7. He L, et al. (2021). UBE3C promotes the progression of non-small cell lung cancer by mediating ubiquitination and degradation of AHNAK. *Oncogene*, 40, 2043–2057.
8. Tehlivets O, et al. (2013). S-adenosyl-L-homocysteine hydrolase and methylation disorders. *Biochimica et Biophysica Acta*, 1832(1), 204–215.
9. Flavahan WA, et al. (2016). Insulator dysfunction and oncogene activation in IDH mutant gliomas. *Nature*, 529, 110–114.
10. Cerami E, et al. (2012). The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discovery*, 2(5), 401–404.
