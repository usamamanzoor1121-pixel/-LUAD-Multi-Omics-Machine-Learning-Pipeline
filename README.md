# 🫁 LUAD Multi-Omics & Machine Learning Pipeline

> **Lung Adenocarcinoma · TCGA-LUAD · RNA-seq + WES + DNA Methylation + ML/DL**

[![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://r-project.org)
[![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)

---

## 📌 Overview

A comprehensive multi-omics study of **TCGA-LUAD** (Lung Adenocarcinoma) integrating:
- **RNA-sequencing**: 513 tumour / 59 normal samples
- **Somatic mutations**: Whole-exome sequencing
- **DNA methylation**: Illumina 450k array

The pipeline identifies master regulatory genes, convergent oncogenic pathways, and clinically actionable molecular targets — extended with a 5-model ML/DL framework for phenotype prediction.

---

## 🔬 Methods

### Multi-Omics Integration
| Layer | Data | Method | Threshold |
|---|---|---|---|
| Transcriptomics | RNA-seq (513T/59N) | DESeq2-equivalent framework | padj<0.05, \|log2FC\|>1.5 |
| Somatic Mutations | Whole-exome sequencing | MAF parsing | — |
| Epigenomics | Illumina 450k methylation | DMP calling | \|dBeta\|≥0.20 |
| PPI Network | STRING database | CytoHubba MCC | Top hub ranking |

### Machine Learning Extension
- **Feature matrix**: 43 features = top 38 DEG expression values + binary mutation status for TP53, KRAS, STK11, EGFR, KEAP1
- **Models**: Random Forest, SVM, Gradient Boosting Machine (GBM), MLP (Deep Learning), LASSO
- **Evaluation**: 5-fold stratified cross-validation + held-out test set (80:20 split; n=457 train / n=115 test)

---

## 📊 Key Results

### Differential Expression
- **12,500** genes expressed; **50 upregulated**, **30 downregulated** (padj<0.05, |log2FC|>1.5)
- Top upregulated: proliferation and cell cycle drivers
- Top downregulated: tumour suppressor and surfactant genes

### Somatic Mutation Landscape
| Gene | Mutation Frequency |
|---|---|
| TP53 | 47.8% |
| KRAS | 32.5% |
| EGFR | 14.1% |
| STK11 | — |
| KEAP1 | — |

### Epigenetic Analysis
- **23,400 hypermethylated CpGs** identified
- Key silenced tumour suppressors: **RASSF1A**, **CDKN2A**, **MLH1**

### PPI Network (STRING + CytoHubba MCC)
- **TP53** = Top hub gene (degree = 318, MCC Rank #1)
- All hub genes converge on **Cell Cycle Regulation** and **RTK/RAS Signalling**

---

## 🛠️ Tools & Environment

```bash
# Core bioinformatics
DESeq2 (R), GATK, Samtools, Picard, FastQC

# Network analysis
STRING PPI, Cytoscape, CytoHubba (MCC method)

# Machine learning
Python: Scikit-Learn, TensorFlow/Keras, Pandas, NumPy, Matplotlib, Seaborn

# Databases
TCGA-LUAD, GDC Data Portal, Illumina 450k, MSigDB
```

---

## 📁 Repository Structure

```
luad-multiomics-pipeline/
├── data/
│   ├── raw/                   # TCGA-LUAD expression, MAF, methylation
│   ├── processed/             # Normalized matrices, filtered DEGs
│   └── external/              # Reference annotations (GTF, CpG islands)
├── notebooks/
│   ├── 01_differential_expression.ipynb
│   ├── 02_somatic_mutation_analysis.ipynb
│   ├── 03_dna_methylation.ipynb
│   ├── 04_ppi_network.ipynb
│   └── 05_ml_dl_models.ipynb
├── scripts/
│   ├── deseq2_analysis.R
│   ├── methylation_analysis.R
│   ├── ppi_network.py
│   └── ml_pipeline.py
├── results/
│   ├── figures/
│   ├── tables/
│   └── models/
├── envs/
│   ├── environment_r.yml
│   └── environment_python.yml
└── README.md
```

---

## 🚀 Getting Started

```bash
git clone https://github.com/usamamanzoor1121-pixel/luad-multiomics-pipeline.git
cd luad-multiomics-pipeline

# Set up environments
conda env create -f envs/environment_r.yml
conda env create -f envs/environment_python.yml

# Download TCGA-LUAD data via GDC Data Portal
# Place in data/raw/

# Run differential expression
conda activate luad_r
Rscript scripts/deseq2_analysis.R

# Run ML pipeline
conda activate luad_python
python scripts/ml_pipeline.py
```

---

## 📚 Data Sources

- [TCGA-LUAD — GDC Data Portal](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)
- [UCSC Xena — TCGA expression + methylation](https://xenabrowser.net)
- [STRING PPI Database](https://string-db.org)
- [MSigDB Gene Sets](https://www.gsea-msigdb.org/gsea/msigdb)

---

## 👤 Author

**Usama Manzoor** · Bioinformatician, JSMU Diagnostic Laboratory, Karachi  
📧 usama.manzoor1121@gmail.com · 🐙 [@usamamanzoor1121-pixel](https://github.com/usamamanzoor1121-pixel)
