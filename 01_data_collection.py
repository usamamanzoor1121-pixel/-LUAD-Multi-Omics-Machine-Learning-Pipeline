"""
=============================================================================
STEP 1: DATA COLLECTION — TCGA-LUAD Multi-Omics Pipeline
=============================================================================
Cancer Type  : Lung Adenocarcinoma (TCGA-LUAD)
Data Sources : GDC Portal (https://portal.gdc.cancer.gov/)
Data Types   : RNA-Seq (HTSeq-FPKM-UQ), Somatic Mutations (MAF),
               DNA Methylation (Illumina 450k)
Author       : USAMA MANZOOR
=============================================================================

INSTRUCTIONS:
1. Install GDC client: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
2. Run this script to generate manifests & download data
3. Alternatively use TCGAbiolinks in R (see R section at bottom)
=============================================================================
"""

import os
import json

# ── Directory structure ──────────────────────────────────────────────────────
DIRS = [
    "data/raw/rnaseq",
    "data/raw/mutations",
    "data/raw/methylation",
    "data/processed",
    "data/results",
    "figures",
    "scripts",
]
for d in DIRS:
    os.makedirs(d, exist_ok=True)
print("✅ Directory structure created.")

# ── GDC API Query — RNA-Seq ──────────────────────────────────────────────────
rnaseq_params = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
        ]
    },
    "format": "TSV",
    "fields": "file_id,file_name,cases.case_id,cases.samples.sample_type",
    "size": "700"
}

# ── GDC API Query — Somatic Mutations ───────────────────────────────────────
mutation_params = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
        ]
    },
    "format": "TSV",
    "fields": "file_id,file_name,cases.case_id",
    "size": "10"
}

# ── GDC API Query — DNA Methylation ─────────────────────────────────────────
methylation_params = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["DNA Methylation"]}},
            {"op": "in", "content": {"field": "files.platform", "value": ["Illumina Human Methylation 450"]}},
        ]
    },
    "format": "TSV",
    "fields": "file_id,file_name,cases.case_id",
    "size": "700"
}

# Save manifests
with open("data/raw/rnaseq_manifest.json", "w") as f:
    json.dump(rnaseq_params, f, indent=2)
with open("data/raw/mutations_manifest.json", "w") as f:
    json.dump(mutation_params, f, indent=2)
with open("data/raw/methylation_manifest.json", "w") as f:
    json.dump(methylation_params, f, indent=2)

print("✅ GDC query manifests saved.")

# ── Download commands (run in terminal) ─────────────────────────────────────
download_commands = """
# ============================================================
# GDC DOWNLOAD COMMANDS — Run these in your terminal
# ============================================================

# 1. Download GDC Data Transfer Tool first:
#    https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

# 2. Submit queries to get file manifests:
curl --request POST --header "Content-Type: application/json" \\
     --data @data/raw/rnaseq_manifest.json \\
     "https://api.gdc.cancer.gov/files" \\
     -o data/raw/rnaseq_files.tsv

curl --request POST --header "Content-Type: application/json" \\
     --data @data/raw/mutations_manifest.json \\
     "https://api.gdc.cancer.gov/files" \\
     -o data/raw/mutation_files.tsv

curl --request POST --header "Content-Type: application/json" \\
     --data @data/raw/methylation_manifest.json \\
     "https://api.gdc.cancer.gov/files" \\
     -o data/raw/methylation_files.tsv

# 3. Download using GDC client:
gdc-client download -m data/raw/rnaseq_files.tsv    -d data/raw/rnaseq/
gdc-client download -m data/raw/mutation_files.tsv  -d data/raw/mutations/
gdc-client download -m data/raw/methylation_files.tsv -d data/raw/methylation/

# ============================================================
# ALTERNATIVE: Use TCGAbiolinks in R (RECOMMENDED — easier)
# ============================================================
# See 01b_TCGAbiolinks_download.R for the R script
"""

with open("data/download_commands.sh", "w") as f:
    f.write(download_commands)

# ── R script for TCGAbiolinks ────────────────────────────────────────────────
r_script = '''
# ============================================================
# TCGAbiolinks Download Script — TCGA-LUAD Multi-Omics
# ============================================================
# Install: BiocManager::install("TCGAbiolinks")
# ============================================================

library(TCGAbiolinks)
library(SummarizedExperiment)

# ── 1. RNA-Seq ─────────────────────────────────────────────
query_rna <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_rna, method = "api", files.per.chunk = 10)
rna_data <- GDCprepare(query_rna)
saveRDS(rna_data, "data/raw/TCGA_LUAD_RNAseq.rds")

# ── 2. Somatic Mutations ───────────────────────────────────
query_mut <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query_mut)
maf_data <- GDCprepare(query_mut)
write.csv(maf_data, "data/raw/TCGA_LUAD_mutations.csv", row.names = FALSE)

# ── 3. DNA Methylation ─────────────────────────────────────
query_meth <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_meth, method = "api", files.per.chunk = 5)
meth_data <- GDCprepare(query_meth)
saveRDS(meth_data, "data/raw/TCGA_LUAD_Methylation.rds")

cat("✅ All data downloaded successfully!\\n")
cat("Samples — RNA:", ncol(rna_data), "| Methylation:", ncol(meth_data), "\\n")
'''

with open("data/01b_TCGAbiolinks_download.R", "w") as f:
    f.write(r_script)

print("✅ R download script saved.")
print("\n📋 DATASET SUMMARY (Expected):")
print("   Project       : TCGA-LUAD (Lung Adenocarcinoma)")
print("   Tumor samples : ~513")
print("   Normal samples: ~59")
print("   Mutations     : ~230,000+ somatic variants")
print("   CpG probes    : ~485,000")
print("   Shared patients (all 3 types): ~300+")
print("\n✅ Step 1 complete. Run download commands or TCGAbiolinks R script.")
