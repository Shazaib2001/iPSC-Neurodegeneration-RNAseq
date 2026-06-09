# iPSC Neurological Disease Modelling — Transcriptomic Signature Discovery in ALS and Parkinson's Disease

## Overview

This project analyses RNA-seq gene expression data from iPSC-derived neurons modelling two neurodegenerative diseases — Amyotrophic Lateral Sclerosis (ALS) and Parkinson's Disease. The pipeline identifies disease-associated transcriptional signatures, uncovers disrupted biological pathways, and builds machine learning classifiers to validate the predictive power of those signatures.

The analysis is built entirely in R and covers differential expression analysis, pathway enrichment, cross-disease comparison, and random forest classification.

---

## Biological Background

### ALS (Amyotrophic Lateral Sclerosis)
ALS is a fatal neurodegenerative disease that destroys motor neurons — the cells that control voluntary movement. A defining molecular feature of ALS is the mislocalisation and aggregation of the RNA-binding protein TDP-43 (encoded by *TARDBP*), which disrupts RNA metabolism across hundreds of downstream genes. This dataset models ALS by knocking down TDP-43 in iPSC-derived spinal motor neurons using siRNA.

### Parkinson's Disease
Parkinson's is a progressive neurodegenerative disease caused by the death of dopaminergic neurons in the substantia nigra. The protein alpha-synuclein (encoded by *SNCA*) forms toxic aggregates known as Lewy bodies in affected neurons. This dataset models Parkinson's using iPSC-derived neurons carrying an SNCA triplication — a natural genetic mutation that overproduces alpha-synuclein and causes early-onset Parkinson's disease.

---

## Datasets

Both datasets were obtained from NCBI GEO (Gene Expression Omnibus).

| Dataset | GEO Accession | Disease | Model | Samples |
|---|---|---|---|---|
| ALS | GSE304550 | ALS | TDP-43 knockdown in iPSC-derived spinal motor neurons | 8 (4 SCR control, 4 TDP knockdown) |
| Parkinson's | GSE299913 | Parkinson's Disease | SNCA triplication vs isogenic corrected control | 6 (3 Triplication, 3 Isogenic) |

### Downloading the data

1. Go to [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/)
2. Search for each accession number (GSE304550, GSE299913)
3. Download the processed count files
4. Place both files in a `data/` folder in the project directory

---

## Project Structure

```
iPSC-Disease-Modelling/
├── iPSC_Neurological_Disease_Modelling.R   # Main analysis script
├── README.md                               # This file
└── data/
    ├── ALS_MN.samples.counts.processed_GSE304550.txt
    └── Parkinsons_RNA_seq_Norm_Counts_GSE299913.txt
```

---

## Pipeline

The analysis follows this structure for each disease, followed by a cross-disease comparison and machine learning step:

```
Data loading and cleaning
        ↓
Exploratory Data Analysis (EDA)
— Count distribution boxplot
— Library size check (ALS only)
— PCA for sample clustering and batch effect detection
        ↓
Differential Expression Analysis
— ALS: voom + limma (raw counts, batch-corrected design)
— Parkinson's: log2 transform + limma (pre-normalised data)
        ↓
Visualisation
— Volcano plot
— Heatmap of top 50 DEGs
        ↓
Pathway Enrichment
— GO biological process enrichment
— KEGG pathway enrichment
        ↓
Cross-disease Comparison
— Shared DEGs between ALS and Parkinson's
— Shared KEGG pathways
— Overlap visualisation
        ↓
Machine Learning
— Random forest classifiers (one per disease)
— Feature importance analysis
— Comparison of top predictive genes across diseases
```

---

## Key Results

### ALS
- **690 differentially expressed genes** identified (456 downregulated, 234 upregulated)
- More genes downregulated than upregulated — consistent with TDP-43's role in stabilising RNA transcripts
- **204 enriched GO pathways** (exploratory, raw p-values used due to small sample size)
- **3 significant KEGG pathways** — Cytoskeleton in muscle cells, Rap1 signalling, Proteoglycans in cancer
- Random forest classifier correctly predicted **2/2 test samples**

### Parkinson's
- **1187 differentially expressed genes** identified (695 upregulated, 492 downregulated)
- More genes upregulated — consistent with SNCA overexpression driving gene activation
- **275 enriched GO pathways** passing strict FDR correction — top hits include synaptic vesicle cycle, neurotransmitter transport, dopaminergic signalling
- **32 significant KEGG pathways** — top hit is dopaminergic synapse, directly relevant to Parkinson's pathology
- Random forest classifier correctly predicted **2/2 test samples**

### Cross-disease Comparison
- **29 shared DEGs** between ALS and Parkinson's
- **3 shared KEGG pathways** — Cytoskeleton in muscle cells, Proteoglycans in cancer, Rap1 signalling
- No overlap in top ML feature importance genes
- The two diseases are largely molecularly distinct but converge on shared cytoskeletal and signalling pathway disruption

---

## R Packages Required

### Bioconductor
```r
BiocManager::install(c(
  "limma",
  "edgeR",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ReactomePA",
  "AnnotationDbi"
))
```

### CRAN
```r
install.packages(c(
  "tidyverse",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "pROC",
  "caret",
  "glmnet",
  "randomForest",
  "e1071"
))
```

---

## How to Run

1. Clone the repository
2. Download both datasets from GEO and place in the `data/` folder
3. Open `iPSC_Neurological_Disease_Modelling.R` in RStudio
4. Update file paths at the top of the script if needed
5. Run the script top to bottom

All plots are generated inline. No additional setup required.

---

## Limitations

- **Small sample sizes** — ALS: n=4 per group, Parkinson's: n=3 per group. Statistical power is limited.
- **ALS GO enrichment** — strict FDR correction was too conservative given the small sample size. Raw p-values were used as an exploratory measure and results should be interpreted cautiously.
- **ML models** — with fewer than 10 samples per dataset, classification results are proof-of-concept rather than robust predictive models. Results would require validation on larger independent datasets.
- **ALS model type** — the ALS dataset uses a TDP-43 knockdown model rather than patient-derived cells. This simulates a key disease mechanism but does not capture the full genetic complexity of ALS.

---

## Author

Elvis  
MSc Artificial Intelligence in the Biosciences — Queen Mary University of London  
BSc Biological Sciences — Royal Holloway, University of London

---

## Data Sources

- Alessandrini F, Wright M, Kurosaki T et al. TDP-43 dysfunction compromises UPF1-dependent mRNA metabolism in ALS. *Neuron* 2026; 114(4):640-660. PMID: 41389796
- GSE299913 — SNCA triplication iPSC neuronal model. NCBI GEO.
