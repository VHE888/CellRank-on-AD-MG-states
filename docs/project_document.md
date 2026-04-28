# Project Documentation

This directory contains the full documentation, code, and interpretation notes for the project:

## "Understanding Microglial State Transitions in Alzheimer's Disease Using CellRank Trajectory Analysis"

Our analysis aims to integrate RNA velocity and single-cell transcriptomic data to uncover transcriptional and regulatory features associated with microglial state changes during Alzheimer's disease progression. Using CellRank, we model the dynamic transitions from homeostatic microglia to disease-associated states (DAM, proliferative, inflammatory) in the context of AD pathology.

---

## Contents Overview

This directory contains a single comprehensive script that implements the complete CellRank analysis pipeline for microglial state transitions in Alzheimer's disease.

**Script:** `pipeline_cellrank.py`

The script is organized into sequential modules that cover the entire analytical workflow from data loading to final interpretation.

---

## Pipeline Modules

### 1. Setup and Imports
Loads all required libraries: CellRank, Scanpy, scVelo, and visualization tools. Sets global parameters and plotting themes.

### 2. Data Loading and Quality Control
Loads pre-processed AnnData object containing microglial single-cell data. Performs basic quality checks and filtering.

### 3. RNA Velocity Computation
Computes RNA velocity using scVelo dynamical model to infer directional dynamics of microglial state transitions. Parameters: mode="dynamical", n_top_genes=2000, n_pcs=30.

### 4. CellRank Kernel Construction
Combines velocity kernel (50%) and connectivity kernel (50%) to build transition matrix defining cell-to-cell movement probabilities.

### 5. Macrostate Identification (GPCCA)
Applies GPCCA to identify discrete microglial macrostates representing distinct biological conditions. Predicts terminal states (sinks) and initial states (sources).

### 6. Fate Probability Analysis
Calculates absorption probabilities for each cell toward each terminal state. Identifies cells committed to specific fates (DAM, proliferative, inflammatory).

### 7. Gene Expression Trends
Extracts genes that change expression along each microglial trajectory. Ranks genes by differential expression magnitude along pseudotime.

### 8. Regulatory Network Analysis (Optional)
Identifies upstream transcription factors driving microglial state transitions. Requires metabolic labeling data.

### 9. Results Saving and Summary
Saves updated AnnData object with all CellRank results. Generates summary statistics and output file locations.

---

## AD-Specific Interpretations

### Microglial States Identified

| State | Markers | Role in AD |
|-------|---------|------------|
| **Homeostatic** | *P2ry12*, *Tmem119*, *Cx3cr1* | Initial/resting state |
| **Early DAM** | *Apoe*, *Trem2*, *Tyrobp* | Transitional disease-associated |
| **Late DAM** | *Gpnmb*, *Spp1*, *Lpl* | Terminal disease-associated |
| **Proliferative** | *Mki67*, *Top2a*, *Cdk1* | Cell cycle active |
| **Inflammatory** | *Il1b*, *Tnf*, *Nos2* | Pro-inflammatory |


### Key Transcription Factors
- *MITF* - DAM regulation
- *KLF12* - Transcriptional repressor
- *SPI1/PU.1* - Microglial master regulator
- *CEBPB* - Inflammatory response

---

## Required Input Data

| File | Format | Description |
|------|--------|-------------|
| `microglia_adata.h5ad` | .h5ad | Pre-processed AnnData with spliced/unspliced counts, UMAP, and cell type annotations |

**Required `.obs` columns:**
- `cell_type` - Cluster/cell type annotations
- `batch` - Batch or donor information (optional)
- `metabolic_label` - For regulatory analysis (optional)
- `experimental_time` - For regulatory analysis (optional)

**Required `.layers`:**
- `spliced` - Spliced count matrix
- `unspliced` - Unspliced count matrix

**Required `.obsm`:**
- `X_umap` - UMAP embedding coordinates

---

## Outputs Generated

### Figures (saved to `./figures/`)

| File | Description |
|------|-------------|
| `velocity_stream.png` | RNA velocity streamlines on UMAP |
| `macrostates.png` | UMAP with colored macrostate regions |
| `fate_probabilities.png` | Fate probability plots for each terminal state |
| `gene_trends.png` | Gene expression dynamics along trajectories |
| `dam_regulators.csv` | TF activity scores (if regulatory analysis run) |

### Data Outputs

| File | Description |
|------|-------------|
| `microglia_adata_cellrank.h5ad` | Updated AnnData with CellRank results |

---

## Reproducibility

- The script is self-contained and can be executed from start to finish.
- All random seeds are fixed (`random_state = 42` where applicable).
- Outputs are saved to organized subdirectories.
- Embedded comments and section headers guide interpretation.
- Verbose output provides real-time progress updates.
