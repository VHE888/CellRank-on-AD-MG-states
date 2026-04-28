# CellRank-on-AD-MG-states

**CellRank trajectory analysis for microglial state transitions in Alzheimer's disease**

---

## Overview

This project applies CellRank to study how microglial cells change their state during Alzheimer's disease progression. Using single-cell RNA sequencing data and RNA velocity, we model the dynamic transitions from homeostatic (resting) microglia to disease-associated states.

---

## Pipeline Steps

### 1. Data Loading & Preprocessing
Load the pre-processed AnnData object containing microglial single-cell data with spliced and unspliced count matrices required for RNA velocity computation.

### 2. RNA Velocity Computation
Calculate RNA velocity using scVelo's dynamical model. This infers the future transcriptional state of each cell based on the ratio of unspliced to spliced mRNA, providing directional information on cellular transitions.

### 3. Kernel Construction
Build a transition matrix by combining two components:
- **Velocity kernel** (50%): Based on RNA velocity vectors
- **Connectivity kernel** (50%): Based on cellular neighborhood similarity

The combined kernel defines the probability of moving from one cell to another.

### 4. Macrostate Identification (GPCCA)
Apply Generalized Perron Cluster Cluster Analysis to compress the complex cellular manifold into discrete macrostates. Each macrostate represents a distinct microglial biological condition (e.g., homeostatic, early DAM, late DAM, proliferative, inflammatory). The method also predicts:
- **Initial states** (sources): Where cells begin their trajectory
- **Terminal states** (sinks): Final states cells transition toward

### 5. Fate Probability Calculation
Compute the probability of each individual cell transitioning toward each terminal state. This reveals which cells are "committed" to becoming DAM, proliferative, or inflammatory microglia versus those still "undecided."

### 6. Gene Expression Trends
Identify genes that change expression along each trajectory (e.g., homeostatic → DAM). Genes are ranked by how strongly their expression correlates with pseudotime, revealing key drivers of state transitions.

### 7. Regulatory Analysis (Optional)
If metabolic labeling data is available, identify upstream transcription factors driving microglial fate decisions. This pinpoints potential master regulators of AD-associated microglial states.

---

## AD Microglial States Studied

| State | Description |
|-------|-------------|
| **Homeostatic** | Resting/surveying microglia; initial starting state |
| **Early DAM** | Transitional disease-associated microglia |
| **Late DAM** | Terminal disease-associated microglia with phagocytic markers |
| **Proliferative** | Cell cycle-active microglia |
| **Inflammatory** | Pro-inflammatory cytokine-producing microglia |

**Expected trajectory direction:** Homeostatic → Early DAM → Late DAM (with alternative branches toward proliferative and inflammatory states)

---

## Key Outputs

| Output | What it tells you |
|--------|-------------------|
| **Macrostate UMAP** | Which microglial states exist and how they relate spatially |
| **Fate probability maps** | The likelihood of each cell becoming each terminal state |
| **Gene trends** | Which genes drive each trajectory (e.g., *Apoe*, *Trem2*, *Gpnmb* for DAM) |
| **Regulatory networks** | Which transcription factors control state transitions |

---

## Requirements

- Python 3.8+
- CellRank 2.0+
- Scanpy
- scVelo

---

## References

- CellRank: Lange et al., *Nature Methods* (2022)
- scVelo: Bergen et al., *Nature Biotechnology* (2020)
- AD microglial states: Lee et al., *Cell* (2024)

---