===============================
1. Load Libraries
===============================
import cellrank as cr
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set global settings
cr.settings.verbosity = 2
scv.settings.set_figure_params(dpi=100, facecolor="white")
sc.set_figure_params(dpi=100, figsize=(8, 6))

===============================
2. Load Input Data
===============================
# Load pre-processed AnnData object
# Replace with your file path
adata = sc.read_h5ad("data/microglia_adata.h5ad")

# Ensure required data is present
print(f"Dataset: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Cell types: {adata.obs['cell_type'].unique()}")

# Check for required layers (spliced/unspliced for velocity)
if 'spliced' not in adata.layers:
    print("Warning: 'spliced' layer not found. RNA velocity may fail.")
if 'unspliced' not in adata.layers:
    print("Warning: 'unspliced' layer not found. RNA velocity may fail.")

===============================
3. Preprocessing for RNA Velocity
===============================
# Filter and normalize for velocity analysis
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# Calculate moments for velocity estimation
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Optional: visualize variance explained
scv.pl.utils.get_n_neighbors(adata)

===============================
4. Compute RNA Velocity
===============================
# Dynamical model for velocity estimation
scv.tl.velocity(adata, mode="dynamical")

# Compute velocity graph for transitions
scv.tl.velocity_graph(adata)

# Optional: visualize velocity stream on UMAP
scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell_type', 
                                  title='RNA Velocity - Microglial States')
plt.savefig("figures/velocity_stream.png", dpi=300, bbox_inches="tight")

===============================
5. Build CellRank Transition Kernel
===============================
# Initialize velocity kernel (50% weight)
vk = cr.kernels.VelocityKernel(adata)

# Initialize connectivity kernel (50% weight)
ck = cr.kernels.ConnectivityKernel(adata)

# Combine kernels
kernel = 0.5 * vk + 0.5 * ck

# Compute transition matrix
kernel.compute_transition_matrix()

===============================
6. Identify Macrostates (GPCCA)
===============================
# Initialize GPCCA estimator
estimator = cr.estimators.GPCCA(kernel)

# Fit model to identify macrostates
# n_states: number of macrostates (adjust based on your data)
estimator.fit(n_states=12, cluster_key="cell_type")

# Predict terminal states (sinks) - where cells end up
estimator.predict_terminal_states(method="top_n", n_states=3)

# Predict initial states (sources) - where cells start
estimator.predict_initial_states()

print(f"Terminal states: {estimator.terminal_states}")
print(f"Initial states: {estimator.initial_states}")

===============================
7. Calculate Fate Probabilities
===============================
# Compute absorption probabilities for each terminal state
estimator.compute_fate_probabilities()

# Optional: visualize fate probabilities on UMAP
estimator.plot_fate_probabilities(same_plot=False, n_cols=3)
plt.savefig("figures/fate_probabilities.png", dpi=300, bbox_inches="tight")

===============================
8. Extract Gene Trends Along Trajectories
===============================
# Define lineages (initial -> terminal)
lineages = [f"{estimator.initial_states[0]} -> {term}" for term in estimator.terminal_states]
print(f"Analyzing lineages: {lineages}")

# Initialize gene trends model
trends = cr.estimators.GeneTrends(adata, estimator)

# Prepare trends along pseudotime
trends.prepare(lineages=lineages, time_key="pseudotime")

# Visualize top genes per lineage
trends.plot_top_genes(n_genes=5, lineages=lineages)
plt.savefig("figures/gene_trends.png", dpi=300, bbox_inches="tight")

===============================
9. Optional: Regulatory Analysis
===============================
# Only run if metabolic labeling data is available
if "metabolic_label" in adata.obs.columns:
    # Compute transcription and degradation rates
    estimator.compute_rates(data_key="X_spliced", time_key="experimental_time")
    
    # Identify transcription factors driving transitions
    regulators = cr.ul.regulatory_network.estimate_transcriptional_regulation(
        adata, lineage=estimator.terminal_states[0], regulator_key="tf_activity"
    )
    
    # Save results
    pd.DataFrame(regulators).to_csv("figures/dam_regulators.csv")
    print("Regulatory analysis complete!")
else:
    print("Skipping regulatory analysis: No metabolic labeling data found.")

===============================
10. Save Results
===============================
# Save updated AnnData with all CellRank results
adata.write("data/microglia_adata_cellrank.h5ad")
print("Results saved to data/microglia_adata_cellrank.h5ad")

# Print summary
print("\n" + "="*30)
print("PIPELINE COMPLETE")
print("="*30)
print(f"Cells analyzed: {adata.n_obs}")
print(f"Macrostates found: {len(estimator.macrostates)}")
print(f"Terminal states: {estimator.terminal_states}")
print(f"Output figures saved to: figures/")