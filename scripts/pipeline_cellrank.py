# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
import cellrank as cr
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set verbosity and plotting styles
cr.settings.verbosity = 2
scv.settings.figdir = "./figures/"
scv.settings.set_figure_params("scvelo", dpi=100, facecolor="white")

# Load your pre-processed AnnData object
# Replace with your own microglial dataset (e.g., .h5ad or .loom)
adata = sc.read_h5ad("data/microglia_adata.h5ad")
print(adata)

# ------------------------------------------------------------
# Phase 1: RNA Velocity & Kernel Construction
# ------------------------------------------------------------
print("Computing RNA velocity and transition kernel...")

# 1. Velocity dynamics
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# 2. Initialize CellRank Kernel (Velocity + Connectivity)
vk = cr.kernels.VelocityKernel(adata)
ck = cr.kernels.ConnectivityKernel(adata)
combined_kernel = 0.5 * vk + 0.5 * ck

# 3. Compute the transition matrix
combined_kernel.compute_transition_matrix()

# ------------------------------------------------------------
# Phase 2: Identifying Microglial States (Macrostates)
# ------------------------------------------------------------
print("Identifying microglial macrostates...")

# 1. Initialize the GPCCA estimator
estimator = cr.estimators.GPCCA(combined_kernel)

# 2. Fit the model
# n_states: Start with ~10-15, then refine based on microglial subtypes
estimator.fit(n_states=12, cluster_key="cell_type") 

# 3. Predict Terminal & Initial States
# In AD, terminal states include 'GPNMB+ DAM' or 'Proliferative' microglia
estimator.predict_terminal_states(method="top_n", n_states=3)
estimator.predict_initial_states()

# 4. Visualization: Macrostates on UMAP
fig, ax = plt.subplots(figsize=(8, 6))
estimator.plot_macrostates(which="all", size=40, title="Microglial Macrostates in AD", ax=ax)
plt.savefig("./figures/microglia_macrostates.png", dpi=300, bbox_inches="tight")

# ------------------------------------------------------------
# Phase 3: Fate Probability & Trajectory-Specific Dynamics
# ------------------------------------------------------------
print("Computing fate probabilities...")

# 1. Compute fate probabilities for every cell
estimator.compute_fate_probabilities()

# 2. Visualize fate probabilities on UMAP
# This shows how cells commit to specific AD trajectories
estimator.plot_fate_probabilities(same_plot=False, n_cols=3)
plt.savefig("./figures/fate_probabilities.png", dpi=300, bbox_inches="tight")

# 3. Extract gene expression trends along the microglial trajectory
print("Extracting gene trends along trajectories...")

# Example lineages (adjust based on your identified states)
lineages = ["Homeostatic -> DAM", "Homeostatic -> Proliferative"]

trends_model = cr.estimators.GeneTrends(adata, estimator)
trends_model.prepare(
    lineages=lineages, 
    time_key="pseudotime"
)

# ------------------------------------------------------------
# Phase 4: Differential Regulation & Hierarchy (AD Specific)
# ------------------------------------------------------------
# This step requires metabolic labeling or temporal data (e.g., ATAC-seq)
if "metabolic_label" in adata.obs.columns:
    print("Computing regulatory strategies...")
    
    # 1. Compute transcription and degradation rates
    estimator.compute_rates(data_key="X_spliced", time_key="experimental_time")
    
    # 2. Identify upstream transcription factors driving AD microglial states
    # e.g., MITF, KLF12, and others implicated in human AD
    regulator_dict = cr.ul.regulatory_network.estimate_transcriptional_regulation(
        adata,
        lineage="DAM",
        regulator_key="tf_activity"
    )
    
    # Optional: Save regulator results
    pd.DataFrame(regulator_dict).to_csv("./figures/dam_regulators.csv")
else:
    print("Skipping regulatory analysis: missing metabolic labeling data.")

print("Pipeline complete! Check the './figures/' directory for outputs.")