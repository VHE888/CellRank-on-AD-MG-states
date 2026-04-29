import cellrank as cr
import scanpy as sc
import matplotlib.pyplot as plt

# ===============================
# Load data (with CellRank results)
# ===============================
adata = sc.read_h5ad("data/microglia_adata_cellrank.h5ad")

# If you don't have the saved file, run pipeline first:
# adata = sc.read_h5ad("data/microglia_adata.h5ad")
# ... then re-run CellRank analysis ...

# ===============================
# Figure 1: RNA Velocity Stream
# ===============================
print("Generating Figure 1: RNA Velocity Stream...")
fig, ax = plt.subplots(figsize=(10, 8))
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    color='cell_type',
    title='RNA Velocity - Microglial States',
    ax=ax,
    size=50,
    alpha=0.6
)
plt.savefig("figures/velocity_stream.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: figures/velocity_stream.png")

# ===============================
# Figure 2: Macrostates on UMAP
# ===============================
print("Generating Figure 2: Macrostates...")
# Note: This requires you have run GPCCA already
fig, ax = plt.subplots(figsize=(10, 8))
estimator.plot_macrostates(
    which='all',
    size=40,
    title='Microglial Macrostates in AD',
    ax=ax,
    legend_loc='right margin'
)
plt.savefig("figures/macrostates.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: figures/macrostates.png")

# ===============================
# Figure 3: Fate Probabilities
# ===============================
print("Generating Figure 3: Fate Probabilities...")
fig = estimator.plot_fate_probabilities(
    same_plot=False,
    n_cols=3,
    figsize=(15, 12),
    size=40
)
plt.savefig("figures/fate_probabilities.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: figures/fate_probabilities.png")

# ===============================
# Figure 4: Gene Trends (Top 5 genes per lineage)
# ===============================
print("Generating Figure 4: Gene Trends...")
lineages = [f"{estimator.initial_states[0]} -> {term}" for term in estimator.terminal_states]
trends = cr.estimators.GeneTrends(adata, estimator)
trends.prepare(lineages=lineages, time_key="pseudotime")
fig = trends.plot_top_genes(
    n_genes=5,
    lineages=lineages,
    figsize=(12, 8)
)
plt.savefig("figures/gene_trends.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: figures/gene_trends.png")

# ===============================
# Figure 5: Simple UMAP with cell types
# ===============================
print("Generating Figure 5: Cell Type UMAP...")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color='cell_type',
    title='Microglial Cell Types',
    ax=ax,
    size=30,
    legend_loc='on data'
)
plt.savefig("figures/celltype_umap.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved: figures/celltype_umap.png")

print("\n✅ All figures generated!")
print("Check the 'figures/' folder for outputs.")