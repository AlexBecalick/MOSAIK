#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ========================
# Imports & Setup
# ========================
import os
import warnings
import logging

# Set the working directory where your data and results are located.
# Replace "/path/to/project/" with the absolute or relative path to your project folder.
path = "/Volumes/SSD/MERFISH/"
os.chdir(path)

# Silence unnecessary warnings and set logging to show only warnings or errors.
logging.basicConfig(level=logging.WARNING)
warnings.filterwarnings("ignore")

from spatialdata.transformations import Affine, set_transformation
import matplotlib.pyplot as plt
from sbf import visualise_fov
import spatialdata as sd
import spatialdata_plot
import seaborn as sns
import scanpy as sc
import pandas as pd
import merscope


# ========================
# Paths & Data Loading
# ========================
# zarr_path:
#   Path to the Zarr-formatted spatial dataset from the current path (established before). 
#   If the code is running for the first time you should give a name to the file 
#   and how to access from current path
zarr_path = "AG07.zarr"

# slide:
#   Path to the slide-specific directory (often a flattened structure used for annotations,
#   cell labels, or per-slide data). This typically corresponds to a single TMA or slide.
slide = "/AG07"

first_run = user_input = input("Is it the first run (0: False, 1: True): ")

if first_run == '1':
    sdata = merscope.merscope(path+slide, z_layers=[2,3])
    sdata.write(zarr_path)

sdata = sd.read_zarr(zarr_path)
adata = sdata.tables["table"]
print(adata.obs.keys())

# ========================
# Basic Spatial Plot
# ========================
xy = adata.obsm['spatial']
plt.scatter(xy[:, 0], xy[:, 1], s=0.0001)

# ========================
# Quality Control Metrics
# ========================
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

# ========================
# QC Plots
# ========================
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs["volume"], kde=False, ax=axs[2])
axs[2].set_title("Volume of segmented cells")
sns.histplot(adata.obs["perimeter_area_ratio"], kde=False, ax=axs[3])
axs[3].set_title("Perimeter area ratio")

# ========================
# Filtering & Normalization
# ========================
print("Original dimension: ", adata.shape)
sc.pp.filter_cells(adata, min_counts = 100)
print("Dimension after filtering cells: ", adata.shape)
sc.pp.filter_genes(adata, min_cells = 100)
print("Dimension after filtering genes: ", adata.shape)


adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# ========================
# Dimensionality Reduction & Clustering
# ========================
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# UMAP & spatial plots
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "leiden"], wspace=0.4, save=True)

adata.obs["x_global_px"] = adata.obsm['spatial'][:,0]
adata.obs["y_global_px"] = adata.obsm['spatial'][:,1]

g = sns.scatterplot(x="x_global_px", y="y_global_px", s=2, marker='.', 
                    data=adata.obs, hue='leiden', palette = "Set2")
sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
handles, labels = g.get_legend_handles_labels()
for h in handles:
    sizes = h.get_markersize()*8
    h.set_markersize(sizes)
plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, ncol=2)
g.set_ylabel("")
g.set_xlabel("")
plt.tight_layout()
plt.savefig('Sample_display_transcripts_cluster.png', format = 'png', dpi = 600)

