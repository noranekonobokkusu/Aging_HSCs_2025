import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF, Preprocess
import anndata as ad
import warnings
import sys
import seaborn as sns

np.random.seed(1)
seed_n = 1

sample = sys.argv[1]
klow = int(sys.argv[2])
khigh = int(sys.argv[3])
dataset = sample.split('_')[1]

# Read h5ad file and subset the dataset
counts_filename = '/broad/vangalenlab/safina/projects/aging_hsc/prepare_for_publication/3_scVI_integration/seu_hsc.var_genes.h5ad'
adata_counts = ad.read_h5ad(counts_filename)
adata = adata_counts[adata_counts.obs['Sample']==sample].copy()

if dataset == 'Zhang2022':
    quantile_value_min = 100
    quantile_value_max = 10000
else:
    total_counts = adata.X.sum(axis=1).A1
    med = np.median(total_counts)
    mad = median_abs_deviation(total_counts)
    quantile_value_min = int(med - 2.5 * mad)
    quantile_value_max = int(med + 2.5 * mad)

sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_cells(adata, min_counts=quantile_value_min)
sc.pp.filter_cells(adata, max_counts=quantile_value_max)
sc.pp.filter_genes(adata, min_cells=10)

adata_counts_name = 'seu.h5ad'
adata.write(adata_counts_name)

ngenes = adata.n_vars
ks = np.arange(klow,khigh+1)

cnmf_obj = cNMF(output_dir= '.', name= 'sample')
cnmf_obj.prepare(counts_fn= adata_counts_name,
                           components= ks,
                           n_iter= 500,
                           seed= seed_n,
                           num_highvar_genes= ngenes)

cnmf_obj.factorize_multi_process(total_workers=32)
cnmf_obj.combine()
cnmf_obj.k_selection_plot()

for k_i in ks:
    for thr in [0.05, 0.1, 0.15]:
        try:
            print(k_i, thr, 'sample')
            cnmf_obj.consensus(k=k_i, density_threshold=thr, show_clustering=True, close_clustergram_fig=True)
            plt.close()

        except Exception as e:
            # Skip to the next iteration without checking the error message
            print(f"Skipping k={k_i}, thr={thr} due to an error: {e}")
            continue  # Skip to the next iteration of `thr`
