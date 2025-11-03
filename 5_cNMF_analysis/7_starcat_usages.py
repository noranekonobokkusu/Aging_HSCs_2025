import sys
import os
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from starcat import starCAT
import json

if not os.path.exists('starcat_usages'):
    os.makedirs('starcat_usages')

print("Inferring meta-program usages for each sample")

ref_file = 'meta_clusters_starcat/metaGEP.Z_score_based_TPMs.tsv'
starcat_ref = starCAT(reference=ref_file)
input_file = '/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/1_dataset_preparation/seu_hsc.var_genes.h5ad'
adata = ad.read_h5ad(input_file)
# samples = adata.obs['Sample'].unique()
# for sample in samples:
#     adata_subset = adata[adata.obs['Sample']==sample].copy()
#     adata_subset.var_names = adata_subset.var_names.str.replace('-', '.')
#     usage, _ = starcat_ref.fit_transform(adata_subset)
#     output_file = 'starcat_usages/' + sample + '.starcat_usage.csv'
#     usage.to_csv(output_file)

datasets = adata.obs['Dataset'].unique()
for dataset in datasets:
    adata_subset = adata[adata.obs['Dataset']==dataset].copy()
    adata_subset.var_names = adata_subset.var_names.str.replace('-', '.')
    usage, _ = starcat_ref.fit_transform(adata_subset)
    output_file = 'starcat_usages/' + dataset + '.starcat_usage.csv'
    usage.to_csv(output_file)
