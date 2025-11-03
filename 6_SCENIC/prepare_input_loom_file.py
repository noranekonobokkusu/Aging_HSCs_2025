import pandas as pd
import numpy as np
import scanpy as sc
import loompy as lp
import os

if not os.path.exists('input_loom_objects'):
    os.makedirs('input_loom_objects')

adata_full = sc.read_h5ad('input_data/seu_hsc.total.h5ad')
print(adata_full.shape)
datasets = adata_full.obs['Dataset'].unique()

for dataset in datasets:
    print(dataset)
    adata = adata_full[adata_full.obs['Dataset'] == dataset].copy()

    row_attrs = {
        'Gene': np.array(adata.var.index) ,
    }
    col_attrs = {
        'CellID': np.array(adata.obs.index) ,
        'nGene': np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
        'nUMI': np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    }

    for key,values in adata.obs.items():
        col_attrs[key]=np.array(values)

    lp.create( 'input_loom_objects/' + dataset + '.seu.loom', adata.X.transpose(), row_attrs, col_attrs )
