import pandas as pd
import numpy as np
import scanpy as sc
import loompy as lp
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
import sys
from scipy.stats import pointbiserialr

dataset = sys.argv[1]

# Read the input and the AUCell output loom files
adata = sc.read_loom( 'input_loom_objects/' + dataset + '.seu.loom' )
f_pyscenic_output = 'aucell_out/' + dataset + '.aucell.thr_0.01.loom'
lf = lp.connect(f_pyscenic_output, mode='r+', validate=False )

# Extract AUC matrix
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
lf.close()

# Save unfiltered (all cohorts present) AUC matrix to a file
auc_mtx_df = pd.DataFrame(auc_mtx)
auc_mtx_df.to_csv('final_loom_files/' + dataset + '.auc_mtx.thr_0.01.csv', index=True)

# Load regulons from ctx output and save regulon matrix with target gene weights
sig = load_signatures('ctx_out/' + dataset + '.thr_0.01.csv')
genes = adata.var.index.tolist()
regulon_dict = { sig[i].name : [ sig[i].gene2weight[gene] if gene in sig[i].gene2weight.keys() else 0 for gene in genes ] for i in np.arange(len(sig)) }
target_gene_weights = pd.DataFrame(regulon_dict, index=genes)
output_file = 'final_loom_files/' + dataset + '.regulons.thr_0.01.csv'
target_gene_weights.to_csv(output_file, index=True)

# Create a pandas object with cell annotations
cellAnnot = pd.concat(
    [
        pd.DataFrame( adata.obs.Cohort, index=adata.obs.Cell ),
        pd.DataFrame( adata.obs.predicted_CellType, index=adata.obs.Cell ),
        pd.DataFrame( adata.obs.Sample, index=adata.obs.Cell ),
        pd.DataFrame( adata.obs.nGene, index=adata.obs.Cell ),
        pd.DataFrame( adata.obs.nUMI, index=adata.obs.Cell ),
    ],
    axis=1
)

# Subset two relevant cohorts to compute association of AUC scores with Young/Aged status
condition = cellAnnot['Cohort'].isin(['Young', 'Aged'])
adata_subset = adata[adata.obs['Cohort'].isin(['Young', 'Aged'])].copy()
adata_subset = adata_subset[adata_subset.obs.sort_values(by='Cohort', key=lambda x: x.map({'Young': 0, 'Aged': 1})).index, :]
cells_subset = adata_subset.obs_names.tolist()
subset_auc_mtx = auc_mtx.loc[cells_subset]
subset_auc_mtx.index.tolist()==cells_subset

# Compute point-biserial (Pearson) correlations
binary_cohort_labels = adata_subset.obs['Cohort'].map({'Young': 0, 'Aged': 1})
# Convert to DataFrame
auc_df = pd.DataFrame(subset_auc_mtx, index=adata_subset.obs_names)
correlations = auc_df.corrwith(binary_cohort_labels, axis=0)

# Select top 10 by absolute correlation
top10_genes = correlations.dropna().sort_values(ascending=False).index.tolist()[0:10]

# Define AP-1 members
ap1_complex = [i + '(+)' for i in ['FOS', 'FOSB', 'JUN', 'JUNB', 'JUND', 'FOSL1', 'FOSL2', 'ATF2', 'ATF3', 'ATF4', 'ATF6', 'ATF7', 'BATF', 'BATF2', 'BATF3', 'MAFA', 'MAFB', 'MAF', 'MAFG', 'MAFF', 'MAFK']]

# Compute Z-scores of AUC matrix
auc_mtx_Z = pd.DataFrame( index=subset_auc_mtx.index )
auc_mtx_Z = (subset_auc_mtx - subset_auc_mtx.mean()) / subset_auc_mtx.std(ddof=0)

# Save Z-scores computed for cells in Young and Aged cohorts
auc_mtx_Z_df = pd.DataFrame(auc_mtx_Z)
auc_mtx_Z_df.to_csv('final_loom_files/' + dataset + '.Z_score_auc_matrix.thr_0.01.csv', index=True)

# Define legend-plotting function
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation='nearest', aspect='auto')
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

# Plot legends
cats_sm = sorted(list(set(adata_subset.obs['Sample'])))
colors = sns.color_palette('Spectral',n_colors=len(cats_sm) )
colorsd_sm = dict( zip( cats_sm, colors ))
sns.set()
sns.set(font_scale=0.5)
fig = palplot( colors, [i.replace('_' + dataset, '') for i in cats_sm], size=1.0)
plt.savefig('figures/' + dataset + '.sample_legend.pdf', dpi=600, bbox_inches = 'tight')

cats_ch = sorted(list(set(adata_subset.obs['Cohort'])))
colors = sns.color_palette('Paired',n_colors=len(cats_ch) )
colorsd_ch = dict( zip( cats_ch, colors ))
sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats_ch, size=1.0)
plt.savefig('figures/Cohort_legend.pdf', dpi=600, bbox_inches = 'tight')

# Define cell colors for heatmaps
row_colors_sample = [ colorsd_sm[x] for x in adata_subset.obs['Sample'] ]
row_colors_cohort = [ colorsd_ch[x] for x in adata_subset.obs['Cohort'] ]

# Plot Z-scores of AUC values for AP-1 members
data_ap1 = auc_mtx_Z[[i for i in ap1_complex if i in auc_mtx_Z.columns]]
sns.set()
sns.set(font_scale=1.2)
g = sns.clustermap(data_ap1, annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=[row_colors_sample, row_colors_cohort], row_cluster=False,col_cluster=True,
    cmap='YlGnBu', figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig('figures/' + dataset + '.ap1.hm.thr_0.01.pdf', dpi=600, bbox_inches = 'tight')

# Plot Z-scores of AUC values for top-10 aging regulons
data_top10 = auc_mtx_Z[[i for i in top10_genes if i in auc_mtx_Z.columns]]
sns.set()
sns.set(font_scale=1.2)
g = sns.clustermap(data_top10, annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=[row_colors_sample, row_colors_cohort], row_cluster=False,col_cluster=True,
    cmap='YlGnBu', figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig('figures/' + dataset + '.top10_aging.hm.thr_0.01.pdf', dpi=600, bbox_inches = 'tight')
