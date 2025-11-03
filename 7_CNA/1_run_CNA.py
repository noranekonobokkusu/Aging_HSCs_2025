import pandas as pd
import scanpy as sc
import cna
import anndata as ad
# from multianndata import MultiAnnData

# Read input object
d = ad.read_h5ad("/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/3_scVI_integration/integrated_seu.h5ad")

# Subset Young and Aged samples
d = d[d.obs.Cohort.isin(['Young', 'Aged'])].copy()

# Convert relevant variables to numeric
col = 'Cohort'
d.obs[col + '_numeric'] = pd.Categorical(d.obs[col], categories=['Young', 'Aged'], ordered=True).codes

for col in ['Sample', 'Dataset']:
    d.obs[col + '_numeric'] = pd.Categorical(d.obs[col]).codes

# Ensure neighborhood graph is computed
sc.pp.neighbors(d, use_rep='X_scVI')
sc.tl.umap(d)

# Create sample-level metadata by combining across all cells from each sample
samplem = cna.ut.obs_to_sample(d,
                               ['Cohort_numeric','Dataset_numeric'],     # the metadata columns to use
                               'Sample_numeric')                        # the column of d.obs that contains sample IDs

# Run CNA
p = cna.tl.association(d,
            samplem.Cohort_numeric,                 #sample-level attribute of intest (case/control status)
            'Sample_numeric',                         #the column in d.obs that contains sample ids
            batches=samplem.Dataset_numeric,        #batch assignments for each sample so that cna can account for batch effects
            key_added='cohort_coef',
            nsteps=100000)

print('\nglobal association p-value:', p)
print((d.obs.cohort_coef_fdr <= 0.1).mean(), 'neighborhoods passed significance at FDR 10%')

# Save results as a file
cna_results = d.obs[['cohort_coef', 'cohort_coef_fdr']]
cna_results.to_csv('/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/7_CNA/cna_results.csv', index=True)
