import scvi
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import pandas as pd

scvi.settings.num_threads = 1

# Load and preprocess data
adata_path = '../seu_full_object.h5ad'
adata = sc.read(adata_path)
print(adata.shape)
nlayers = 4

adata.obs['data_type'] = 'Single-cell'
# Specifically mark the single-nucleus dataset
adata.obs.loc[adata.obs['Dataset'] == 'Chen2024', 'data_type'] = 'Single-nucleus'
print(adata.obs['data_type'].unique())

# Avoid identical sample names in different datasets
adata.obs['Sample'] =  adata.obs['Sample'] + '_' + adata.obs['Dataset']

# Set up scVI
scvi.model.SCVI.setup_anndata(adata, batch_key='Sample', categorical_covariate_keys=['Dataset', 'data_type'])

# Train the scVI model
model = scvi.model.SCVI(adata, n_layers=int(nlayers), n_latent=30)
model.train(max_epochs=60)

# Get batch-corrected latent space
adata.obsm['X_scVI'] = model.get_latent_representation()

# Generate UMAP after integration
sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=15)

resolutions = [0.5, 0.8, 1, 1.2, 1.5]  # Adjust based on dataset

# Run Leiden clustering for each resolution and save in adata.obs
for res in resolutions:
    cluster_key = f'leiden_{res}'  # Create a key based on resolution
    sc.tl.leiden(adata, resolution=res, key_added=cluster_key)

sc.tl.umap(adata)
plt.figure(figsize=(8, 6))
sc.pl.umap(adata, color=['Dataset', 'predicted_CellType'], title='UMAP After Integration', show=False)
plt.savefig('umap_after_integration.png')
plt.close()

adata.write('integrated_seu.h5ad')
print('Integration completed.')
print(f'Number of genes used: {adata.shape[1]}')

# Extract metrics
reconstruction_loss_train = model.history['reconstruction_loss_train']
elbo_train = model.history['elbo_train']

# Create a plot
plt.figure(figsize=(8, 6))
plt.plot(reconstruction_loss_train, label='Reconstruction Loss (Train)', color='blue')
plt.plot(elbo_train, label='ELBO (Train)', color='green')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Reconstruction Loss and ELBO (Training)')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('reconstruction_and_elbo_train_plot.png', dpi=300)
plt.close()
