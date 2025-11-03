import numpy as np
import pandas as pd
import sys
import glob
import re
import os

dir = sys.argv[1]
print(dir)

# Pattern to match the files
file_pattern = '*.gene_spectra_tpm.k_*.dt_*.txt'

# Get a list of files matching the pattern
spectra_tpm_all_files = glob.glob(os.path.join(dir, file_pattern))
tpm_stats_file = glob.glob(os.path.join(dir, 'cnmf_tmp/*.tpm_stats.df.npz'))[0]

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj

def renormalize_tpm(spectra_tpm_all_files, tpm_stats_file, output_dir):
    # Load the common standard deviations
    stds = load_df_from_npz(tpm_stats_file)['__std']

    for spectra_file in spectra_tpm_all_files:
        # Load spectra_tpm_all
        tpm = pd.read_csv(spectra_file, index_col=0, sep='\t')
        # Normalize the TPM values by sum and scale to 1e6
        normalized_tpm = tpm.div(tpm.sum(axis=1), axis=0) * 1e6
        # Apply standard deviation scaling using the common `stds`
        renorm_tpm = normalized_tpm.div(stds)
        # print(renorm_tpm)
        k = re.search(r'\.k_(\d+)\.dt', spectra_file).group(1)
        dt = re.search(r'\.dt_(.*)\.txt', spectra_file).group(1)

        # Save the renormalized TPM DataFrame
        output_file = f'{output_dir}/{k}.{dt}.txt'
        renorm_tpm.to_csv(output_file, sep='\t')
        print(f'Renormalized TPM file saved to {output_file}')

# Directory for output files
output_dir = '.'

# Call the function
renormalize_tpm(spectra_tpm_all_files, tpm_stats_file, output_dir)
