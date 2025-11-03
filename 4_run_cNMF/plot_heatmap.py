import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import glob
import re
import os

# Pattern to match the files
dir = './'
file_pattern = 'sample.usages.k_*.dt_0*.consensus.txt'

# Get a list of files matching the pattern
usage_files = glob.glob(file_pattern)

# Plot usage heatmaps
for i in usage_files:
    df = pd.read_table(i, sep='\t', header=0, index_col=0)
    df.fillna(0, inplace=True)
    df.replace('', 0, inplace=True)
    if df.sum().sum() == 0:
        continue
    df = df.div(df.sum(axis=1), axis=0)
    df = df.T
    plt.figure(figsize=(5, 4))
    sns.heatmap(df, cmap='Reds', vmin=0, vmax=1, xticklabels=False, yticklabels=True, cbar=False)
    plt.tight_layout()
    plt.savefig('postplots/' + i.replace('.consensus.txt', '').replace('sample.usages.', '') + '.heatmap.png', dpi=300)
