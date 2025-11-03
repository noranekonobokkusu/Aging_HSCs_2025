# The functions below (correlate_geps, cluster_geps, pairwise_cluster_geps) are based on the code from *CAT repository (https://github.com/immunogenomics/starCAT/blob/main/src/starcat/build_consensus_reference.py) with slight modifications. correlate_geps() has been significantly simplified; pairwise_cluster_geps() now allows clustering of two programs from the same cNMF run (for some samples, more than one aging signature are identified, which cluster together)

import pandas as pd
import numpy as np
import os

def correlate_geps(m):
    m = m.drop(['Sample'], axis=1)
    R = m.T.corr()
    return(R)

def cluster_geps(R, corr_thresh, pct_thresh, nhits):
    geps = R.index
    # Define adjacency matrix between GEP pairs
    A = pd.DataFrame(np.zeros((len(geps), len(geps))), index = geps, columns = geps)
    dataset_ind = pd.Series([gep.split(':')[0] for gep in geps], index=geps)
    other_dataset_map = {x:dataset_ind.index for x in dataset_ind.unique()}
    # Define an edge for correlated GEPs not from the same cNMF result
    for gep in R.columns:
        ds = gep.split(':')[0]
        top_geps = R.loc[other_dataset_map[ds], gep].sort_values(ascending = False).head(nhits)
        top_geps = top_geps.loc[top_geps > corr_thresh]
        # Remove GEPs from the same dataset
        tokeep = top_geps.drop(gep).index
        A.loc[gep, tokeep] = top_geps.loc[tokeep]
    # Order by correlation
    ord_A = pd.DataFrame(A.unstack().sort_values(ascending = False)).reset_index()
    ord_A.columns = ['gep1', 'gep2', 'corr']
    ord_A = ord_A[ord_A['corr']!=0]
    # Cluster GEPs if edge exists between them
    gep_clusters = pd.DataFrame({'gep' : A.index, 'cluster' : 0})
    gep_clusters.index = gep_clusters['gep']
    for i in range(0, ord_A.shape[0]):
        gep = ord_A.loc[i, 'gep1']
        gep_alt = ord_A.loc[i, 'gep2']
        if A.loc[gep_alt, gep]>0:
            gep_clusters = pairwise_cluster_geps(gep, gep_alt, gep_clusters, A, pct_thresh)
    # Renumber clusters and assign singletons to unique clusters
    clus_dict = {clus_num:sorted(gep_clusters.loc[gep_clusters['cluster']==clus_num, 'gep'].values)
             for clus_num in sorted(gep_clusters['cluster'].unique())}
    if 0 in clus_dict:
        del clus_dict[0]
    clus_dict_clean = {}
    for clus in clus_dict.keys():
        new_clus = next(i for i, e in enumerate(sorted(clus_dict_clean.keys()) + [ None ], 1) if i != e)
        clus_dict_clean[new_clus] = clus_dict[clus]

    clustered = set(member for members in clus_dict_clean.values() for member in members)
    unclustered = [gep for gep in geps if gep not in clustered]

    data = []
    for cluster, members in clus_dict_clean.items():
        for member in members:
            data.append({'Cluster': cluster, 'Member': member})

    # Add unclustered items with a special 'Cluster' label
    for member in unclustered:
        data.append({'Cluster': 'Unclustered', 'Member': member})

    # Convert to DataFrame
    df = pd.DataFrame(data)

    # Recalculate cluster sizes excluding 'Unclustered'
    cluster_sizes = df[df['Cluster'] != 'Unclustered']['Cluster'].value_counts()

    # Map old cluster numbers to new ones based on size ranking
    sorted_clusters = cluster_sizes.sort_values(ascending=False).index.tolist()
    cluster_mapping = {old_cluster: str(new_cluster + 1) for new_cluster, old_cluster in enumerate(sorted_clusters)}

    # Update the DataFrame with new cluster numbering
    df['Cluster'] = df['Cluster'].replace(cluster_mapping)

    # Include 'Unclustered' as the last category for sorting
    all_clusters_sorted = sorted(cluster_mapping.values(), key=int) + ['Unclustered']
    df['Cluster'] = pd.Categorical(df['Cluster'], categories=all_clusters_sorted, ordered=True)

    # Sort the DataFrame
    df = df.sort_values(by='Cluster')
    return(df)

def pairwise_cluster_geps(gep, gep_alt, gep_clusters, A, pct_thresh):
    # Find which GEPs have been clustered already
    gep_num = gep_clusters.loc[gep, 'cluster']
    alt_gep_num = gep_clusters.loc[gep_alt, 'cluster']
    clustered = [gepX for (gepX, num) in [(gep, gep_num), (gep_alt, alt_gep_num)] if num > 0 ]
    # Merge 2 non-clusterd GEPs
    if len(clustered)==0:
        clus_num = gep_clusters['cluster'].max() + 1
        gep_clusters.loc[[gep, gep_alt], 'cluster'] = clus_num
    # Merge 2 different clusters if enough shared edges
    elif len(clustered)==2 and (gep_num != alt_gep_num):
        gep_list = sorted(gep_clusters.loc[gep_clusters['cluster']==gep_num, 'gep'].values)
        alt_gep_list = sorted(gep_clusters.loc[gep_clusters['cluster']==alt_gep_num, 'gep'].values)
        num_unconnected = (A.loc[gep_list, alt_gep_list]==0).sum().sum() + (A.loc[alt_gep_list, gep_list]==0).sum().sum()
        num_total = len(gep_list) * len(alt_gep_list) * 2
        pct_adj = 1 - (num_unconnected / num_total)
        gep_ds = [x.split(':', 1)[0] for x in gep_list]
        alt_gep_ds = [x.split(':', 1)[0] for x in alt_gep_list]
        if (pct_adj > pct_thresh):
            gep_clusters.loc[alt_gep_list, 'cluster'] = gep_num
    # Add unclustered GEP to a cluster if enough shared edges
    elif len(clustered)==1:
        clus_num = gep_clusters.loc[clustered[0], 'cluster']
        clus_list = sorted(gep_clusters.loc[gep_clusters['cluster']==clus_num, 'gep'].values)
        unclustered = [gepX for gepX in [gep, gep_alt] if gepX not in clustered]
        pct_adj = 1 - np.sum((np.array(A.loc[unclustered, clus_list] == 0).sum(),
                              (np.array(A.loc[clus_list, unclustered] == 0).sum())))/(len(unclustered)*len(clus_list)*2)
        clus_ds = [x.split(':', 1)[0] for x in clus_list]
        unclustered_ds = [x.split(':', 1)[0] for x in unclustered]
        if (pct_adj > pct_thresh):
            gep_clusters.loc[unclustered, 'cluster'] = clus_num
    elif len(clustered)==2 and (gep_num == alt_gep_num):
        pass
    else:
        sys.exit(-1)
    return(gep_clusters)

# Define thresholds
corr_thresh = 0.1
pct_thresh = 0.5

# Create a new directory for results
if not os.path.exists('meta_clusters_starcat'):
    os.makedirs('meta_clusters_starcat')

# Cluster programs, excluding singletons
singletons = pd.read_csv('singletons.txt', header=None)
singletons = singletons[0].tolist()
m = pd.read_csv('cnmf_output/Z_scores.total_cNMF_output.csv')
m = m[~ m.index.isin(singletons)]
R = correlate_geps(m)
nhits = R.shape[0]
df = cluster_geps(R, corr_thresh, pct_thresh, nhits)
df.to_csv('meta_clusters_starcat/clusters.csv', index=False)
