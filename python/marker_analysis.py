import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import episcanpy as epi



def top_markers(adata, ntop=100):
    """
    Find the top biomarkers in a list of differentially distributed
    features in an AnnData object from scanpy.tl.rank_genes_groups().
    - adata: AnnData object.
    - ntop: Number of top markers to keep
    """

    markers = None

    # If processed with Scanpy
    if 'rank_genes_groups' in adata.uns:
        markers = pd.DataFrame(
            adata.uns['rank_genes_groups']['names']
        ).head(ntop)

    # Else if processed with Episcanpy
    elif 'rank_features_groups' in adata.uns:
        markers = pd.DataFrame(
            adata.uns['rank_features_groups']['names']
        ).head(ntop)
    else:
        print('No rank_genes_groups variable found in the submitted ',
              'anndata. Please run scanpy.tl.rank_genes_groups() to ',
              'find group markers.')

    return(markers)



def plot_jaccardindex(data, xlab='x', ylab='y'):
    """
    Plot the jaccard index matrix as a heatmap.
    """

    data = data.round(decimals=2)
    # This will use the same order for the labels
    data = data.sort_index(0).sort_index(1)
    # This will reverse the matrix so the diagonal looks good
    data = data[::-1]

    plot = sns.heatmap(data, annot = True, vmin = 0.,
                       vmax = 1., cmap = 'Reds')
    plt.xticks(rotation = 70)
    plt.xlabel(xlab, fontsize = 30)
    plt.ylabel(ylab, fontsize = 30)
    plt.title('Jaccard Index scores', fontsize = 40, pad = 20)
    plt.show()

    return(plot)



def matchscore(ref_markers, obs_markers, plot=True,
              xlab='My dataset', ylab='Reference dataset'):
    """
    Analysis of similarity between a set of markers found in
    a reference dataset and a set of markers from the observed
    query dataset.
    - ref_markers: List of markers from the reference
    - obs_markers: List of markers from the observed query
    - plot: Whether or not to plot the jaccard index matrix
    """

    score = 0
    mat = []

    for ref_cluster in ref_markers:
        ref_len = len(ref_markers[ref_cluster])

        jacc_index = []

        for obs_cluster in obs_markers:

            obs_len = len(obs_markers[obs_cluster])
            shared = len(set(ref_markers[ref_cluster]).intersection(
                obs_markers[obs_cluster]))
            jacc = shared / (ref_len + obs_len - shared)
            jacc_index.append(jacc)

        score += np.max(jacc_index)
        mat.append(jacc_index)

    mat = pd.DataFrame(mat, index = list(ref_markers.columns),
                       columns = list(obs_markers.columns))

    score = score / ref_markers.shape[1]
    max_ji = mat.max()
    anno_lab = mat.idxmax()
    if plot:
        plot_jaccardindex(mat, xlab=xlab, ylab=ylab)

    return((score, anno_lab, max_ji, mat))



def find_common_genes(ref, sample, target_n_genes=2000):
    """
    Find a common set of features between two scRNA datasets
    using the highly variable genes.
    - ref: Reference AnnData object
    - sample: Query AnnData object
    - target_n_genes: Target number of common genes to find
    """

    print(f'Dimensions before filtering: {ref.shape} and {sample.shape}\n')

    n = 5000
    common = []
    while len(common) < target_n_genes:

        print(f'Looking for {n} HVG')
        sc.pp.highly_variable_genes(ref, n_top_genes=n)
        sc.pp.highly_variable_genes(sample, n_top_genes=n)

        ref_genes = ref.var.highly_variable[ref.var.highly_variable == True].index
        sample_genes = sample.var.highly_variable[sample.var.highly_variable == True].index

        common = list(set(sample_genes).intersection(ref_genes))
        print(f'Found {len(common)} genes in common')
        n += 1000

    print('Filtering the data to these features and scaling')

    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print(f'\nDimensions after filtering: {ref.shape} and {sample.shape}')

    return ref, sample



def find_common_variable_peaks(ref, sample, target_n_peaks=5000):
    """
    Find a common set of features between two scATAC datasets
    using the highly variable peaks.
    - ref: Reference AnnData object
    - sample: Query AnnData object
    - target_n_genes: Target number of common peaks to find
    """

    print('Dimensions before filtering: ',
          f'{ref.shape} and {sample.shape}\n')

    n = target_n_peaks
    common = []
    while len(common) < target_n_peaks:

        print(f'Looking for {n} HVP')
        ref_var = epi.pp.select_var_feature(
            ref, nb_features=n, show=False, copy=True)
        sample_var = epi.pp.select_var_feature(
            sample, nb_features=n, show=False, copy=True)

        ref_genes = list(ref_var.var.index)
        sample_genes = list(sample_var.var.index)

        common = list(set(sample_genes).intersection(ref_genes))
        print(f'Found {len(common)} peaks in common')
        n += 1000

    print('Filtering the data to these features and scaling')

    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print('Dimensions after filtering: ',
          f'{ref.shape} and {sample.shape}\n')

    return ref, sample



def find_refmarkers_in_variable_peaks(ref, sample, target_n_peaks=5000):
    """
    Find a common set of features between two scATAC datasets
    using the marker peaks from the reference and the highly
    variable peaks from the query sample.
    - ref: Reference AnnData object
    - sample: Query AnnData object
    - target_n_genes: Target number of common peaks to find
    """

    print('Dimensions before filtering: ',
          f'{ref.shape} and {sample.shape}\n')

    dfn = ref.obs['predicted.id'].value_counts() > 20
    keep_celltypes = dfn[dfn == True].index.to_list()
    ref = ref[ref.obs['predicted.id'].isin(keep_celltypes)]
    sc.pp.normalize_total(ref)
    epi.pp.log1p(ref)

    sc.pp.normalize_total(sample)
    epi.pp.log1p(sample)

    n_cat = len(ref.obs['predicted.id'].cat.categories)
    n = int(round(target_n_peaks*2/n_cat))
    common = []

    while len(common) < target_n_peaks:

        print(f'Looking for {n} marker features from each label')
        epi.tl.rank_features(ref, 'predicted.id', omic="ATAC",
                             use_raw=False, n_features=n)
        marker_peaks = list(ref.uns['rank_features_groups']['names'])
        marker_peaks = list(itertools.chain.from_iterable(marker_peaks))
        print(len(marker_peaks))

        sample_peaks = epi.pp.select_var_feature(
            sample, nb_features=len(marker_peaks), show=False, copy=True)
        sample_peaks = list(sample_peaks.var.index)

        common = list(set(marker_peaks).intersection(sample_peaks))
        print(f'Found {len(common)} peaks in common')
        n += 100

    print('Filtering the data to these features and scaling')

    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print('Dimensions after filtering: ',
          f'{ref.shape} and {sample.shape}\n')

    return ref, sample
