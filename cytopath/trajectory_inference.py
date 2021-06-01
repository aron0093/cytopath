from .trajectory_estimation.sample_clustering import sample_clustering
from .trajectory_estimation.cyto_path import cytopath, estimate_cell_data
from .trajectory_estimation.cytopath_merger import cytopath_merger

def trajectories(data, smoothing=False, alignment_cutoff=0.0, basis='umap', neighbors_basis='pca', fill_cluster=True, n_neighbors_cluster=30, cluster_freq=0.1,
                 groupby='median', weighted=False, surrogate_cell=False, n_neighbors_alignment='auto', cluster_num=1, method = 'kmeans', distance = 'cosine', num_cores=1, copy=False):
    
    adata = data.copy() if copy else data
    
    if 'X_pca' not in adata.obsm.keys(): raise ValueError('Compute PCA using cytopath.scvelo.pp.pca')

    # Clustering of cell sequences to obtain trajectories
    sample_clustering(adata, basis=basis, smoothing=smoothing, 
                          cluster_num=cluster_num, method=method, distance=distance)

    # Align cells along trajectories
    cytopath(adata, basis=basis, neighbors_basis=neighbors_basis, surrogate_cell=surrogate_cell, fill_cluster=fill_cluster, cluster_freq=cluster_freq, n_neighbors_cluster=n_neighbors_cluster,
             n_neighbors=n_neighbors_alignment, cut_off=alignment_cutoff, num_cores=num_cores)
    
    # Estimate pseudotime, cell fate prob and alignment score at cell level
    estimate_cell_data(adata, groupby=groupby, weighted=weighted)

    return adata if copy else None

def merge_steps(data, overlap=0.0, num_cores=1, copy=False):

    adata = data.copy() if copy else data
    cytopath_merger(adata, overlap=overlap, num_cores=num_cores)   
    return adata if copy else None
