from .trajectory_estimation.sample_clustering import sample_clustering
from .trajectory_estimation.cyto_path import cytopath, estimate_cell_data
from .trajectory_estimation.cytopath_merger import cytopath_merger

def trajectories(data, smoothing=False, alignment_cutoff=0.0, basis='umap', neighbors_basis='pca', fill_cluster=True, n_neighbors_cluster='auto', cluster_freq=0.05,
                 groupby='mean', weighted=True, surrogate_cell=False, n_neighbors_alignment='auto', cluster_num=None, method = 'kmeans', distance='euclidean', num_cores=1, copy=False):

    """Trajectory inference using simulations on the transition proabability matrix.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    smoothing: Boolean (default:False)
        Whether or not to smooth over the trajectories.
    alignment_cutoff: float (default:0.0)
        Minimum alignment score to consider a cell-step alignment valid
    basis: str/list (default: umap)
        The space in which the neighboring cells should be searched. If None imputed expression is used. If list, imputed expression from supplied genes is used.
    neighbors_basis: str (default: pca)
        The space used to find neighboring cells of trajectory steps.
    fill_cluster: `Boolean` (default:True)
        Enforce only cells in compostional clusters are assigned score.
    n_neighbors_cluster: `integer` (Default:30)
        Number of cells to consider for determining compositional clusters
    cluster_freq: `float` (Default: 0.05)
        Frequency of cell cluster cells per step to consider as compositonal cluster
    groupby:  `str` (default: mean)
        One of max, min, median or mean. Grouping method for determing alignment score per cell per trajectory
    weighted: `Boolean` (default:False)
        If groupby is mean, reports mean pseduotime weighted by alignment score if true.
    surrogate_cell: `Boolean` (default:False)
        Whether or not a surrogate cell should be used for the neighborhood search
    n_neighbors_alignment:  `str/int/float` (default:'auto')
        Number of neighbors to searched along the average trajectory.
    cluster_num:  list (default:None)
        Number of trajectories (clusters) to be expected for each terminal region.
    method: str (default:'kmeans'):
        Which clustering method to use.
    distance: str (default:'cosine'):
        Which distabce metric to use (see py-hausdorff for details)
    num_cores: 'integer' (default:1)
        Number of cpu cores to use.
    copy: 'Boolean' (default: False)
        Create a copy of the anndata object or work inplace.
    Returns
    -------
    adata.uns["trajectories_2"]: Dictionary of the average trajectories
    adata.uns["trajectories"]: Long format of the average trajectories
    adata.uns["trajectories_count"]: Number of trajectories contributing to each cluster
    adata.uns['trajectories']["cells_along_trajectories"]: List of arrays, which denote the average step for cell.
    adata.uns['trajectories']["cells_along_trajectories_each_step"]: List of arrays containing the cell indexes for each step
    """

    adata = data.copy() if copy else data
    
    if 'X_pca' not in adata.obsm.keys(): raise ValueError('Compute PCA using cytopath.scvelo.pp.pca')

    # Clustering of cell sequences to obtain trajectories
    sample_clustering(adata, basis=basis, cluster_num=cluster_num, 
                      distance=distance, num_cores=num_cores)

    # Align cells along trajectories
    cytopath(adata, basis=basis, neighbors_basis=neighbors_basis, surrogate_cell=surrogate_cell, fill_cluster=fill_cluster, 
             cluster_freq=cluster_freq, n_neighbors_cluster=n_neighbors_cluster,
             n_neighbors=n_neighbors_alignment, cut_off=alignment_cutoff, num_cores=num_cores)
    
    # Estimate pseudotime, cell fate prob and alignment score at cell level
    estimate_cell_data(adata, groupby=groupby, weighted=weighted)

    return adata if copy else None

def merge_steps(data, overlap=0.0, num_cores=1, copy=False):

    adata = data.copy() if copy else data
    cytopath_merger(adata, overlap=overlap, num_cores=num_cores)   
    return adata if copy else None
