Usage
====

The entire process of inferring trajectories with Cytopath is divided into two steps. In the first step, a transition probability matrix is used to simulate possible differentiation paths. In the second step, the simulations are used to estimate trajectories from root states to each terminal region (cluster containing terminal states).  

Pre-processing
--------------

We assume that the user has already performed `velocity analysis <https://scvelo.readthedocs.io/>`_ at this point and has the processed anndata object in memory. Following is an example of an anndata processed with scvelo.

.. code-block:: console

   AnnData object with n_obs × n_vars = n_cells × n_genes
       obs: 'root_cells', 'end_points', 'louvain'
       var: 'velocity_gamma', 'velocity_qreg_ratio', 'velocity_r2', 'velocity_genes'
       uns: 'neighbors', 'pca', 'umap', 'velocity_graph', 'velocity_graph_neg', 'velocity_params'
       obsm: 'X_pca', 'X_umap', 'velocity_umap'
       varm: 'PCs'
       layers: 'Ms', 'Mu', 'ambiguous', 'matrix', 'spliced', 'unspliced', 'variance_velocity', 'velocity'
       obsp: 'connectivities', 'distances'
       
Store the transition probability matrix,
 
.. code-block:: console
 
   adata.uns['T_forward'] = scv.utils.get_transition_matrix(adata)
 
Markov chain sampling
---------------------

The sampling procedure requires a transition probability matrix (parameter: matrix_key), a clustering that will be used to define terminal regions (parameter: cluster_key) and root/terminal state probabilities (stored under root_cells and end_points obs keys). If a clustering is not specified then Louvain is used by default. 

The auto_adjust parameter will automatically and dynamically select values for technical parameters of the sampling process based on the properties of the dataset. We advise keeping this parameter set to True.

Additionally, set the number of vCPUs (parameter: num_cores) for significantly faster performance. If you wish to return a copy of the anndata object then set copy parameter to True.

.. code-block:: console
  
   import os
   import cytopath
   
   cytopath.sampling(adata, auto_adjust=True, matrix_key = 'T_forward', cluster_key = 'louvain', 
                     num_cores=os.cpu_count()-1, copy=False)
   
If the user wishes to manually specify root and terminal cell states, then use the following,

.. code-block:: console

   end_points = np.array(adata.obs.end_points > 0.99)[0] # Numerical index; use any selection criteria
   root_cells = np.array(adata.obs.root_cells > 0.99)[0] # Numerical index; use any selection criteria
   
   cytopath.sampling(adata, auto_adjust=True, matrix_key = 'T_forward', cluster_key = 'louvain', 
                     end_points=end_points, root_cells=root_cells, num_cores=os.cpu_count()-1, copy=False)
   
Finally, entire clusters can be designated as root or terminal regions as a means of incorporating biological know-how not reflected in velocity-based selection of root/terminal states

.. code-block:: console

   # Example designation of root/terminal clusters
   end_clusters = ['end_cluster_1', 'end_cluster_2'] 
   root_clusters = ['root_cluster_1']
   
   cytopath.sampling(adata, auto_adjust=True, matrix_key = 'T_forward', cluster_key = 'louvain', 
                     end_clusters=end_clusters, root_clusters=root_clusters, num_cores=os.cpu_count()-1, copy=False)
                                                     
Trajectory inference
--------------------

Trajectory inference with default parameters will require the anndata to contain a PCA embedding (uns: 'X_pca') and UMAP (uns: 'X_umap'). The latter can be substituted for any 2D embedding being used to visualise the data (parameter: basis).

.. code-block:: console

   cytopath.trajectories(adata, num_cores=os.cpu_count()-1)
   
Plotting
--------

.. code-block:: console
   
   cytopath.plot_trajectories(adata, basis='umap')
   
Inference output
----------------

The trajectories inferred by Cytopath are composed of segments. Cells are aligned to these trajectory segments to determine their relative position along the trajectory (pseudotime) and relative association with multiple trajectories (cell-fate). 

The complete inference output containing all cell-trajectory alignments is stored under the following key,

.. code-block:: console
   
   adata.uns['trajectories']['cells_along_trajectories_each_step']
   
Inference output summarised at the single cell level is stored under.

.. code-block:: console
   
   adata.uns['trajectories']['cells_along_trajectories']
   
Simulations are stored under,

.. code-block:: console
   
   adata.uns['samples']
   
Trajectory coordinates are stored under,

.. code-block:: console
   
   adata.uns['trajectories']
    
