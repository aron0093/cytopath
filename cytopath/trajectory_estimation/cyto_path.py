import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import spatial
from scvelo.preprocessing.neighbors import get_connectivities
from sklearn.metrics.pairwise import cosine_similarity
from joblib import Parallel, delayed

def surrogate_cell_neighborhood_finder(adata, end_point, map_state, fill_cluster=True, n_neighbors=30, cluster_freq=0.1, mode='distances', recurse_neighbors=True):
    """Finds surrogate cell and uses the scvelo neighborhood graph to finds its recursive neighbors. Applicable when mean is used to compute trajectory coordinates. Default is median.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    end_point: `string`
        End point cluster to which the trajectories belong
    map_state: matrix
        Coordinates of cells in selected basis.
    n_neighbors: `integer`(default: 30)
        Number of neighbors to consider

    Returns
    -------
    Returns list of arrays containing the sequence of neighborhoods along the trajectory and the sequence of surrogate cells.
    """
    # Locate surrogate cells to anchor trajectory
    final_cluster=adata.uns['trajectories']["trajectories_coordinates"][end_point]
    final_sequences=[]
    for i in tqdm(range(len(final_cluster))):
        final_cluster_sequence=final_cluster['trajectory_'+str(i)+'_coordinates']
        cell_sequences=[]
        for j in range(len(final_cluster_sequence)):
            cell_sequences.append(spatial.KDTree(map_state).query(final_cluster_sequence[j])[1])
        final_sequences.append(cell_sequences)

    # Neighborhood of surrogate cells
    neighborhood=get_connectivities(adata, mode=mode, n_neighbors=n_neighbors, recurse_neighbors=recurse_neighbors)
    neighborhood_sequence=[]
    for i in tqdm(range(len(final_sequences))):
        local_neighs=[]
        for j in range(len(final_sequences[i])):
            local_neighs.append(neighborhood[final_sequences[i][j]].nonzero()[1])
            
            if fill_cluster == True:
                print('Fill cluster not implemented.')

        neighborhood_sequence.append(local_neighs)
    return neighborhood_sequence, cell_sequences

def cell_neighborhood_finder(adata, map_state, end_point, neighbors_basis='pca', fill_cluster=True, n_neighbors_cluster=30, cluster_freq = 0.5, n_neighbors='auto'):
    """Finds the nearest neighbors along the average trajectory.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    end_point: `string`
        End point cluster to which the trajectories belong
    fill_cluster: Boolean (default:True)
        Enforce only cells in compostional clusters are assigned score.
    n_neighbors_cluster: `integer` (Default:30)
        Number of cells to consider for determining compositional clusters
    cluster_freq: `float` (Default: 0.25)
        Frequency of cell cluster cells per step to consider as compositonal cluster
    n_neighbors: `` (default: 'auto')
        Number of cells to consider for alignment scoring.

    Returns
    -------
    Returns list of arrays containing the sequence of neighborhoods along the trajectory.
    """

    # Cell neighborhood of trajectory
    final_cluster=adata.uns['trajectories']["trajectories_coordinates"][end_point]
    
    # Find clusters composing the trajectory
    compositional_clusters = []
    for i in tqdm(range(len(final_cluster))):
        final_cluster_sequence=final_cluster['trajectory_'+str(i)+'_coordinates']
        for j in range(len(final_cluster_sequence)):
            cell_sequences_ = spatial.KDTree(map_state).query(final_cluster_sequence[j], k=n_neighbors_cluster)[1]
            compositional_clusters_ = adata.obs.loc[adata.obs.index.values[cell_sequences_], adata.uns['run_info']['cluster_key']].astype(str)
            compositional_clusters_ = compositional_clusters_.value_counts()
            compositional_clusters_ /= sum(compositional_clusters_)
            compositional_clusters_ = compositional_clusters_.loc[compositional_clusters_ >= cluster_freq]
            compositional_clusters.extend(compositional_clusters_.index.tolist())
    adata.uns['trajectories']['compositional_clusters'][end_point] = np.unique(compositional_clusters)

    # Find neighbors in PCA space
    if neighbors_basis=='pca':

        # Use pca to find neigbhors if it exists
        map_state_ = adata.obsm['X_pca']

        final_cluster_ = {}
        for i in range(len(final_cluster)):
            final_cluster_sequence=final_cluster['trajectory_'+str(i)+'_coordinates']
            final_cluster_sequence_ = []
            for j in range(len(final_cluster_sequence)):
                final_cluster_sequence_.append(adata.obsm['X_pca'][spatial.KDTree(map_state).query(final_cluster_sequence[j])[1]])
            final_cluster_['trajectory_'+str(i)+'_coordinates'] = final_cluster_sequence_
        final_cluster = final_cluster_
    else:
        map_state_ = map_state

    
    # Neighborhood size
    if n_neighbors == 'auto':
        n_neighbors = adata.obs.loc[adata.obs[adata.uns['run_info']['cluster_key']].astype(str).isin(adata.uns['trajectories']['compositional_clusters'][end_point]), adata.uns['run_info']['cluster_key']].value_counts().max()

    # Create neighborhood
    neighborhood_sequence=[]
    for i in tqdm(range(len(final_cluster))):
        final_cluster_sequence=final_cluster['trajectory_'+str(i)+'_coordinates']
        cell_sequences=[]
        for j in range(len(final_cluster_sequence)):
            if fill_cluster == False:
                cell_sequences.append(spatial.KDTree(map_state_).query(final_cluster_sequence[j], k=n_neighbors)[1])
            elif fill_cluster == True:
                cell_sequences_ = spatial.KDTree(map_state_).query(final_cluster_sequence[j], k=n_neighbors)[1]
                cell_sequences.append(cell_sequences_[np.where(adata.obs.loc[adata.obs.index.values[cell_sequences_], adata.uns['run_info']['cluster_key']].astype(str).isin(adata.uns['trajectories']['compositional_clusters'][end_point]))[0]])
        neighborhood_sequence.append(cell_sequences)
    return neighborhood_sequence
    
def directionality_score(adata, neighborhood_sequence, end_point, map_state, num_cores=1): 

    """Scores the directionality of the cells along the path:
    Evaluation using the cosine angle between trajectory and neighboring cell to cell transitions and
    the transition probability
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    end_point: `string`
        End point cluster to which the trajectories belong
    map_state: matrix
        Coordinates of cells in selected basis.
    neighborhood_sequence: 
        List of the cells along the trajectory.
        
    Returns
    -------
    Returns list of arrays containing the scores of all cells at each step.
    """

    final_cluster=adata.uns['trajectories']["trajectories_coordinates"][end_point]
    
    transmat=adata.uns["T_forward"].copy()
    transmat[transmat.nonzero()]=1
    
    if 'velocity_graph' and 'velocity_graph_neg' in adata.uns.keys():
        velocity_graph=adata.uns["velocity_graph"]+adata.uns["velocity_graph_neg"]
        
    velocity_graph=velocity_graph.multiply(transmat)
    all_scores=[]
    def parallelize(j):
        score=np.zeros((len(neighborhood_sequence[i][j])))

        if j==0:
            avg_difference=final_cluster['trajectory_'+str(i)+'_coordinates'][j+1]-final_cluster['trajectory_'+str(i)+'_coordinates'][j]
            neighborhood_transitions=velocity_graph[neighborhood_sequence[i][j]]

            for k in range(len(neighborhood_sequence[i][j])):
                cell_jumps=neighborhood_transitions[k,:].nonzero()[1]
                origin_cell=map_state[neighborhood_sequence[i][j][k]]
                end_cells=map_state[cell_jumps]
                difference=origin_cell-end_cells

                # Score based on alignment and distance 
                score[k]=np.average(cosine_similarity(avg_difference.reshape(1,-1),difference)*np.exp(neighborhood_transitions[k,:].data))

        elif j!=(len(neighborhood_sequence[i])-1) and j!=0: 
            avg_difference_forward=final_cluster['trajectory_'+str(i)+'_coordinates'][j+1]-final_cluster['trajectory_'+str(i)+'_coordinates'][j]
            avg_difference_backward=final_cluster['trajectory_'+str(i)+'_coordinates'][j]-final_cluster['trajectory_'+str(i)+'_coordinates'][j-1]

            neighborhood_transitions=velocity_graph[neighborhood_sequence[i][j]]
            for k in range(len(neighborhood_sequence[i][j])):
                cell_jumps=neighborhood_transitions[k,:].nonzero()[1]
                origin_cell=map_state[neighborhood_sequence[i][j][k]]
                end_cells=map_state[cell_jumps]
                difference=origin_cell-end_cells

                # Score based on alignment and distance 
                score_forward=np.average(cosine_similarity(avg_difference_forward.reshape(1,-1), difference)*np.exp(neighborhood_transitions[k,:].data))
                score_backward=np.average(cosine_similarity(avg_difference_backward.reshape(1,-1), difference)*np.exp(neighborhood_transitions[k,:].data))
                
                score[k]=max(score_forward, score_backward)

        elif j==(len(neighborhood_sequence[i])-1):
            avg_difference_backward=final_cluster['trajectory_'+str(i)+'_coordinates'][j]-final_cluster['trajectory_'+str(i)+'_coordinates'][j-1]

            neighborhood_transitions=velocity_graph[neighborhood_sequence[i][j]]
            for k in range(len(neighborhood_sequence[i][j])):
                cell_jumps=neighborhood_transitions[k,:].nonzero()[1]
                origin_cell=map_state[neighborhood_sequence[i][j][k]]
                end_cells=map_state[cell_jumps]
                difference=origin_cell-end_cells

                # Score based on alignment and distance  
                score_backward=np.average(cosine_similarity(avg_difference_backward.reshape(1,-1), difference)*np.exp(neighborhood_transitions[k,:].data))
                score[k]=score_backward

        return (j, score)

    # Calculation of the directionality score for each cell in each trajectory
    for i in range(len(final_cluster)):
        cluster_scores = Parallel(n_jobs=num_cores)(delayed(parallelize)(j) for j in tqdm(range(0, len(neighborhood_sequence[i]))))
        cluster_scores = dict((index, score) for index, score in cluster_scores)
        all_scores.append([cluster_scores[j] for j in range(0, len(neighborhood_sequence[i]))])

    return all_scores

def cutoff_score(adata, end_point, neighborhood_sequence, all_scores, cut_off=0.0):
    """Deletes all cells with a score lower than the cutoff.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    map_state: matrix 
        Which projection to use for the data.
    end_point: `integer`
        End point cluster to which the trajectories belong.
    neighborhood_sequence: 
        List of the cells along the trajectory.
    all_scores:
        List of arrays with the score for each cell allong the trajectory
    Returns
    -------
    Returns list of arrays containing the scores of all cells at each step.
    """
    final_cluster = adata.uns['trajectories']["trajectories_coordinates"][end_point]
    directional_neighborhood_sequence = neighborhood_sequence
    if cut_off != None:
        for i in tqdm(range(len(final_cluster))):
            cutoff_score = cut_off
            for j in range(0, len(neighborhood_sequence[i])):
                 directional_neighborhood_sequence[i][j] = np.delete(directional_neighborhood_sequence[i][j], np.where(all_scores[i][j]<=cutoff_score))
                 all_scores[i][j] = np.delete(all_scores[i][j], np.where(all_scores[i][j]<=cutoff_score))
    return directional_neighborhood_sequence, all_scores
    
def time_step_cells(adata, directional_neighborhood_sequence, map_state):
    """Calculates the average time step for each cell.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    map_state: matrix 
        Which projection to use for the data.
    directional_neighborhood_sequence: 
        List of the cells along the trajectory after the cutoff
    
    Returns
    -------
    Returns list of vectors, where for each cell index the average step is denoted.
    """

    step_counter = np.zeros((len(directional_neighborhood_sequence), len(map_state)))
    occurance_counter = np.zeros((len(directional_neighborhood_sequence), len(map_state)))
    for i in range(len(directional_neighborhood_sequence)):
        for j in range(len(directional_neighborhood_sequence[i])):
            step_counter[i,directional_neighborhood_sequence[i][j]]+=j    
            occurance_counter[i,directional_neighborhood_sequence[i][j]]+=1  
    average_step = np.divide(step_counter, occurance_counter, out=np.zeros_like(step_counter), where=occurance_counter!=0)
    average_step[occurance_counter==0] = np.NaN
    return average_step
    
def cytopath(adata, basis="umap", neighbors_basis='pca', surrogate_cell=False, fill_cluster=True, n_neighbors_cluster=30, cluster_freq=0.25, n_neighbors='auto', cut_off=0.0, num_cores=1):
    """Calculates the average time step for each cell.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    basis: str (default: umap)
        The space in which the neighboring cells should be searched
    surrogate_cell: Boolean (default:False)
        Whether or not a surrogate cell should be used for the neighborhood search
    n_neighbors:  str/int/float (default:'auto')
        Number of neighbors to searched along the average trajectory.
    cut_off_percentile: float (default:0.0)
        Cuttof for the directionality score along the average trajectory for the cells.
    Returns
    -------
    Returns adata.uns['trajectories']["cells_along_trajectories"]: List of arrays, which denote the average step for cell.
    adata.uns['trajectories']["cells_along_trajectories_each_step"]: List of arrays containing the cell indexes for each step
    """
    cell_score=[]
    step_ordering_trajectory=[]
    cells_trajectory=[]
    map_state=adata.obsm['X_'+basis]
    for end_point_cluster in adata.uns['trajectories']["trajectories_coordinates"].keys():
        adata.uns['trajectories']['compositional_clusters'] = {}
        if surrogate_cell==True:
            print('Anchoring trajectories for end point  ' + end_point_cluster +' to cells in dataset and computing neighborhoods')
            neighborhood_sequence, cell_sequences = surrogate_cell_neighborhood_finder(adata, end_point_cluster, map_state, mode='distances', fill_cluster=fill_cluster, n_neighbors_cluster=n_neighbors_cluster,
                                                                                    cluster_freq=cluster_freq, n_neighbors=n_neighbors, recurse_neighbors=True)
        else:
            print('Computing neighborhoods of trajectories for end point ' + end_point_cluster +' at each step')
            neighborhood_sequence = cell_neighborhood_finder(adata, map_state, end_point=end_point_cluster, neighbors_basis=neighbors_basis, fill_cluster=fill_cluster, n_neighbors_cluster=n_neighbors_cluster, cluster_freq=cluster_freq,
                                                             n_neighbors=n_neighbors)

        print('Computing alignment score of cells in trajectory neighborhood w.r.t. trajectories for end point  ' + end_point_cluster)
        all_scores=directionality_score(adata, neighborhood_sequence=neighborhood_sequence, end_point=end_point_cluster, map_state=map_state, num_cores=num_cores)

        print('Removing cells below cutoff threshold from trajectories for end point  ' + end_point_cluster +' (i.e. cells neighborhood)')
        directional_neighborhood_sequence, all_scores = cutoff_score(adata, end_point=end_point_cluster, neighborhood_sequence=neighborhood_sequence, 
                                                                  all_scores=all_scores, cut_off=cut_off)
        cell_score.append(all_scores)

        print('Computing step based pseudotime in trajectories for end point  ' + end_point_cluster +' (i.e. cells neighborhood)')
        step_ordering_trajectory.append(time_step_cells(adata, directional_neighborhood_sequence=directional_neighborhood_sequence, map_state=map_state))
        cells_trajectory.append(directional_neighborhood_sequence)
    
    cells_along_trajectories = []
    end_point_clusters = list(adata.uns['trajectories']["trajectories_coordinates"].keys())
    for i in range(len(cells_trajectory)):
        for j in range(len(cells_trajectory[i])):
            for k in range(len(cells_trajectory[i][j])):
                for l in range(len(cells_trajectory[i][j][k])):
                    cells_along_trajectories.append((end_point_clusters[i], j, k, cells_trajectory[i][j][k][l], cell_score[i][j][k][l]))
       
    adata.uns['trajectories']["cells_along_trajectories"] = np.concatenate(step_ordering_trajectory)
    adata.uns['trajectories']["cells_along_trajectories_each_step"] = np.rec.fromrecords(cells_along_trajectories, 
                                                                                         dtype=[('End point', 'U32'), ('Trajectory', int), 
                                                                                                ('Step', int), ('Cell', int), ('Allignment Score', float)])

