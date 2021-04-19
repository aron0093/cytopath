import math
import numpy as np
from scipy import sparse
from tqdm import tqdm
from joblib import Parallel, delayed
import warnings

def matrix_prep(trans_matrix):

    # For increased speed: Create and save np array of the cumsum on nonzero entries-->80% faster
    trans_matrix_ = sparse.csr_matrix(trans_matrix)
    trans_matrix_indices = np.split(trans_matrix_.indices, trans_matrix_.indptr)[1:-1]
    trans_matrix_probabilites = np.split(trans_matrix_.data, trans_matrix_.indptr)[1:-1]
    
    for i in range(len(trans_matrix_probabilites)):
        trans_matrix_probabilites[i]=np.cumsum(trans_matrix_probabilites[i])
    
    # For increased speed: Save data and indices of nonzero entries in np array--> 210% faster
    trans_matrix_indice_zeros = np.zeros((len(trans_matrix_probabilites), max(np.count_nonzero(trans_matrix, axis=1))), dtype=np.int32)
    trans_matrix_probabilites_zeros = np.zeros((len(trans_matrix_probabilites), max(np.count_nonzero(trans_matrix, axis=1))), dtype=np.float32)
    
    for i in range(len(trans_matrix_probabilites)):
        for j in range(len(trans_matrix_probabilites[i])):
            trans_matrix_indice_zeros[i, j] = trans_matrix_indices[i][j]
            trans_matrix_probabilites_zeros[i, j] = trans_matrix_probabilites[i][j]
    
    # FIXME: Ensure that cummulative sum reaches 1 precisely (rounding to 3 decimals ensures that this almost never happens)
    return trans_matrix_indice_zeros, np.round(trans_matrix_probabilites_zeros, decimals=3)

# Markov sim function: sample according to the probabilities in the transition matrix.
# Gets used in the parallel function.
def markov_sim(j, sim_number, max_steps, root_cells, clusters, trans_matrix, trans_matrix_indice_zeros, trans_matrix_probabilites_zeros):
    # Create empty arrays for the simulation and sample random numbers for all samples and max steps
    prob_state = np.empty((sim_number, max_steps-1), dtype=np.float32)
    currState = np.empty((sim_number, max_steps), dtype=np.int32)
    clust_state = np.empty((sim_number, max_steps), dtype='U8')

    currState[:, 0] = root_cells[j]
    clust_state[:, 0] = clusters[root_cells[j]]
    random_step = np.random.rand(sim_number, max_steps-1)
    
    for i in range(1, max_steps):
        for k in range(sim_number):
            # Retrieve next cell transition and save the probabilities of the transition and the cluster the cell is in
            currState[k, i] = trans_matrix_indice_zeros[currState[k, i-1], 
                                                        np.where(random_step[k, i-1] <= trans_matrix_probabilites_zeros[currState[k,i-1], :])[0][0]]
            prob_state[k, i-1] = trans_matrix[currState[k, i-1], currState[k, i]]
            clust_state[k, i] = clusters[currState[k, i]]
    return currState, np.sum(-np.log(prob_state), axis=1), clust_state

def sampling(data, matrix_key = 'T_forward', cluster_key = 'louvain', max_steps=200, traj_number=1000, sim_number=2000,
                end_point_probability=0.95, root_cell_probability=0.95, end_points=[], root_cells=[], end_clusters = [], root_clusters=[], min_clusters=3,
                normalize=False, unique=True, num_cores=1, copy=False):
    
    """Markov sampling of cell sequences starting from defined root cells to defined terminal regions based on a cell-cell transition probability matrix.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with transition probabilities, end points and clustering.
    matrix_key: Key for accessing transition matrix stored in adata.uns
    cluster_key: Clustering of cells to use for terminal region annotation.  
    max_steps: `integer` (default: 200)
        Maximum number steps of simulation.
    traj_number: `integer`(default: 1000)
        Minimum number of samples generated for each terminal region.
    sim_number: `integer`(default: 2000)
        Initial number of samples generated for each terminal region.
    end_point_probability: 'float' (default:0.95)
        End point probability threshold at which a cell is considered to be a terminal state for samples.
    root_state_probability: 'float' (default:0.95)
        Start state probability threshold at which a cell is considered to be an origin for samples.
    end_points: `list` (default: [])
        Numerical indices of the cells considered to be end points for manual annotation. Precendence over end_point_probability and end_clusters.
    root_cells: `list` (default: [])
        Numerical indices of the cells considered to be root states for manual annotation. Precendence over root_state_probability and root_clusters.
    end_clusters: `list` (default: [])
        Cluster IDs to be considered end points. Precendence over end_point_probability.
    root_clusters: `list` (default: [])
        Cluster IDs to be considered root states. Precendence over root_state_probability.
    min_clusters: `integer` (default:3)
        Minium number of clusters covered by each simulation. (cluster_key)
    normalize: 'Boolean' (default: False)
        Toggle row sum normalization.
    unique: 'Boolean' (default:True)
        Only keep unique samples.
    num_cores: 'integer' (default:1)
        Number of cpu cores to use.
    copy: 'Boolean' (default: False)
        Create a copy of the anndata object or work inplace.
    
    Returns
    -------
    adata.uns["samples"]["cell_sequences"]: List of arrays containing the index of the cells in a sequence.
    adata.uns["samples"]["transition_score"]: List with all of the transition probabilities.
    adata.uns["samples"]["cluster_sequences"]: The sequence of clusters to which the cells belong.
    """

    adata = data.copy() if copy else data

    ### Run checks

    # Check if transition matrix exists
    try:
        adata.uns[matrix_key]
    except NameError:
        raise ValueError('Transition matrix must be provided in adata.uns[matrix_key].')
    
    # Check if a square matrix is provided
    try:
        assert len(adata.uns[matrix_key].shape) == 2
        assert adata.uns[matrix_key].shape[0] == adata.uns[matrix_key].shape[1]
        assert adata.uns[matrix_key].shape[0] == adata.shape[0]
    except:
        raise ValueError('Transition matrix must be of shape num_cells x num_cells')

    # Convert sparse matrx to numpy array    
    if sparse.issparse(adata.uns[matrix_key]):
        trans_matrix = np.float64(np.asarray(adata.uns[matrix_key].toarray()))
    else:
        trans_matrix = np.float64(np.asarray(adata.uns[matrix_key]))
    
    if sparse.issparse(adata.X):
        adata.X = np.asarray(adata.X.toarray())
    
    # Check if matrix has been normalized
    if normalize == False:
        row_sums = np.sum(trans_matrix, axis=1)
        try: 
            assert max(row_sums) <= 1.0 + 1e-3
            assert min(row_sums) >= 1.0 - 1e-3
        except AssertionError:
            raise ValueError('Transition matrix has not been normalized. Either pass normalize==True or provide a row sum normalized matrix.')

    # Row sum normalization 
    if normalize == True:
        new_array = trans_matrix
        normalization_v = new_array.sum(axis=1)
        for i in range(len(normalization_v)):
            trans_matrix[i,:]=trans_matrix[i,:]/normalization_v[i]
        
    # Check if cluster labels are provided
    try:
        adata.obs[cluster_key]
        assert adata.obs[cluster_key].dtype == 'category'
    except:
        raise ValueError('Categorical cluster labels must be provided in adata.obs[cluster_key]')
    
    # Use louvain warning
    if cluster_key != 'louvain':
        warnings.warn('A finely resolved clustering is required to identify terminal regions. Consider using louvain clustering.')

    ### Checks Done

    # Define root and end points, clusters
    clusters = adata.obs[cluster_key].astype(str).values
    adata.uns['run_info'] = {'cluster_key': cluster_key}

    trans_matrix_indice_zeros, trans_matrix_probabilites_zeros = matrix_prep(trans_matrix)

    # Retrieve all root cells
    if len(root_cells) == 0:
        if len(root_clusters) == 0:
            try:
                root_cells = np.asarray(np.where(np.asarray(adata.obs["root_cells"]) >= root_cell_probability))[0]
            except:
                raise ValueError('Root cells must be provided in adata.obs["root_cells"]. Run cytopath.terminal_states')
        else:
            root_cells = np.asarray(np.where(np.asarray(adata.obs[cluster_key].isin(root_clusters))))[0]
    adata.uns['run_info']['root_cells'] = root_cells

    # Retrieve all end point cells
    if len(end_points) == 0:
        if len(end_clusters) == 0:
            try:
                end_points = np.asarray(np.where(np.asarray(adata.obs["end_points"]) >= end_point_probability))[0]
            except:
                raise ValueError('End points must be provided in adata.obs["end_points"]. Run cytopath.terminal_states')
        else:
            end_points = np.asarray(np.where(np.asarray(adata.obs[cluster_key].isin(end_clusters))))[0]
    adata.uns['run_info']['end_points'] = end_points
                
    end_point_clusters = adata.obs[cluster_key][end_points].astype(str).values
    end_clusters_ = np.unique(end_point_clusters)
    
    # Initialize all empty lists
    traj_num = [0]
    glob_all_seq = []
    glob_prob_all_seq = []
    glob_cluster_seq_all = []
    sim_number_ = sim_number

    count = 1
    # Resample samples until the number of trajectories for each end point has been reached.
    while min(traj_num) < traj_number:

        print("Sampling round " + str(count))
        count+=1

        # TODO: Add transition probability matrix update that removes transitions exclusive to end points sampled adequately 
              # to increase sampling efficiency for larger datasets.

        # Parallelisation over the cells: 1 Thread for the samples of each cell.
        if len(root_cells) == 1:
            all_seq, prob_all_seq, cluster_seq_all = markov_sim(0, sim_number_, max_steps, root_cells, clusters, trans_matrix,
                                                                trans_matrix_indice_zeros, trans_matrix_probabilites_zeros)
        else:
            results = Parallel(n_jobs=num_cores)(delayed(markov_sim)(j, sim_number_, max_steps, root_cells, clusters, trans_matrix,
                                                                                trans_matrix_indice_zeros, 
                                                                                trans_matrix_probabilites_zeros) for j in tqdm(range(len(root_cells))))

            # Results are stored in np.arrays
            all_seq = np.empty((len(root_cells)*sim_number_, max_steps), dtype=np.int32)
            prob_all_seq = np.empty((len(root_cells)*sim_number_), dtype=np.float32)
            cluster_seq_all = np.empty((len(root_cells)*sim_number_, max_steps), dtype='U32')

            for i in range(len(root_cells)):
                all_seq[(i*sim_number_):((i*sim_number_)+sim_number_), :] = results[i][0]
                prob_all_seq[(i*sim_number_):((i*sim_number_)+sim_number_)] = results[i][1]
                cluster_seq_all[(i*sim_number_):((i*sim_number_)+sim_number_), :] = results[i][2]
            
            del results
               
        # Filter samples that visit less than min_clusters clusters
        indices = [idx for idx in np.arange(all_seq.shape[0]) if np.unique(cluster_seq_all[idx]).shape[0] >= min_clusters]
        if len(indices) == 0:
            raise ValueError('Zero samples obtained. Try increasing max_steps.')
        all_seq = [all_seq[i] for i in sorted(indices)]
        prob_all_seq = [prob_all_seq[i] for i in sorted(indices)]
        cluster_seq_all = [cluster_seq_all[i] for i in sorted(indices)]

        # If unique, then only keep unique samples
        if unique == True:
            indices = np.unique(all_seq, axis=0, return_index=True)[1]
            all_seq = [all_seq[i] for i in sorted(indices)]
            prob_all_seq = [prob_all_seq[i] for i in sorted(indices)]
            cluster_seq_all = [cluster_seq_all[i] for i in sorted(indices)]
        
        glob_all_seq.append(all_seq)
        glob_prob_all_seq.append(prob_all_seq)
        glob_cluster_seq_all.append(cluster_seq_all)
        
        # Count the number of samples to each end point
        all_seq = np.vstack(glob_all_seq)
        traj_num=[]
        for i in range(len(end_clusters_)):
            end_cluster_indices = np.concatenate(np.where(end_point_clusters == end_clusters_[i]), axis=0)
            # Retrieve only sequences, which end in the end points of the cluster
            indices_cluster=[]
            sample_end_points = np.asarray([item[-1] for item in all_seq])
            
            for l in range(len(end_cluster_indices)):
                indices_cluster.append(np.where(sample_end_points == end_points[end_cluster_indices[l]]))

            indices_cluster = np.concatenate(indices_cluster, axis=1)
            all_seq_cluster = all_seq[indices_cluster]
            all_seq_cluster = all_seq_cluster[0]
            traj_num.append(len(all_seq_cluster))
            
        sim_number_ = math.ceil(traj_number/(min(traj_num) + 1))*sim_number
    
    # Retrieve traj_number of samples for each end point.
    glob_prob_all_seq = np.concatenate(glob_prob_all_seq)
    all_seq=np.vstack(glob_all_seq)
    cell_sequences = []
    transition_score = []
    cluster_sequences = []

    for i in range(len(end_clusters_)):
        end_cluster_indices = np.concatenate(np.where(end_point_clusters == end_clusters_[i]), axis=0)
        indices_cluster = []
        sample_end_points = np.asarray([item[-1] for item in all_seq])
        
        for l in range(len(end_cluster_indices)):
            indices_cluster.append(np.where(sample_end_points == end_points[end_cluster_indices[l]]))
        indices_cluster = np.concatenate(indices_cluster, axis=1)
        all_seq_cluster = all_seq[indices_cluster]
        all_seq_cluster = all_seq_cluster[0]
        idx_range = np.random.randint(all_seq_cluster.shape[0], size=traj_number)
        
        cell_sequences.append(all_seq_cluster[idx_range])
        
        loc_tran_score = glob_prob_all_seq[indices_cluster[0]]
        transition_score.append(loc_tran_score[idx_range])

        loc_cluster_seq = np.vstack(glob_cluster_seq_all)[indices_cluster[0]]
        cluster_sequences.append(loc_cluster_seq[idx_range])
    
    # Save results
    adata.uns['run_info']['didrun'] = True

    adata.uns['samples'] = {}
    adata.uns['samples']["cell_sequences"] = np.vstack(cell_sequences)
    adata.uns['samples']["transition_score"] = np.concatenate(transition_score)
    adata.uns['samples']["cluster_sequences"]=np.vstack(cluster_sequences)

    return adata if copy else None