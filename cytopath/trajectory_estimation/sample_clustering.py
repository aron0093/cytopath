import numpy as np

from fastdtw import fastdtw
from hausdorff import hausdorff_distance

from sklearn.cluster import HDBSCAN

from tqdm.auto import tqdm

# Function to define terminal regions for trajectories
def end_point_cluster(adata):

    cluster_key = adata.uns['run_info']['cluster_key']
    sel_end_points = adata.uns['run_info']['end_points']
    sel_end_clusters_list = adata.obs[cluster_key][sel_end_points]
    sel_end_clusters_unique = np.unique(sel_end_clusters_list)

    return sel_end_clusters_unique, sel_end_clusters_list, sel_end_points

# Function to retrieve markov chains terminating at a terminal region  
def end_point_trajectories(adata, end_clusters, sel_end_clusters, sel_end_points, cluster):
    
    all_seq = adata.uns["samples"]["cell_sequences"]
    end_cluster_indexes = np.concatenate(np.where(sel_end_clusters==end_clusters[cluster]), axis=0)
    # Retrieve only sequences, which end in the end points of the cluster
    indexes_cluster=[]
    end_points=np.asarray([item[-1] for item in all_seq])
    
    for l in range(len(end_cluster_indexes)):
        indexes_cluster.append(np.where(end_points==sel_end_points[end_cluster_indexes[l]]))
    indexes_cluster = np.concatenate(indexes_cluster, axis=1)
    all_seq_cluster=all_seq[indexes_cluster]
    all_seq_cluster=all_seq_cluster[0]

    return all_seq_cluster

# Function to assign coordiantes to markov chains in PCA space or another embedding
def coordinate_assigner(adata, all_seq_cluster, basis="umap"):

    if basis is None:
        map_state = adata.layers['Ms']
    elif type(basis) == list:
        map_state = adata[:, basis].layers['Ms']
    else:
        map_state = adata.obsm['X_'+basis]
        
    all_chains = []
    for r in range(len(all_seq_cluster)):
        x_y_chain = []
        for j in range(len(all_seq_cluster[r])):
            x_y_chain.append(map_state[int(all_seq_cluster[r][j]),:])
        all_chains.append(np.vstack(x_y_chain))
    all_chains=np.stack(all_chains, axis=0)

    return all_chains

def distance_matrix(sequence_coordinates, distance='euclidean', num_cores=1):

    # Not neccessarily same as coords for DTW
    all_chains=sequence_coordinates

    # Calculate hausdorff distance
    print('Calculating hausdorff distances')
    distances=np.zeros((len(all_chains), len(all_chains)))
    for i in tqdm(range(len(distances))):
        for j in  range(len(distances)-i):
            haus = hausdorff_distance(all_chains[j+i], all_chains[i], distance=distance)
            distances[i,j+i] = haus
            distances[j+i,i] = haus   
    
    return distances

# First of two stage clustering of markov chains 
def clustering(sequence_coordinates, distances, num_cores=1):

    all_chains=sequence_coordinates

    # Perform clustering using hausdorff distance
    print('Clustering using hausdorff distances')
    hdb = HDBSCAN(min_cluster_size=max(2, int(len(all_chains)*0.01)), 
                  metric='precomputed', n_jobs=num_cores, 
                  allow_single_cluster=True)
    cluster_labels = hdb.fit_predict(distances)
    clusters = np.unique(cluster_labels)

    # Pairwise alignment of chains within clusters using DTW
    cluster_chains=[]
    cluster_strength=[]

    print('Forming trajectory by aligning clusters') 
    # Pairwise allignment until only 1 avg. trajectory remains per cluster
    # TODO: Replace with multivariate DTW alignment
    for i in tqdm(range(len(clusters))):

        if i>=0:

            index_cluster = np.where(cluster_labels==clusters[i])[0]
            aver_tr_cluster = all_chains[index_cluster]
            cluster_strength.append(len(index_cluster))

            if len(aver_tr_cluster)>1:
                average_trajectory = aver_tr_cluster

                while len(average_trajectory)>1:
                    pair_range = range(len(average_trajectory))
                    pairs=[]
                    for j in range(int((len(average_trajectory)-(len(average_trajectory) % 2))/2)):
                        pairs.append([pair_range[j*2], pair_range[((j)*2)+1]])
                    if (len(average_trajectory) % 2) != 0:
                        pairs.append([pair_range[-2], pair_range[-1]])
                    
                    average_traj_while=[]
                    for l in range(len(pairs)):

                        alligned_trajectory = fastdtw(average_trajectory[pairs[l][0]],average_trajectory[pairs[l][1]])[1]
                        alligned_trajectory = np.asarray(alligned_trajectory)
                        alligned_tr = np.zeros((2,len(average_trajectory[0][0])))
                        alligned_av = np.zeros((len(alligned_trajectory),len(average_trajectory[0][0])))

                        for n in range(len(alligned_trajectory)):
                            alligned_tr[0,:] = average_trajectory[pairs[l][0]][alligned_trajectory[n,0]]
                            alligned_tr[1,:] = average_trajectory[pairs[l][1]][alligned_trajectory[n,1]]
                            alligned_av[n,:] = np.mean(alligned_tr,axis=0)
                        average_traj_while.append(alligned_av)
                    average_trajectory=average_traj_while
                average_alligned_tr=average_trajectory[0]
            else: 
                average_alligned_tr=aver_tr_cluster[0]
                
            cluster_chains.append(average_alligned_tr)

    return cluster_chains, cluster_strength, cluster_labels
    
def sample_clustering(adata, neighbors_basis='pca', basis="umap", cluster_num=None, 
                      support=0.10, distance='euclidean', num_cores=1):
    """Clusters samples for each terminal region and estimates trajectories.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    neighbors_basis: str/list (default: pca)
        The space in which the distances and neighbors are computed. If None expression is used. 
        If list, imputed expression from supplied genes is used.
    basis: str (default: umap)
        The space in which DTW is performed. Use for projection.
    cluster_num:  list (default:None)
        Number of trajectories (clusters) to be expected for each terminal region.
    support: float(default:0.10)
        Minimum ratio of samples that must support a trajectory
    distance: str (default:'cosine'):
        Which distabce metric to use (see py-hausdorff for details)
    num_cores: 'integer' (default:1)
        Number of cpu cores to use.
    Returns
    -------
    adata.uns["trajectories"]: Creates and updates dictionary with Markov chains.
    """
    
    # Check if samples has been run and information is complete
    try:
        adata.uns['run_info']['didrun']
    except:
        raise ValueError("Run cytopath.sampling before trajectory inference.")

    if neighbors_basis is None:
        adata.uns['run_info']['trajectory_basis'] = str(adata.X.shape[1]) + 'D_expression'
    elif type(neighbors_basis) == list:
        adata.uns['run_info']['trajectory_basis'] = str(len(neighbors_basis)) + 'D_custom_geneset_expression'
        adata.uns['run_info']['trajectory_basis_geneset'] = neighbors_basis
    else:
        adata.uns['run_info']['trajectory_basis'] = str(adata.obsm['X_'+neighbors_basis].shape[1]) + 'D_'+ neighbors_basis

    adata.uns['run_info']['projection_basis'] = basis
        
    end_clusters, sel_end_clusters, end_points = end_point_cluster(adata)
    adata.uns['run_info']['end_point_clusters'] = end_clusters

    final_trajectories = []
    final_cluster_count = []
    for i in range(len(end_clusters)):
        trajectories = end_point_trajectories(adata, end_clusters, sel_end_clusters, end_points, cluster=i)
        if trajectories.shape[0] > 0:
            # Compute coordinates for DTW alignment
            sequence_coordinates = coordinate_assigner(adata, trajectories, basis=basis)

            # Compute distances between simulations for clustering
            sequence_coordinates_distance = coordinate_assigner(adata, trajectories, basis=neighbors_basis)
            distances=distance_matrix(sequence_coordinates_distance, distance=distance, num_cores=num_cores)

            if cluster_num is not None:
                # TODO: Allow different number of trajectories per terminal region
                raise NotImplementedError('Custom selection of number of trajectories has not been implemented.')

            elif cluster_num is None:

                final_trajectory, final_cluster_strength, cluster_labels = clustering(sequence_coordinates=sequence_coordinates,
                                                                                      distances=distances, num_cores=num_cores)
                print("Sample clustering done. Aligning clusters for end point " + str(end_clusters[i]))

            # Discard trajectories that contain less than 20 % of all samples (traj_num)
            traj_num=adata.uns['samples']['cell_sequences'].shape[0]/adata.uns['run_info']['end_point_clusters'].shape[0]

            indexes = np.where(np.array(final_cluster_strength)>(support*traj_num))[0]

            final_trajectories.append(np.array(final_trajectory)[indexes])
            final_cluster_count.append(np.array(final_cluster_strength)[indexes])   
        
    trajectory_dict = {}
    trajectory_sample_count_dict = {}
    adata.uns['run_info']['trajectory_count'] = {}

    for i in range(len(final_trajectories)):
        trajectory_dict[end_clusters[i]] = {}
        trajectory_sample_count_dict[end_clusters[i]] = []
        adata.uns['run_info']['trajectory_count'][end_clusters[i]] = len(final_trajectories[i])

        for j in range(len(final_trajectories[i])):
            trajectory_dict[end_clusters[i]]['trajectory_'+str(j)+'_coordinates'] = []
            trajectory_sample_count_dict[end_clusters[i]].append(final_cluster_count[i][j])

            for k in range(len(final_trajectories[i][j])):
                trajectory_dict[end_clusters[i]]['trajectory_'+str(j)+'_coordinates'].append(final_trajectories[i][j][k])
                
    adata.uns['trajectories'] = {}
    adata.uns['trajectories']["trajectories_coordinates"] = trajectory_dict
    adata.uns['run_info']["trajectories_sample_counts"] = trajectory_sample_count_dict


    