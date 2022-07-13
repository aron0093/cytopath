import numpy as np
from scipy import stats, spatial, sparse

#from dtw import dtw
from fastdtw import fastdtw
from hausdorff import hausdorff_distance

import networkit
from networkit.community import ParallelLeiden, PLM

from sklearn.cluster import AffinityPropagation, KMeans

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

def get_graph(dm):
    """
    Builds a networkit graph from the input.
    :param dm:
    :return: networkit graph
    """

    m, _ = dm.shape
    g = networkit.Graph(m, weighted=True)
    mask_x, mask_y = np.mask_indices(m, np.tril, -1)
    masking_zip = zip(mask_x, mask_y, dm[mask_x, mask_y])

    for nodeA, nodeB, weight in masking_zip:
        if weight == 0:
            continue
        g.addEdge(nodeA, nodeB, weight)

    return g

# First of two stage clustering of markov chains 
def preclustering(adata, all_seq_cluster, sequence_coordinates, basis="umap", distance='euclidean'):

   if basis is None:
       map_state = adata.layers['Ms']
   elif type(basis) == list:
       map_state = adata[:, basis].layers['Ms']
   else:
       map_state = adata.obsm['X_'+basis]
        
   all_chains=sequence_coordinates

   # Calculate hausdorff distance
   print('Calculating hausdorff distances')
   distances=np.zeros((len(all_chains), len(all_chains)))
   for i in tqdm(range(len(distances))):
        for j in  range(len(distances)-i):
            haus = hausdorff_distance(all_chains[j+i], all_chains[i], distance=distance)
            distances[i,j+i] = haus
            distances[j+i,i] = haus   
   
   # Perform clustering using hausdorff distance
   print('Clustering using hausdorff distances')

   #cluster_labels = AffinityPropagation(affinity="precomputed").fit_predict(affinity)
   graph = get_graph(distances)
   cluster = networkit.community.detectCommunities(graph, algo=ParallelLeiden(graph, gamma=1.25))
   cluster.compact()
   cluster_labels = np.array(cluster.getVector())
   clusters = list(cluster.getSubsetIds())

    # Pairwise alignment of chains within stage 1 clusters using DTW

   cluster_chains=[]
   cluster_strength=[]

   print('Forming trajectory by aligning clusters') 
   # Pairwise allignment until only 1 avg. trajectory remains per cluster
   for i in tqdm(range(len(clusters))):

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
                    #dtw_obj = dtw(average_trajectory[pairs[l][0]],average_trajectory[pairs[l][1]])
                    #alligned_trajectory = np.stack((dtw_obj.index1, dtw_obj.index2)).T

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

# Second clustering of aligned stage 1 clusters that output trajectories
def clustering(adata, sequence_coordinates, cluster_chains, cluster_strength, cluster_labels_1, 
               n_clusters=None, method=None, distance='euclidean', smoothing=False):
        
    average_cl_d=np.zeros((len(cluster_chains), len(cluster_chains)))

    for i in tqdm(range(len(average_cl_d))):
        for j in  range(len(average_cl_d)):
            average_cl_d[i,j] = hausdorff_distance(cluster_chains[i], cluster_chains[j], distance=distance)
    
    print('Clustering using hausdorff distances')
    if n_clusters is None:
        graph = get_graph(average_cl_d)
        cluster = networkit.community.detectCommunities(graph, algo=ParallelLeiden(graph, gamma=0.75))
        cluster.compact()
        cluster_labels = np.array(cluster.getVector())
        clusters = list(cluster.getSubsetIds())
    else:
        if method=='kmeans':
            cluster_labels = KMeans(n_clusters=n_clusters, precompute_distances=True).fit_predict(average_cl_d_affinity)
            clusters=np.unique(cluster_labels)
        
    all_trajectories_labels = np.zeros((len(cluster_labels_1)))
    
    for i in range(len(cluster_labels)):
        all_trajectories_labels[np.where(cluster_labels_1==i)]=cluster_labels[i]
    
    cluster_strength=np.asarray(cluster_strength,dtype=int)
    
    final_cluster_strength= np.empty([len(clusters)],dtype=int)
    for i in range(len(clusters)):
        final_cluster_strength[i]=len(np.where(all_trajectories_labels==i)[0].astype(int))
    
    print('Forming trajectory by aligning clusters') 
    
    # pairwise allignment until only 1 avg. trajectory remains per cluster
    final_cluster=[]
    all_chains=sequence_coordinates
    for i in tqdm(range(len(clusters))):
        index_cluster=np.where(all_trajectories_labels==clusters[i])
        aver_tr_cluster=all_chains[index_cluster]
        if len(aver_tr_cluster)>1:
            average_trajectory=aver_tr_cluster             
            while len(average_trajectory)>1:
                pair_range=range(len(average_trajectory))
                pairs=[]
                for j in range(int((len(average_trajectory)-(len(average_trajectory) % 2))/2)):
                    pairs.append([pair_range[j*2],pair_range[((j)*2)+1]])
                if (len(average_trajectory) % 2) != 0:
                    pairs.append([pair_range[-2],pair_range[-1]])
                
                average_traj_while=[]
                for l in range(len(pairs)):
                    #dtw_obj = dtw(average_trajectory[pairs[l][0]],average_trajectory[pairs[l][1]])
                    #alligned_trajectory = np.stack((dtw_obj.index1, dtw_obj.index2)).T

                    alligned_trajectory=fastdtw(average_trajectory[pairs[l][0]],average_trajectory[pairs[l][1]])[1]
                    alligned_trajectory=np.asarray(alligned_trajectory)
                    alligned_tr =np.zeros((2,len(average_trajectory[0][0])))
                    alligned_av =np.zeros((len(alligned_trajectory),len(average_trajectory[0][0])))
                    for n in range(len(alligned_trajectory)):
                        alligned_tr[0,:]=average_trajectory[pairs[l][0]][alligned_trajectory[n,0]]
                        alligned_tr[1,:]=average_trajectory[pairs[l][1]][alligned_trajectory[n,1]]
                        alligned_av[n,:]=np.mean(alligned_tr,axis=0)
                    average_traj_while.append(alligned_av)
                average_trajectory=average_traj_while
            average_alligned_tr=average_trajectory[0]
        else: 
            average_alligned_tr=aver_tr_cluster[0]
            
        final_cluster.append(average_alligned_tr)
           
    # TODO: Moving average concept needs development

    if smoothing:
        for i in range(len(final_cluster)):
            window_width=4
            cumsum_vec = np.cumsum(np.insert(final_cluster[i][:,0], 0, 0)) 
            ma_vec_Y = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width    
            cumsum_vec = np.cumsum(np.insert(final_cluster[i][:,1], 0, 0)) 
            ma_vec_X = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            final_cluster[i]=np.zeros((len(ma_vec_X),2))
            final_cluster[i][:,0]=ma_vec_Y
            final_cluster[i][:,1]=ma_vec_X

    return final_cluster, final_cluster_strength
    
def sample_clustering(adata, basis="umap", smoothing=False, cluster_num=None, method='kmeans', distance='euclidean'):
    """Clusters samples for each terminal region and estimates trajectories.
    
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    basis: str/list (default: umap)
        The space in which the neighboring cells should be searched. If None imputed expression is used. If list, imputed expression from supplied genes is used.
    smoothing: Boolean (default:False)
        Whether or not to smooth over the trajectories.
    cluster_num:  list (default:None)
        Number of trajectories (clusters) to be expected for each terminal region.
    method: str (default:'kmeans'):
        Which clustering method to use.
    distance: str (default:'cosine'):
        Which distabce metric to use (see py-hausdorff for details)
    Returns
    -------
    adata.uns["trajectories_2"]: Dictionary of the average trajectories
    adata.uns["trajectories"]: Long format of the average trajectories
    adata.uns["trajectories_count"]: Number of trajectories contributing to each cluster
    """
    
    # Check if samples has been run and information is complete
    try:
        adata.uns['run_info']['didrun']
    except:
        raise ValueError("Run cytopath.sampling before trajectory inference.")

    if basis is None:
        adata.uns['run_info']['trajectory_basis'] = str(adata.X.shape[1]) + 'D_expression'
    elif type(basis) == list:
        adata.uns['run_info']['trajectory_basis'] = str(len(basis)) + 'D_custom_geneset_expression'
        adata.uns['run_info']['trajectory_basis_geneset'] = basis
    else:
        adata.uns['run_info']['trajectory_basis'] = str(adata.obsm['X_'+basis].shape[1]) + 'D_'+ basis
        
    end_clusters, sel_end_clusters, end_points = end_point_cluster(adata)

    adata.uns['run_info']['end_point_clusters'] = end_clusters

    final_trajectories = []
    final_cluster_count = []

    subtrajectory_dict={}
    subtrajectory_dict['trajectory_samples'] = {}
    subtrajectory_dict['subtrajectory_labels'] = {}
    subtrajectory_dict['subtrajectory_coordinates'] = {}

    for i in range(len(end_clusters)):
        trajectories = end_point_trajectories(adata, end_clusters, sel_end_clusters, end_points, cluster=i)
        if trajectories.shape[0] > 0:
            sequence_coordinates = coordinate_assigner(adata, trajectories, basis=basis)  
            subtrajectory_dict['trajectory_samples'][end_clusters[i]] = trajectories

        # TODO: Allow different number of trajectories per terminal region
            if cluster_num is not None:
                print("Clustering and aligning samples for end point " + str(end_clusters[i]))
                cluster_chains, cluster_strength, cluster_labels = preclustering(adata, sequence_coordinates=sequence_coordinates,
                                                                                         basis=basis, all_seq_cluster=trajectories, 
                                                                                         distance=distance)
                
                print("Final clustering done. Aligning clusters for end point " + str(end_clusters[i]))
                final_trajectory, final_cluster_strength = clustering(adata, sequence_coordinates, cluster_chains, cluster_strength, 
                                                                            cluster_labels_1=cluster_labels, n_clusters=cluster_num, method=method, 
                                                                            distance=distance, smoothing=smoothing)
            elif cluster_num is None:

                print("Stage 1 clustering done. Aligning clusters for end point " + str(end_clusters[i]))
                cluster_chains, cluster_strength, cluster_labels = preclustering(adata, sequence_coordinates=sequence_coordinates,
                                                                                 basis=basis, all_seq_cluster=trajectories, 
                                                                                 distance=distance)
                
                
                print("Final clustering done. Aligning clusters for end point " + str(end_clusters[i]))
                final_trajectory, final_cluster_strength = clustering(adata, sequence_coordinates, cluster_chains, cluster_strength, 
                                                                            cluster_labels_1=cluster_labels, n_clusters=None, method=None,
                                                                            distance=distance, smoothing=smoothing)
                #final_trajectory = cluster_chains
                #final_cluster_strength = cluster_strength
                

            subtrajectory_dict['subtrajectory_labels'][end_clusters[i]] = cluster_labels
            cluster_coordinates = []
            for k in range(len(cluster_chains)):
                subtraj_coords = [np.array(cluster_chains[k]).T[j] for j in range(np.array(cluster_chains[k]).shape[1])]
                dtypes_ = [('Coordinate_'+str(j), float) for j in range(np.array(cluster_chains[k]).shape[1])]
                
                subtraj_coords.extend([np.arange(np.array(cluster_chains[k]).shape[0]), np.array([cluster_labels[k]]*np.array(cluster_chains[k]).shape[0])])
                dtypes_.extend([('Step', int), ('Subtrajectory', int)])
                
                subtraj_coords = np.rec.fromarrays(subtraj_coords, dtype=dtypes_)
                cluster_coordinates.append(subtraj_coords)
      
            subtrajectory_dict['subtrajectory_coordinates'][end_clusters[i]] = np.hstack(cluster_coordinates)

            # Discard trajectories that contain less than 10 % of all samples
            indexes = np.where(final_cluster_strength>(0.1*np.sum(final_cluster_strength)))[0]

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
                

    adata.uns["subtrajectories"] = subtrajectory_dict
    adata.uns['trajectories'] = {}
    adata.uns['trajectories']["trajectories_coordinates"] = trajectory_dict
    adata.uns['run_info']["trajectories_sample_counts"] = trajectory_sample_count_dict


    