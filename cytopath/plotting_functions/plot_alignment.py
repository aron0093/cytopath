import numpy as np

from scipy import spatial

import matplotlib.pyplot as plt
      
def fft_smoothing(coords):

    #TODO: More relevant procedure required
    signal = coords[:,0] + 1j*coords[:,1]

    # FFT and frequencies
    fft = np.fft.fft(signal)
    freq = np.fft.fftfreq(signal.shape[-1])

    # filter
    cutoff = 0.1
    fft[np.abs(freq) > cutoff] = 0

    # IFFT
    signal_filt = np.fft.ifft(fft)

    coords[:,0] = signal_filt.real
    coords[:,1] = signal_filt.imag

    return coords

def pl_cytopath_alignment(adata, basis="umap", smoothing=False, figsize=(15,4), size = 3, 
                          show=True, save=False,save_type='png', folder=""):
    
    map_state = adata.obsm['X_'+basis]
    av_allign_score_glob=[]
    std_allign_score_glob=[]
    
    step_time = adata.uns['trajectories']['step_time']
    fate_prob = adata.uns['trajectories']['cell_fate_probability']

    sequence=0

    # TODO: Separate per step average alignment score calculation from plotting
    for end_point_cluster in adata.uns['run_info']["end_point_clusters"]:
        trajectories = adata.uns['trajectories']["cells_along_trajectories_each_step"]\
                       [np.where(adata.uns['trajectories']["cells_along_trajectories_each_step"]["End point"]==end_point_cluster)[0]]
        for i in range(adata.uns['run_info']['trajectory_count'][end_point_cluster]):

            av_trajectories=trajectories[np.where(trajectories["Trajectory"]==i)[0]]
            av_allign_score=np.zeros((len(np.unique(av_trajectories["Step"]))))
            std_allign_score=np.zeros((len(np.unique(av_trajectories["Step"]))))
            for l in range(len(np.unique(av_trajectories["Step"]))):
                av_allign_score[l]=np.average((av_trajectories[np.where(av_trajectories["Step"]==l)[0]]["Allignment Score"]))
                std_allign_score[l]=np.std((av_trajectories[np.where(av_trajectories["Step"]==l)[0]]["Allignment Score"]))
            
            # Plotting
            path = folder+"_end_point_"+end_point_cluster+"_cytopath_"+str(i)+\
                 "occurance"+str(adata.uns['run_info']["trajectories_sample_counts"][end_point_cluster][i])+"."+save_type

            fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=figsize) 
            ax1.plot(range(len(np.unique(av_trajectories["Step"]))), av_allign_score, color='black')
            ax1.fill_between(range(len(np.unique(av_trajectories["Step"]))), 
                             av_allign_score+std_allign_score, av_allign_score-std_allign_score, facecolor='grey', alpha=0.6)
            ax1.set_ylabel('Mean/std. of alignment scores per step')
            ax1.set_xlabel('Steps')
            
            # Plot step size for aligned cells
            sc_step = ax2.scatter(map_state[:,0], map_state[:,1], alpha=0.6, s=size, color="whitesmoke")
            sc_step = ax2.scatter(map_state[:,0], map_state[:,1], alpha=0.9, s=size, 
                                  vmin=0, vmax=np.nanmax(step_time), c=step_time[sequence,:], cmap='YlGnBu')
            fig.colorbar(sc_step, ax=ax2, label='Step time')

            ax2.set_ylabel(basis.upper()+' 2')
            ax2.set_xlabel(basis.upper()+' 1')

            ax2.set_title('End point: {} {} Support: {}'.format(end_point_cluster, i, 
                                                                adata.uns['run_info']['trajectories_sample_counts'][end_point_cluster][i]))

            # Plot alignment score
            sc_score = ax3.scatter(map_state[:,0], map_state[:,1], alpha=0.6, s=size, color="whitesmoke")
            sc_score = ax3.scatter(map_state[:,0], map_state[:,1], alpha=0.9, s=size, 
                                   vmin=0, vmax=1, c=fate_prob[sequence,:], cmap='Reds')
            fig.colorbar(sc_score, ax=ax3, label='Cell fate probability')

            ax3.set_ylabel(basis.upper()+' 2')
            ax3.set_xlabel(basis.upper()+' 1')

            # Plot trajectory
            if basis in adata.uns['run_info']['trajectory_basis']:
                coords = np.array(adata.uns['trajectories']['trajectories_coordinates'][end_point_cluster]['trajectory_'+str(i)+'_coordinates'])

            elif ('pca' in adata.uns['run_info']['trajectory_basis']) and (basis != 'pca'):
                coords_ = np.array(adata.uns['trajectories']['trajectories_coordinates'][end_point_cluster]['trajectory_'+str(i)+'_coordinates'])

                cell_sequences=[]
                for j in range(len(coords_)):
                    cell_sequences.append(spatial.KDTree(adata.obsm['X_pca']).query(coords_[j])[1])
                
                coords = map_state[cell_sequences]

            if smoothing == True:
                coords = fft_smoothing(coords)
                
            ax2.plot(coords[:, 0], coords[:, 1], color='black')
            ax3.plot(coords[:, 0], coords[:, 1], color='black')

            plt.tight_layout()

            if save:
                fig.savefig(path, bbox_inches='tight', dpi=300)
            if show:
                plt.show()
            # End plotting

            sequence+=1

            av_allign_score_glob.append(av_allign_score)                                
            std_allign_score_glob.append(std_allign_score)

            