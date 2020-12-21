import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

def cytopath_merger(adata, overlap=0.5, num_cores=1):
    """
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with end points.
    overlap: float (default: 0.5)
        overlap until new bucket of cells is created.
    Returns
    -------
    adata.uns['trajectories']["cells_along_trajectories_each_step_merged"]: List of arrays containing the cell indexes and allignment scores for each step

    """
    cells_along_trajectories = adata.uns['trajectories']["cells_along_trajectories_each_step"].copy()
    cells_along_trajectories_merged = []
    
    # For each terminal region and each trajectory
    for end_point in adata.uns['run_info']['end_point_clusters']:
        print('Merging steps for trajectories for ' + end_point)
        sel_end_point_data = cells_along_trajectories[np.where(cells_along_trajectories["End point"]==end_point)[0]]
        for i in tqdm(range(adata.uns['run_info']["trajectory_count"][end_point])):
            # Find all trajectories
            sel_end_point_trajectories=sel_end_point_data[np.where(sel_end_point_data["Trajectory"]==i)[0]]
            start=0
            steps=0
            global_merge_steps=[]
            # Create list, which steps to merge
            
            # Retrieve all steps with % overlap
            while steps < len(np.unique(sel_end_point_trajectories["Step"])):
                merge_steps = [start]
                steps += 1
                intersect_step = len(np.intersect1d(sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==start)[0]]["Cell"],sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==steps)[0]]["Cell"]))
                total_unqiue = len(np.unique(np.append(sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==start)[0]]["Cell"],sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==steps)[0]]["Cell"])))
                overlap_perc = intersect_step/(total_unqiue+1e-15)
                if overlap_perc > overlap:
                        merge_steps.append(steps)
                while overlap_perc > overlap and steps < len(np.unique(sel_end_point_trajectories["Step"])):
                    steps += 1
                    intersect_step = len(np.intersect1d(sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==start)[0]]["Cell"],sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==steps)[0]]["Cell"]))
                    total_unqiue = len(np.unique(np.append(sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==start)[0]]["Cell"],sel_end_point_trajectories[np.where(sel_end_point_trajectories["Step"]==steps)[0]]["Cell"])))
                    overlap_perc = intersect_step/total_unqiue
                    if overlap_perc > overlap:
                        merge_steps.append(steps)
                    
                start = steps
                global_merge_steps.append(merge_steps)              

            counter=0
            # Merge all adjacent steps
            for h in range(len(global_merge_steps)):
                for k in range(len(global_merge_steps[h])):
                    end_p = cells_along_trajectories["End point"] == end_point
                    average_t = cells_along_trajectories["Trajectory"] == i
                    step_av = cells_along_trajectories["Step"] == k + counter
                    end_p_avg = np.zeros((len(end_p)))

                    for l in range(len(end_p_avg)):
                        end_p_avg[l] = end_p[l] and average_t[l] and step_av[l]
                        cells_along_trajectories["Step"][np.where(end_p_avg)[0]] = h                   
                counter+=len(global_merge_steps[h])

            # Save new bags of cell under separate key
            for m in range(len(global_merge_steps)):
                end_p = cells_along_trajectories["End point"]==end_point
                average_t = cells_along_trajectories["Trajectory"]==i
                step_av = cells_along_trajectories["Step"]==m
                end_p_avg = np.zeros((len(end_p)))

                for l in range(len(end_p_avg)):
                    end_p_avg[l] = end_p[l] and average_t[l] and step_av[l]

                cells = np.unique(cells_along_trajectories["Cell"][np.where(end_p_avg)[0]])
                cell_along_traj = cells_along_trajectories[np.where(end_p_avg)[0]]

                for n in range(len(cells)):
                    max_a = np.max(cell_along_traj[np.where(cell_along_traj["Cell"]==cells[n])]["Allignment Score"])
                    cells_along_trajectories_merged.append([end_point, i, m, cells[n], max_a])
                    
    adata.uns['trajectories']["cells_along_trajectories_each_step_merged"] = np.rec.fromrecords(cells_along_trajectories_merged, 
                                                                                                dtype=[('End point', 'U8'), 
                                                                                                       ('Trajectory', int), 
                                                                                                       ('Step', int), 
                                                                                                       ('Cell', 'U8'), 
                                                                                                       ('Allignment Score', float)])

