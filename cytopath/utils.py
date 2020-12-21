from cytopath import anndata

def read_processed(data_loc):

    adata = anndata.read(data_loc)

    try:
        if type(adata.uns['run_info']['end_point_clusters']) == str:
            adata.uns['run_info']['end_point_clusters'] = [adata.uns['run_info']['end_point_clusters']]
    except:
        pass
    
    try:
        for key in adata.uns['run_info']['trajectories_sample_counts'].keys():
            if type(adata.uns['run_info']['trajectories_sample_counts'][key]) != list:
                adata.uns['run_info']['trajectories_sample_counts'][key] = [adata.uns['run_info']['trajectories_sample_counts'][key]]
    except:
        pass

    try:
        adata.uns['trajectories']['cells_along_trajectories_each_step'] = \
        adata.uns['trajectories']['cells_along_trajectories_each_step'].astype([('End point', 'U8'), ('Trajectory', int), ('Step', int), ('Cell', int), 
                                                                                ('Allignment Score', float)])
    except:
        pass

    try:
        adata.uns['trajectories']['cells_along_trajectories_each_step_merged'] = \
        adata.uns['trajectories']['cells_along_trajectories_each_step_merged'].astype([('End point', 'U8'), ('Trajectory', int), ('Step', int), ('Cell', int), 
                                                                                ('Allignment Score', float)])
    except:
        pass

    return adata