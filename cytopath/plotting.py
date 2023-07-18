from .plotting_functions.plot_alignment import pl_cytopath_alignment
from .plotting_functions.plot_sampling import pl_cytopath_sampling
from .plotting_functions.plot_radial_heatmap import radial_heatmap

def plot_trajectories(data, basis='', smoothing=False, figsize=(15,4), size=3, 
                      show=True, save=False, save_type='png', directory=''):
    pl_cytopath_alignment(data, basis=basis, smoothing=smoothing, figsize=figsize, size=size, 
                          show=show, save=save, save_type=save_type, folder=directory,)

def plot_sampling(data, directory='', basis=''):
    pl_cytopath_sampling(data, folder=directory, basis=basis)
