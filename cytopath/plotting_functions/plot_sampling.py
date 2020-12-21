import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
      
def frame_gen(i):
    plt.scatter(adata.obsm['X_'+basis][np.unique(samples[:, i-1]),0], adata.obsm['X_'+basis][np.unique(samples[:, i-1]),1], c='blue', s=100)

def pl_cytopath_sampling(adata, basis="umap", folder=""):

    if folder[-1] != '/':
        folder+='/'

    fig = plt.figure()
    plt.scatter(adata.obsm['X_'+basis][:,0], adata.obsm['X_'+basis][:,1], c='grey', s=100)
    samples = adata.uns['samples']['cell_sequences']
    animator = ani.FuncAnimation(fig, frame_gen, frames= samples.shape[1], interval = 200)
    animator.save(folder+'sampling.gif')

    