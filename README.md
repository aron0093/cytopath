![PyPI](https://img.shields.io/pypi/v/cytopath?color=informational) [![Downloads](https://pepy.tech/badge/cytopath)](https://pepy.tech/project/cytopath) [![Documentation Status](https://readthedocs.org/projects/cytopath/badge/?version=latest)](https://cytopath.readthedocs.io/en/latest/?badge=latest)

## Cytopath
Cytopath is a method for trajectory inference with single-cell RNA sequencing data. Transcriptional activity information from RNA velocity of single cells is used to define a Markov chain model; simulation of this model yields an ensemble of possible differentiation trajectories that are used to estimate the lineage path.

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree structure) are placed on the inference as each trajectory is modelled independently. The number of trajectories to be inferred can either be defined or estimated in an unsupervised fashion. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process. 

Cytopath can model complex behaviours like cycling and convergence as well as cooccurring combinations of multiple processes.

![method_figure](https://user-images.githubusercontent.com/25486108/166925895-25fde8d1-c25f-4927-93ad-0331871ef319.png)

## Installation
``` pip install cytopath```

cytopath depends on *scvelo* to process the data and may require additional dependencies of scvelo to be installed manually. Other dependencies are the following:
```
python==3.9 numpy==1.23.5 scipy pandas>=1.4.1 anndata scvelo>=0.1.25 joblib fastdtw hausdorff scikit-learn>=1.30.0 numba tqdm
```
## Documentation
Sample notebooks can be found [here](https://github.com/aron0093/cytopath-notebooks). A step-by-step installation and analysis guide can be found [here](https://cytopath.readthedocs.io/en/latest/).

## Citation

Gupta R., Cerletti D., Gut G., Oxenius A., Claassen M. Simulation-based inference of differentiation trajectories from RNA velocity fields. Cell Reports Methods,
Volume 2, Issue 12, 2022, 100359, ISSN 2667-2375. [doi:https://doi.org/10.1016/j.crmeth.2022.100359](https://doi.org/10.1016/j.crmeth.2022.100359).





