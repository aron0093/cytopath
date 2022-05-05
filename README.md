![PyPI](https://img.shields.io/pypi/v/cytopath?color=informational) [![Downloads](https://pepy.tech/badge/cytopath)](https://pepy.tech/project/cytopath) [![Documentation Status](https://readthedocs.org/projects/cytopath/badge/?version=latest)](https://cytopath.readthedocs.io/en/latest/?badge=latest)

## Cytopath
Cytopath is a method for trajectory inference that takes advantage of transcriptional activity information from RNA velocity of single cells by defining a Markov chain model, simulation of an ensemble of possible differentiation trajectories. The preprint can be found [here](https://www.biorxiv.org/content/10.1101/2020.12.21.423801v5).

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree structure) are placed on the inference as each trajectory is modelled independently. Number of trajectories to be inferred can either be defined or estimated in an unsupervised fashion. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process.

![method_figure](https://user-images.githubusercontent.com/25486108/166925895-25fde8d1-c25f-4927-93ad-0331871ef319.png)

## Installation
``` pip install cytopath```

cytopath depends on *scvelo* to process the data and may require additional dependencies of scvelo to be installed manually. Other dependencies are the following:

* python>=3.7
* numpy>=1.20.0
* scipy
* anndata
* scvelo>=0.1.25
* joblib
* fastdtw
* hausdorff
* scitkit-network
* tqdm

## Documentation
Sample notebooks can be found [here](https://github.com/aron0093/cytopath-notebooks). Step by step, installation and analysis guide can be found [here](https://cytopath.readthedocs.io/en/latest/).




