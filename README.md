## cytopath
Cytopath is a method for trajectory inference that takes advantage of transcriptional activity information from RNA velocity of single cells by defining a Markov chain model, simulation of an ensemble of possible differentiation trajectories. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process. The preprint can be found [here](https://www.biorxiv.org/content/10.1101/2020.12.21.423801v1).

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree stucuture) are placed on the inference. Number of trajectories to estimated can either be defined or inferred in an unsupervised fashion.

## Installation
``` pip install cytopath```

cytopath depends on scvelo to process the data and may require additional dependencies of scvelo to be installed manually.

## Sample notebooks
Sample notebooks can be found [here](https://github.com/aron0093/cytopath-notebooks).




