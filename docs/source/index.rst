Simulation-based inference of differentiation trajectories from RNA velocity fields
=============================================================================================

.. note::

   This project is under active development.

Cytopath is a method for trajectory inference with single-cell RNA sequencing data. Transcriptional activity information from RNA velocity of single cells is used to define a Markov chain model; simulation of this model yields an ensemble of possible differentiation trajectories that are used to estimate the lineage path.

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree structure) are placed on the inference as each trajectory is modelled independently. The number of trajectories to be inferred can either be defined or estimated in an unsupervised fashion. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process.

Cytopath can model complex behaviours like cycling and convergence as well as co-occurring combinations of multiple processes. Read more `here <https://doi.org/10.1016/j.crmeth.2022.100359>`_.

Check out the :doc:`notebooks` section for a demonstration of Cytopath on publicly available datasets.

.. image:: https://user-images.githubusercontent.com/25486108/166925895-25fde8d1-c25f-4927-93ad-0331871ef319.png

.. toctree::
   :maxdepth: 1
   :hidden:
   
   installation
   usage
   notebooks
