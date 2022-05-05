Cytopath: Simulation based inference of differentiation trajectories from RNA velocity fields
=============================================================================================

.. note::

   This project is under active development.

Cytopath is a method for trajectory inference that takes advantage of transcriptional activity information from RNA velocity of single cells by defining a Markov chain model, simulating possible differentiation paths and inferring an ensemble trajectory. The preprint can be found `here <https://www.biorxiv.org/content/10.1101/2020.12.21.423801v5>`_.

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree structure) are placed on the inference as each trajectory is modelled independently. Number of trajectories to be inferred can either be defined or estimated in an unsupervised fashion. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process.

Check out the :doc:`notebooks` section for demonstration of cytopath on publicly available datasets.

.. image:: https://user-images.githubusercontent.com/25486108/166925895-25fde8d1-c25f-4927-93ad-0331871ef319.png

.. toctree::
   :maxdepth: 1
   :hidden:
   
   installation
   usage
   notebooks
