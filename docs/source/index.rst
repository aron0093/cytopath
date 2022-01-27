Cytopath: Simulation based inference of differentiation trajectories from RNA velocity fields
===================================

Cytopath is a method for trajectory inference that takes advantage of transcriptional activity information from RNA velocity of single cells by defining a Markov chain model, simulation of an ensemble of possible differentiation trajectories. The preprint can be found `here <https://www.biorxiv.org/content/10.1101/2020.12.21.423801v4.full>`_.

Cytopath can infer trajectories with or without root/terminal state supervision. No topological constraints (e.g. a tree structure) are placed on the inference as each trajectory is modelled independently. Number of trajectories to be inferred can either be defined or estimated in an unsupervised fashion. Subsequent statistical analysis reveals the topological and molecular characteristics of the differentiation process.

Check out the :doc:`notebooks` section for demonstration of cytopath on publicly available datasets.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
   notebooks
