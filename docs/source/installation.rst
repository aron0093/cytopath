Installation
============

Cytopath requires Python 3.7 or later. We recommend using `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ and setting up a new conda environment.

.. code-block:: console

   conda create -y -n cytopath_env python==3.8
   conda activate cytopath_env

PyPI
----

Install cytopath from PyPI using pip

.. code-block:: console

   pip install cytopath
   
Dependencies
------------

Direct dependencies of Cytopath will be installed automatically in the step above however additional dependencies of scvelo will need to be installed manually.

.. code-block:: console

   conda install -c conda-forge python-igraph louvain
   
Cytopath depends upon scikit-network for clustering of simulations. This dependency has been replaced in the `development version <https://github.com/aron0093/cytopath/tree/development>`_. Currently, we reccomend using the following version.


.. code-block:: console

   pip install scikit-network==0.23.1
   
Jupyter notebook
----------------

To run the tutorials install jupyter notebooks.

.. code-block:: console

   conda install notebook



