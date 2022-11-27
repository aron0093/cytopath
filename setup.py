from setuptools import setup, find_packages

setup(
  name='cytopath',
  version='0.1.9',
  description='Simulation based inference of differentiation trajectories from RNA velocity fields.',
  license='BSD 3-Clause License',
  packages=find_packages(),
  author = 'Revant Gupta',
  author_email = 'revant.gupta.93@gmail.com',
  url = 'https://github.com/aron0093/cytopath',
  download_url = 'https://github.com/aron0093/cytopath/archive/v_019.tar.gz',
  keywords = ['Trajectory Inference', 'single-cell RNA sequencing', 'RNA velocity'],
  install_requires=[
          'numpy>=1.20.0',
          'scipy',
          'anndata',
          'scvelo>=0.1.25',
          'joblib',
          'fastdtw',
          'hausdorff',
          'tqdm',
	  'networkit',
	  'numba'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ])
