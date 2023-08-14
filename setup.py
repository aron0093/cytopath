from setuptools import setup, find_packages

setup(
  name='cytopath',
  version='0.2.1',
  description='Simulation based inference of differentiation trajectories from RNA velocity fields.',
  license='BSD 3-Clause License',
  packages=find_packages(),
  author = 'Revant Gupta',
  author_email = 'revant.gupta.93@gmail.com',
  url = 'https://github.com/aron0093/cytopath',
  download_url = 'https://github.com/aron0093/cytopath/archive/v_021.tar.gz',
  keywords = ['Trajectory Inference', 'single-cell RNA sequencing', 'RNA velocity'],
  install_requires=[
          'scvelo>=0.1.25',
          'numpy==1.23.5',
          'joblib',
          'fastdtw',
          'hausdorff',
          'scikit-learn>=1.3.0',
          'tqdm',
	        'numba'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.9',
  ])
