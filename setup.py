from setuptools import setup, find_packages

setup(
  name='cytopath',
  version='0.1.3',
  description='Trajectory inference for single cell RNA seq. data using transition probabiloities derived from RNA velocity of single cells.',
  license='BSD 3-Clause License',
  packages=find_packages(),
  author = 'Revant Gupta',                   # Type in your name
  author_email = 'revant.gupta.93@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/aron0093/cytopath',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/aron0093/cytopath/archive/v_01_03.tar.gz',    # I explain this later on
  keywords = ['Trajectory Inference', 'single-cell RNA sequencing', 'RNA velocity'],   # Keywords that define your package best
  install_requires=[
          'numpy',
          'scipy',
          'anndata',
          'scvelo>=0.1.25',
          'joblib',
          'fastdtw',
          'hausdorff',
          'tqdm'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: BSD License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
  ])