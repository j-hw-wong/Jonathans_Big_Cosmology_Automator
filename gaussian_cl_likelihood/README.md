## Code used to produce all figures in [Upham, Brown & Whittaker 2021, “Sufficiency of a Gaussian power spectrum likelihood for accurate cosmology from upcoming weak lensing surveys”, arXiv:2012.06267](https://arxiv.org/abs/2012.06267)
[NOTE: Updated and managed by J. WONG].

### Documentation

Documentation is available at [Read the Docs](https://gaussian-cl-likelihood.readthedocs.io/), including:

* Steps to produce all figures from the paper.

* Full reference of all files and functions.

### Prerequisites

* Python 3
* numpy
* scipy
* matplotlib
* [NPEET](https://github.com/gregversteeg/NPEET) (only needed for `mutual_info` module)
* [healpy](https://healpy.readthedocs.io/en/latest/install.html) (only needed for `simulation` module)
* [pseudo_cl_likelihood](https://github.com/robinupham/pseudo_cl_likelihood) (only needed for `pcl_like` module)
* [NaMaster](https://namaster.readthedocs.io/en/latest/installation.html) (only needed for flat-sky functionality within `simulation` module)
