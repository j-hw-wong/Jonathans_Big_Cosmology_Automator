import os
import numpy as np
import healpy as hp

import sys
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/angular_binning')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood')
#sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt')

from gaussian_cl_likelihood.python import cosmosis_utils, simulation
from angular_binning import loop_likelihood_nbin, posterior, mask, param_grids


#Step a)
params = {
    'cosmological_parameters--w': {
        'min': -1.0218,
        'max': -0.9782,
        'steps': 49
    },
    'cosmological_parameters--wa': {
        'min': -0.072,
        'max': 0.072,
        'steps': 49
    }
}
upper_diag = [(-0.99, 0.072), (-0.9782, 0.02)]
lower_diag = [(-1.0218, -0.02), (-1.01, -0.072)]
n_chains = 11
output_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/chains/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


param_grids.reduced_grid_chains(params, upper_diag, lower_diag, n_chains, output_dir=output_dir)


'''
input_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/chains/' 
chain_subdir_mask = 'chain{i}/'
filemask = '_{n}.tgz'
clean = True

cosmosis_utils.combine_chain_output(input_dir, chain_subdir_mask, filemask, clean)
'''

'''
input_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/chains/'
filemask = '_*.tgz'
params_filename = 'cosmological_parameters/values.txt'
nodelete = False # or True if you want to retain the .tgz files in case you want to extract anything else later
nbin_3x2 = 5

cosmosis_utils.extract_data(input_dir, filemask=filemask, params_filename=params_filename, nodelete=nodelete, nbin_3x2=nbin_3x2)
'''

'''
n_zbin = 5
gals_per_sq_arcmin = 30 / n_zbin
sigma_e = 0.3
#lmin = 2
#lmax = 2000
lmin = 2
lmax = 2000
save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/noise/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

simulation.noise_cls(gals_per_sq_arcmin, sigma_e, lmin, lmax, save_dir)
'''


'''
wmap_mask_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits' # available from https://lambda.gsfc.nasa.gov/product/map/dr5/masks_get.cfm
nside = 512
target_fsky = 0.3
mask_save_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/mask.fits.gz'

mask.generate_mask(wmap_mask_path, nside, target_fsky, mask_save_path)

mask_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/mask.fits.gz'
save_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/mask.pdf'

mask.plot_mask(mask_path, save_path)
'''