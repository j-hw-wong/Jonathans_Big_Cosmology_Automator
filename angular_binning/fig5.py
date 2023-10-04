import os
import numpy as np
import healpy as hp

import sys
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/angular_binning')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood')
#sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt')

from gaussian_cl_likelihood.python import cosmosis_utils, simulation
from angular_binning import loop_likelihood_nbin, posterior, mask, param_grids, covariance, error_vs_nbin

'''
#Step a)
params = {
    'cosmological_parameters--w': {
        'min': -1.0218,
        'max': -0.9782,
        'steps': 8
    },
    'cosmological_parameters--wa': {
        'min': -0.072,
        'max': 0.072,
        'steps': 8
    }
}
'''
'''
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
output_dir = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

param_grids.reduced_grid_chains(params, upper_diag, lower_diag, n_chains, output_dir)
'''

'''
cosmosis_utils.combine_chain_output(
    input_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/',
    chain_subdir_mask='chain{i}/',
    filemask='_{n}.tgz',
    )

cosmosis_utils.extract_data(
    input_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/',
    filemask='_*.tgz',
    params_filename='cosmological_parameters/values.txt',
    nodelete=False,
    nbin_3x2=3)
'''


#Step b)

n_zbin=3
noise_save_dir = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/noise/'
if not os.path.exists(noise_save_dir):
    os.makedirs(noise_save_dir)

simulation.noise_cls(
    gals_per_sq_arcmin=30 / n_zbin,
    sigma_e=0.3,
    lmin=2,
    lmax=2000,
    save_dir=noise_save_dir)


#Step c)

mask.generate_mask(
    wmap_mask_path='/raid/scratch/wongj/mywork/3x2pt/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits',
    nside=512,
    target_fsky=0.3,
    mask_save_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/mask.fits.gz')

mask.plot_mask(
    mask_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/mask.fits.gz',
    save_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/mask.pdf')



#Step d)
mask.get_3x2pt_mixmats(
    #mask_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/wmap_64.fits',
    mask_path='/raid/scratch/wongj/mywork/3x2pt/Euclid_Masks/Euclid_DR1_1024.fits',
    nside=2048,
    lmin=2,
    lmax_mix=2000,
    lmax_out=2000,
    save_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/mixmats.npz')


#Step e)

theory_cl_dir = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/_42/'

covariance.get_3x2pt_cov(
    n_zbin=3,
    pos_pos_filemask=theory_cl_dir + 'galaxy_cl/bin_{hi_zbin}_{lo_zbin}.txt',
    pos_she_filemask=theory_cl_dir + 'galaxy_shear_cl/bin_{pos_zbin}_{she_zbin}.txt',
    she_she_filemask=theory_cl_dir + 'shear_cl/bin_{hi_zbin}_{lo_zbin}.txt',
    lmax_in=2000,
    lmin_in=2,
    lmax_out=2000,
    lmin_out=2,
    pos_nl_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/noise/pos_nl.txt',
    she_nl_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/noise/she_nl.txt',
    noise_lmin=2,
    mask_path='/raid/scratch/wongj/mywork/3x2pt/Euclid_Masks/Euclid_DR1_1024.fits',
    nside=2048,
    save_filemask='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz')


#Step f)

covariance.bin_combine_cov(
    n_zbin=3,
    lmin_in=2,
    lmin_out=2,
    lmax=2000,
    n_bp_min=1,
    n_bp_max=10,
    input_filemask='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
    save_filemask='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/cov_{n_bp}bp.npz')


#Step g)

fid_base_dir = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/_42/'

loop_likelihood_nbin.like_bp_gauss_mix_loop_nbin(
    grid_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/',
    n_bps=np.arange(1,11),
    n_zbin=3,
    lmax_like=2000,
    lmin_like=2,
    lmax_in=2000,
    lmin_in=2,
    fid_pos_pos_dir=fid_base_dir + 'galaxy_cl/',
    fid_she_she_dir=fid_base_dir + 'shear_cl/',
    fid_pos_she_dir=fid_base_dir + 'galaxy_shear_cl/',
    pos_nl_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/noise/pos_nl.txt',
    she_nl_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/noise/she_nl.txt',
    mixmats_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/mixmats.npz',
    bp_cov_filemask='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/cov_{n_bp}bp.npz',
    binmixmat_save_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/',
    varied_params= ['w', 'wa'],
    like_save_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/',
    create_obs_bp_save_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/',
    obs_bandpowers_dir=None)

'''
posterior.cl_post(
    log_like_filemask=cl_like_filemask,
    contour_levels_sig=[1, 2, 3],
    bp=5,
    colour='C0',
    linestyle='-',
    plot_save_path=None)    
'''

'''
#Step h)

loop_likelihood_nbin.obs_from_fid(
input_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/_42/', 
output_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/_42/obs-cf.npz', 
n_zbin=3, 
lmax=2000, 
lmin=2)

n_zbin = 5

# Low nbin run
#grid_dir_lo = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/'
#n_theta_bins_lo = np.arange(1, 4)

fid_base_dir = '/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/_42/'
n_zbin = 3

loop_likelihood_nbin.like_cf_gauss_loop_nbin(
grid_dir='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/chains/', 
n_theta_bins== np.arange(1,11),
n_zbin = n_zbin,
lmin = 2, 
lmax = 2000, 
theta_min = np.radians(0.1), 
theta_max = np.radians(10), 
fid_pos_pos_dir=fid_base_dir + 'galaxy_cl/', 
fid_she_she_dir=fid_base_dir + 'shear_cl/', 
fid_pos_she_dir=fid_base_dir + 'galaxy_shear_cl/', 
obs_path='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/obs-cf.npz', 
survey_area_sqdeg=15000, 
gals_per_sqarcmin_per_zbin=30 / n_zbin, 
sigma_e=0.3, 
varied_params=['w', 'wa'], 
like_save_dir='/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/')

#Step i)

error_vs_nbin.area_vs_nbin(
cl_like_filemask='/raid/scratch/wongj/mywork/3x2pt/TEST_3_BINS_NEW/inference_default/like_lmaxlike2000_{n_bp}bp.txt', 
cf_like_filemask='/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5_/like_thetamin0.1_{n_bin}bins.txt', 
contour_levels_sig=np.arange(1, 3, 100), 
n_bps=np.arange(1, 11), 
n_theta_bins=np.arange(1, 11), 
save_path=None)


'''


