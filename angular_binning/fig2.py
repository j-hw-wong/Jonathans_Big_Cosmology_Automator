import os
import numpy as np
import healpy as hp

import sys
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/angular_binning')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood')
#sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt')

from gaussian_cl_likelihood.python import cosmosis_utils, simulation
from angular_binning import loop_likelihood_nbin, posterior

'''
#Step a)
params = {
    'cosmological_parameters--w': {
        'min': -1.0119,
        'max': -0.9881,
        'steps': 49
    },
    'cosmological_parameters--wa': {
        'min': -0.039,
        'max': 0.039,
        'steps': 45
    }
}

n_chains = 11 # for a 12-CPU machine, to leave one free
output_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

cosmosis_utils.generate_chain_input(params, n_chains, output_dir)
'''

'''
input_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/chains/'
chain_subdir_mask = 'chain{i}/'
filemask = '_{n}.tgz'
clean = False

cosmosis_utils.combine_chain_output(input_dir, chain_subdir_mask, filemask, clean)
'''


'''
input_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/chains/'
filemask = '_*.tgz'
params_filename = 'cosmological_parameters/values.txt'
nodelete = True # or True if you want to retain the .tgz files in case you want to extract anything else later
nbin_3x2 = 5

cosmosis_utils.extract_data(input_dir, filemask=filemask, params_filename=params_filename, nodelete=nodelete, nbin_3x2=nbin_3x2)
'''



'''
#Step b)
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
n_chains = 24
output_dir = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/b/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


param_grids.reduced_grid_chains(params, upper_diag, lower_diag, n_chains, output_dir=output_dir)
'''

'''
input_dir = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/b/chains/' 
chain_subdir_mask = 'chain{i}/'
filemask = '_{n}.tgz'
clean = True

cosmosis_utils.combine_chain_output(input_dir, chain_subdir_mask, filemask, clean)
'''

'''
input_dir = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/b/chains/'
filemask = '_*.tgz'
params_filename = 'cosmological_parameters/values.txt'
nodelete = False # or True if you want to retain the .tgz files in case you want to extract anything else later
nbin_3x2 = 5

cosmosis_utils.extract_data(input_dir, filemask=filemask, params_filename=params_filename, nodelete=nodelete, nbin_3x2=nbin_3x2)
'''

'''
#Step c)
n_zbin = 5
gals_per_sq_arcmin = 30 / n_zbin
sigma_e = 0.3
lmin = 2
lmax = 2000
#lmin = 2
#lmax = 75
save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/noise/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

simulation.noise_cls(gals_per_sq_arcmin, sigma_e, lmin, lmax, save_dir)
'''

'''
#Step d)
input_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/chains/_1000/'
output_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/mock.npz'
n_zbin = 5
lmax = 2000
lmin =  2

loop_likelihood_nbin.obs_from_fid(input_dir, output_path, n_zbin, lmax, lmin)
'''

'''
#Step e)
grid_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/chains/'
#n_bps = np.arange(2, 5)
#n_bps = np.array([5,10,15])
n_bps = [1, 5, 10, 20]
#n_bps = [20]
n_zbin = 5
lmax = 2000
#lmin_like = 10
lmin_like = 2
lmin_in = 2
fid_base_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/a/chains/_1000/'
fid_pos_pos_dir = fid_base_dir + 'galaxy_cl/'
fid_she_she_dir = fid_base_dir + 'shear_cl/'
fid_pos_she_dir = fid_base_dir + 'galaxy_shear_cl/'
pos_nl_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/noise/pos_nl.txt'
she_nl_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/noise/she_nl.txt'
noise_ell_path = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/noise/noise_ell.txt'
pbl_save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/'
obs_bp_save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/'
inv_cov_save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/'
varied_params = ['w', 'wa']
like_save_dir = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/log_like_ps/'

if not os.path.exists(like_save_dir):
    os.makedirs(like_save_dir)

loop_likelihood_nbin.like_bp_gauss_loop_nbin(grid_dir, n_bps, n_zbin, lmax, lmin_like, lmin_in, fid_pos_pos_dir, fid_she_she_dir, fid_pos_she_dir, pos_nl_path, she_nl_path, noise_ell_path, pbl_save_dir, obs_bp_save_dir, inv_cov_save_dir, varied_params, like_save_dir)
'''


log_like_filemask = '/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/log_like_ps/like_lmax2000_{n_bp}bp.txt'
contour_levels_sig = [1, 2, 3]
n_bps = [[20, 1], [20, 5], [20, 10]]
#n_bps = [[5,5]]
colours = ['C0', 'C1']
linestyles = [['-'], [(0, (2, 2))]]
#plot_save_path = 'path-to-save-figure.pdf' # or None to just display it
plot_save_path = None

posterior.cl_posts(log_like_filemask, contour_levels_sig, n_bps, colours, linestyles, plot_save_path=plot_save_path)
#posterior.cl_post(log_like_filemask, contour_levels_sig, n_bps, colours, linestyles, plot_save_path=plot_save_path)
posterior.cl_post(
    log_like_filemask='/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/log_like_ps/like_lmax2000_5bp.txt',
    contour_levels_sig=[1, 2, 3],
    bp=5,
    colour='C0',
    linestyle='-',
    plot_save_path='/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig2/s_plot.pdf')
    #plot_save_path=None)



'''
# It is necessary to split this into two different grids
# because the performance deteriorates so much for low nbin

# Common parameters
n_zbin = 5
lmin = 2
lmax = 2000
theta_min = np.radians(0.1) # 0.1 deg
theta_max = np.radians(10)  # 10 deg
fid_base_dir = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/a/chains/_1000/'
fid_pos_pos_dir = fid_base_dir + 'galaxy_cl/'
fid_she_she_dir = fid_base_dir + 'shear_cl/'
fid_pos_she_dir = fid_base_dir + 'galaxy_shear_cl/'
obs_path = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/mock.npz'
survey_area_sqdeg = 15000
gals_per_sqarcmin_per_zbin = 30 / n_zbin
sigma_e = 0.3
varied_params = ['w', 'wa']
like_save_dir = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/log_like/'

# Low nbin run
grid_dir_lo = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/b/chains/'
n_theta_bins_lo = np.arange(1, 6)
loop_likelihood_nbin.like_cf_gauss_loop_nbin(grid_dir_lo, n_theta_bins_lo, n_zbin, lmin, lmax, theta_min, theta_max, fid_pos_pos_dir, fid_she_she_dir, fid_pos_she_dir, obs_path, survey_area_sqdeg, gals_per_sqarcmin_per_zbin, sigma_e, varied_params, like_save_dir)

# High nbin run
grid_dir_hi = '/scratch/nas_mberc2/wongj/mywork/ang_bin_fig2/a/chains/'
n_theta_bins_hi = np.arange(6, 31)
loop_likelihood_nbin.like_cf_gauss_loop_nbin(grid_dir_hi, n_theta_bins_hi, n_zbin, lmin, lmax, theta_min, theta_max, fid_pos_pos_dir, fid_she_she_dir, fid_pos_she_dir, obs_path, survey_area_sqdeg, gals_per_sqarcmin_per_zbin, sigma_e, varied_params, like_save_dir)
'''

