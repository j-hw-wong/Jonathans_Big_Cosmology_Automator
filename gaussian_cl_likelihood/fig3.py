from python import cosmosis_utils
from python import posteriors
from python import run_likelihoods
from python import simulation

import os
import configparser
import numpy as np

pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']

config = configparser.ConfigParser()
config.read(pipeline_variables_path)

#Step a) 1
params = {
    'cosmological_parameters--w': {
        'min': -1.005,
        'max': -0.995,
        'steps': 5
    },
    'cosmological_parameters--wa': {
        'min': -0.01,
        'max': 0.01,
        'steps': 5
    },
    'cosmological_parameters--omega_m': {
        'min': 0.3133,
        'max': 0.3143,
        'steps': 5
    }
}

'''
n_chains = 11 # for a 12-CPU machine, to leave one free
#output_dir can be changed to $SAVE_DIR from run_3x2pt_tomo.sh (?)
#output_dir = '/raid/scratch/wongj/mywork/3x2pt/JWONG_LEGEND'
output_dir = os.environ['SAVE_DIR'] + os.environ['INFERENCE_OUTPUT_DIR']
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
cosmosis_utils.generate_chain_input(params, n_chains, output_dir)
'''


'''
#Step a) 4
input_dir = os.environ['SAVE_DIR'] + os.environ['INFERENCE_OUTPUT_DIR'] + 'chains/'
chain_subdir_mask = 'chain{i}/'
filemask = '_{n}.tgz'
clean = False

cosmosis_utils.combine_chain_output(input_dir, chain_subdir_mask, filemask, clean)
'''

'''
#Step a) 5
input_dir = os.environ['SAVE_DIR'] + os.environ['INFERENCE_OUTPUT_DIR'] + 'chains/'
filemask = '_*.tgz'
params_filename = 'cosmological_parameters/values.txt'
nodelete = False # or True if you want to retain the .tgz files in case you want to extract anything else later
nbin_3x2 = 5

cosmosis_utils.extract_data(input_dir, filemask=filemask, params_filename=params_filename, nodelete=nodelete, nbin_3x2=nbin_3x2)
'''

'''
#Step b)
#n_zbin = 5
n_zbin = int(config['global']['NBINS'])
gals_per_sq_arcmin = 30 / n_zbin
sigma_e = 0.3
lmin = float(config['global']['ELL_MIN'])
lmax = float(config['global']['ELL_MAX'])
#lmin = 2
#lmax = 2000
save_dir = os.environ['SAVE_DIR'] + os.environ['INFERENCE_OUTPUT_DIR'] + 'noise_cls/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

simulation.noise_cls(gals_per_sq_arcmin, sigma_e, lmin, lmax, save_dir)
'''

#NOTE TO FUTURE JWONG -
#GOT UP TO HERE. PROBABLY WILL HAVE TO POINT THE OBS_CLS TO WHERE MY RECOV_CAT_CLS OUTPUT IS
#WILL HAVE TO DOUBLE CHECK FORMATTING ON OBS_CLS.TXT TO MAKE SURE MY PIPELINE IS IN SAME FORMAT
#ALSO CHECK MASK AND HOW THIS AFFECTS GAUSSIAN_CL_LIKELIHOOD
#PAIN.LOL.
'''
n_zbin = 5
fid_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/chains/_0'
pos_pos_in_dir = fid_dir + '/galaxy_cl'
she_she_in_dir = fid_dir + '/shear_cl'
pos_she_in_dir = fid_dir + '/galaxy_shear_cl'
lmax = 2000
lmin_in = 2
pos_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/noise_cls/pos_nl.txt'
she_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/noise_cls/she_nl.txt'
nside = 1024
lmin_out = 2
mask_path = None # full sky
obs_cls_save_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/obs_cls.txt'

simulation.single_obs_cls(n_zbin, pos_pos_in_dir, she_she_in_dir, pos_she_in_dir, lmax, lmin_in, pos_nl_path, she_nl_path, nside, lmin_out, obs_cls_save_dir, mask_path)

obs_cl_arr = np.loadtxt(obs_cls_save_dir)
simulation.save_cls_nob(obs_cl_arr,n_zbin=5,save_dir='/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/obs_cls/',lmin=lmin_out)
'''

'''
grid_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/chains/'
varied_params = ['w', 'wa', 'omega_m']
save_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/output_likelihood.txt'
n_zbin = 5
obs_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/obs_cls'
obs_pos_pos_dir = obs_dir + '/galaxy_cl'
obs_she_she_dir = obs_dir + '/shear_cl'
obs_pos_she_dir = obs_dir + '/galaxy_shear_cl'
pos_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/noise_cls/pos_nl.txt'
she_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/noise_cls/she_nl.txt'
noise_ell_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/noise_cls/noise_ell.txt'
fid_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/chains/_0'
fid_pos_pos_dir = fid_dir + '/galaxy_cl'
fid_she_she_dir = fid_dir + '/shear_cl'
fid_pos_she_dir = fid_dir + '/galaxy_shear_cl'
lmax = 2000

run_likelihoods.run_likes_cl_wishart_gauss(grid_dir, varied_params, save_path, n_zbin, obs_pos_pos_dir, obs_she_she_dir, obs_pos_she_dir, pos_nl_path, she_nl_path, noise_ell_path, fid_pos_pos_dir, fid_she_she_dir, fid_pos_she_dir, lmax)
'''

'''
#Step e)
log_like_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/output_likelihood.txt'
save_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/likelihood_grid.npz'

posteriors.grid_3d(log_like_path, save_path)
'''

'''
#Step f)
grid_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig3/likelihood_grid.npz'
contour_levels_sig = [1, 2, 3]
# All the optional parameters can be tweaked as required

posteriors.plot_3d(grid_path, contour_levels_sig, like_label_1='Wishart',like_label_2='Gaussian')
'''