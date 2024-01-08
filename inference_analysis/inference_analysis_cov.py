import configparser
import numpy as np
import os
import sys

angular_binning_path = os.environ['ANGULAR_BINNING_PATH']
gaussian_cl_likelihood_path = os.environ['GAUSSIAN_CL_LIKELIHOOD_PATH']
pipeline_dir = os.environ['PIPELINE_DIR']

sys.path.insert(1, angular_binning_path)
sys.path.insert(1, gaussian_cl_likelihood_path)
sys.path.insert(1, pipeline_dir)

from gaussian_cl_likelihood.python import cosmosis_utils, simulation, posteriors
from angular_binning import loop_likelihood_nbin, posterior, mask, param_grids, covariance, error_vs_nbin, \
    like_bp_gauss_mix, covariance_1x2pt

pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']

config = configparser.ConfigParser()
config.read(pipeline_variables_path)

save_dir = str(config['inference_analysis_params']['MEASUREMENT_SAVE_DIR'])
nbins = int(config['inference_analysis_params']['N_ZBIN'])
nside = int(config['inference_analysis_params']['NSIDE'])
input_lmin = int(float(config['inference_analysis_params']['INPUT_ELL_MIN']))
input_lmax = int(float(config['inference_analysis_params']['INPUT_ELL_MAX']))
output_lmin = int(float(config['inference_analysis_params']['OUTPUT_ELL_MIN']))
output_lmax = int(float(config['inference_analysis_params']['OUTPUT_ELL_MAX']))
n_bandpowers = int(config['inference_analysis_params']['N_BANDPOWERS'])

zmin = float(config['inference_analysis_params']['ZMIN'])
zmax = float(config['inference_analysis_params']['ZMAX'])
obs_type = str(config['inference_analysis_params']['OBS_TYPE'])
ngal = str(config['inference_analysis_params']['NGAL'])

map_lmin = 0
map_lmax = (3*nside)-1

inference_analysis_output_dir = save_dir+'inference_analysis/'

chains_input_dir = inference_analysis_output_dir + 'chains/'

mask_dir = str(config['inference_analysis_params']['PATH_TO_MASK'])
mix_mats_save_path = save_dir + 'inference_analysis/mixmats.npz'

cov_fromsim_path = save_dir + 'cov_fromsim/'

theory_cl_dir = save_dir + 'theory_cls/'
noise_save_dir = save_dir + 'measured_noise_cls/'

binmixmat_save_dir = inference_analysis_output_dir + 'bin_mix_mats/'
if not os.path.exists(binmixmat_save_dir):
    os.makedirs(binmixmat_save_dir)

cl_like_filemask = inference_analysis_output_dir + 'like_lmaxlike%s_{n_bp}bp.txt' % (output_lmax)


# Unpack CosmoSIS data from grid
cosmosis_utils.combine_chain_output(input_dir=chains_input_dir,
                                    chain_subdir_mask='chain{i}/',
                                    filemask='_{n}.tgz'
                                    )

cosmosis_utils.extract_data(input_dir=chains_input_dir,
                            filemask='_*.tgz',
                            params_filename='cosmological_parameters/values.txt',
                            nodelete=False,
                            nbin_3x2=nbins)


# Calculate mixing matrices from mask
mask.get_3x2pt_mixmats(mask_path=mask_dir,
                       nside=nside,
                       lmin=input_lmin,
                       lmax_mix=map_lmax,
                       lmax_out=output_lmax,
                       input_lmax=input_lmax,
                       save_path=mix_mats_save_path)


# Bin PCl covariance matrix to form Bandpower covariance matrix

if obs_type == '1X2PT':

    field = str(config['inference_analysis_params']['FIELD'])


    # Bin the covariance measured from simulations
    covariance.bin_combine_cov(n_zbin=nbins,
                               lmin_in=output_lmin,
                               lmin_out=output_lmin,
                               lmax=output_lmax,
                               n_bp_min=n_bandpowers,
                               n_bp_max=n_bandpowers,
                               input_filemask=cov_fromsim_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                               save_filemask=cov_fromsim_path + 'cov_{n_bp}bp.npz',
                               bandpower_spacing='log',
                               obs_type='1x2pt')
    

    # Calculate covariance of signal part of measurement
    save_signal_cov_path = inference_analysis_output_dir + 'signal_cov/'
    if not os.path.exists(save_signal_cov_path):
        os.makedirs(save_signal_cov_path)

    covariance_1x2pt.get_1x2pt_cov_signal_pcl(
        n_zbin=nbins,
        pos_pos_filemask=theory_cl_dir + 'galaxy_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        pos_she_filemask=theory_cl_dir + 'galaxy_shear_cl/bin_{pos_zbin}_{she_zbin}.txt',
        she_she_filemask=theory_cl_dir + 'shear_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        field=field,
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        lmax_out=output_lmax,
        lmin_out=output_lmin,
        mask_path=mask_dir,
        nside=nside,
        save_filemask=save_signal_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz')

    # Calculate covariance of noise part of measurement

    save_noise_cov_path = inference_analysis_output_dir + 'noise_cov/'
    if not os.path.exists(save_noise_cov_path):
        os.makedirs(save_noise_cov_path)

    covariance_1x2pt.get_1x2pt_cov_noise(
        n_zbin=nbins,
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        lmax_out=output_lmax,
        lmin_out=output_lmin,
        pos_nl_path=noise_save_dir + 'galaxy_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        she_nl_path=noise_save_dir + 'shear_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        pos_she_nl_path=noise_save_dir + 'galaxy_shear_cl/bin_{pos_zbin}_{she_zbin}.txt',
        field=field,
        noise_lmin=output_lmin,
        nside=nside,
        save_filemask=save_noise_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz')

    # Step f)
    # Bin the signal covariance into bandpowers
    covariance_1x2pt.bin_combine_cov(n_zbin=nbins,
                                     lmin_in=output_lmin,
                                     lmin_out=output_lmin,
                                     lmax=output_lmax,
                                     n_bp_min=n_bandpowers,
                                     n_bp_max=n_bandpowers,
                                     input_filemask=save_signal_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                                     save_filemask=save_signal_cov_path + 'cov_{n_bp}bp.npz')

    # Bin the noise covariance into bandpowers
    covariance_1x2pt.bin_combine_cov(n_zbin=nbins,
                               lmin_in=output_lmin,
                               lmin_out=output_lmin,
                               lmax=output_lmax,
                               n_bp_min=n_bandpowers,
                               n_bp_max=n_bandpowers,
                               input_filemask=save_noise_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                               save_filemask=save_noise_cov_path + 'cov_{n_bp}bp.npz')

    save_combined_cov_path = inference_analysis_output_dir + 'combined_signal_noise_cov/'
    if not os.path.exists(save_combined_cov_path):
        os.makedirs(save_combined_cov_path)

    # Combine the individual signal and noise covariance blocks
    covariance_1x2pt.combine_signal_noise_cov_spec(
        n_zbin=nbins,
        field=field,
        signal_cov_path=save_signal_cov_path,
        noise_cov_path=save_noise_cov_path,
        combined_cov_save_path=save_combined_cov_path)

    # Combine the overall covariance matrix
    covariance.combine_signal_noise_cov(
        n_zbin=nbins,
        lmin_in=output_lmin,
        lmin_out=output_lmin,
        lmax=output_lmax,
        n_bp_min=n_bandpowers,
        n_bp_max=n_bandpowers,
        signal_cov_path=save_signal_cov_path + 'cov_{n_bp}bp.npz',
        noise_cov_path=save_noise_cov_path + 'cov_{n_bp}bp.npz',
        combined_cov_save_path=save_combined_cov_path + 'cov_{n_bp}bp.npz')

    save_combined_cov_path = inference_analysis_output_dir + 'combined_signal_noise_cov/'
    if not os.path.exists(save_combined_cov_path):
        os.makedirs(save_combined_cov_path)

    loop_likelihood_nbin.like_bp_gauss_mix_loop_nbin_1x2pt(
        grid_dir=chains_input_dir,
        n_bps=np.array([n_bandpowers]),
        n_zbin=nbins,
        lmax_like=output_lmax,
        lmin_like=output_lmin,
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        field=field,
        #fid_she_she_dir=theory_cl_dir + 'shear_cl/',
        noise_path=noise_save_dir,
        mixmats_path=mix_mats_save_path,
        #bp_cov_filemask=cov_fromsim_path + 'cov_{n_bp}bp.npz',
        bp_cov_filemask=save_combined_cov_path + 'cov_{n_bp}bp.npz',
        binmixmat_save_dir=binmixmat_save_dir,
        varied_params=['w', 'wa'],
        like_save_dir=inference_analysis_output_dir,
        obs_bandpowers_dir=save_dir + 'measured_3x2pt_bps/',
        bandpower_spacing='log',
        bandpower_edges=None,
        #cov_blocks_path=cov_fromsim_path)
        cov_blocks_path=save_combined_cov_path)

if obs_type == '3X2PT':

    covariance.bin_combine_cov(n_zbin=nbins,
                               lmin_in=output_lmin,
                               lmin_out=output_lmin,
                               lmax=output_lmax,
                               n_bp_min=n_bandpowers,
                               n_bp_max=n_bandpowers,
                               input_filemask=cov_fromsim_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                               save_filemask=cov_fromsim_path + 'cov_{n_bp}bp.npz',
                               bandpower_spacing='log',
                               obs_type='3x2pt')

    # Let's try to make an analytical covariance matrix...

    # Calculate covariance of signal part of measurement
    save_signal_cov_path = inference_analysis_output_dir + 'signal_cov/'
    if not os.path.exists(save_signal_cov_path):
        os.makedirs(save_signal_cov_path)

    covariance.get_3x2pt_cov_signal_pcl(
        n_zbin=nbins,
        pos_pos_filemask=theory_cl_dir + 'galaxy_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        pos_she_filemask=theory_cl_dir + 'galaxy_shear_cl/bin_{pos_zbin}_{she_zbin}.txt',
        she_she_filemask=theory_cl_dir + 'shear_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        lmax_out=output_lmax,
        lmin_out=output_lmin,
        mask_path=mask_dir,
        nside=nside,
        save_filemask=save_signal_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz')

    # Calculate covariance of noise part of measurement

    save_noise_cov_path = inference_analysis_output_dir + 'noise_cov/'
    if not os.path.exists(save_noise_cov_path):
        os.makedirs(save_noise_cov_path)

    covariance.get_3x2pt_cov_noise(
        n_zbin=nbins,
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        lmax_out=output_lmax,
        lmin_out=output_lmin,
        pos_nl_path=noise_save_dir + 'galaxy_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        she_nl_path=noise_save_dir + 'shear_cl/bin_{hi_zbin}_{lo_zbin}.txt',
        pos_she_nl_path=noise_save_dir + 'galaxy_shear_cl/bin_{pos_zbin}_{she_zbin}.txt',
        noise_lmin=output_lmin,
        nside=nside,
        save_filemask=save_noise_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz')

    # Step f)
    # Bin the signal covariance into bandpowers
    covariance.bin_combine_cov(n_zbin=nbins,
                               lmin_in=output_lmin,
                               lmin_out=output_lmin,
                               lmax=output_lmax,
                               n_bp_min=n_bandpowers,
                               n_bp_max=n_bandpowers,
                               input_filemask=save_signal_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                               save_filemask=save_signal_cov_path + 'cov_{n_bp}bp.npz')

    # Bin the noise covariance into bandpowers
    covariance.bin_combine_cov(n_zbin=nbins,
                               lmin_in=output_lmin,
                               lmin_out=output_lmin,
                               lmax=output_lmax,
                               n_bp_min=n_bandpowers,
                               n_bp_max=n_bandpowers,
                               input_filemask=save_noise_cov_path + 'cov_spec1_{spec1_idx}_spec2_{spec2_idx}.npz',
                               save_filemask=save_noise_cov_path + 'cov_{n_bp}bp.npz')

    save_combined_cov_path = inference_analysis_output_dir + 'combined_signal_noise_cov/'
    if not os.path.exists(save_combined_cov_path):
        os.makedirs(save_combined_cov_path)

    # Combine the individual signal and noise covariance blocks
    covariance.combine_signal_noise_cov_spec(
        n_zbin=nbins,
        signal_cov_path=save_signal_cov_path,
        noise_cov_path=save_noise_cov_path,
        combined_cov_save_path=save_combined_cov_path)

    # Combine the overall covariance matrix
    covariance.combine_signal_noise_cov(
        n_zbin=nbins,
        lmin_in=output_lmin,
        lmin_out=output_lmin,
        lmax=output_lmax,
        n_bp_min=n_bandpowers,
        n_bp_max=n_bandpowers,
        signal_cov_path=save_signal_cov_path + 'cov_{n_bp}bp.npz',
        noise_cov_path=save_noise_cov_path + 'cov_{n_bp}bp.npz',
        combined_cov_save_path=save_combined_cov_path + 'cov_{n_bp}bp.npz')

    #save_combined_cov_path = inference_analysis_output_dir + 'combined_signal_noise_cov/'

    # Execute the likelihood analysis
    loop_likelihood_nbin.like_bp_gauss_mix_loop_nbin(
        grid_dir=chains_input_dir,
        n_bps=np.array([n_bandpowers]),
        n_zbin=nbins,
        lmax_like=output_lmax,
        lmin_like=output_lmin,
        lmax_in=input_lmax,
        lmin_in=input_lmin,
        fid_pos_pos_dir=theory_cl_dir+'galaxy_cl/',
        fid_she_she_dir=theory_cl_dir+'shear_cl/',
        fid_pos_she_dir=theory_cl_dir+'galaxy_shear_cl/',
        noise_path=noise_save_dir,
        mixmats_path=mix_mats_save_path,
        #bp_cov_filemask=cov_fromsim_path+'cov_{n_bp}bp.npz',
        bp_cov_filemask=save_combined_cov_path+'cov_{n_bp}bp.npz',
        binmixmat_save_dir=binmixmat_save_dir,
        varied_params=['w', 'wa'],
        like_save_dir=inference_analysis_output_dir,
        obs_bandpowers_dir=save_dir+'measured_3x2pt_bps/',
        bandpower_spacing='log',
        bandpower_edges=None,
        #cov_blocks_path=cov_fromsim_path)
        cov_blocks_path=save_combined_cov_path)

posterior.cl_post(
    log_like_filemask=cl_like_filemask,
    contour_levels_sig=[1, 2, 3],
    bp=n_bandpowers,
    colour='C0',
    linestyle='-',
    zrange=[zmin, zmax],
    lrange=[output_lmin, output_lmax],
    ngals=ngal,
    nside=nside,
    n_bandpowers=n_bandpowers,
    obs_type=obs_type,
    plot_save_path=inference_analysis_output_dir+'contours_l{}-{}_{}.png'.format(output_lmin, output_lmax, obs_type)
)
