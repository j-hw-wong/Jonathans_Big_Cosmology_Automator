import os
import configparser
import numpy as np
import healpy as hp
import pymaster as nmt
from collections import defaultdict
import shutil
import sys

angular_binning_path = os.environ['ANGULAR_BINNING_PATH']
gaussian_cl_likelihood_path = os.environ['GAUSSIAN_CL_LIKELIHOOD_PATH']

sys.path.insert(1, angular_binning_path)
sys.path.insert(1, gaussian_cl_likelihood_path)

import gaussian_cl_likelihood
from gaussian_cl_likelihood.python import simulation


def measure_bps_config(pipeline_variables_path):
    config = configparser.ConfigParser()
    config.read(pipeline_variables_path)

    save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
    nbins = int(config['create_nz']['N_ZBIN'])
    nside = int(float(config['measurement_setup']['NSIDE']))
    no_iter = int(config['measurement_setup']['REALISATIONS'])
    mask_path = str(config['measurement_setup']['PATH_TO_MASK'])

    output_lmin = int(float(config['measurement_setup']['OUTPUT_ELL_MIN']))
    output_lmax = int(float(config['measurement_setup']['OUTPUT_ELL_MAX']))

    map_lmin = 0
    map_lmax = (3*nside)-1

    input_lmin = int(float(config['measurement_setup']['INPUT_ELL_MIN']))
    input_lmax = int(float(config['measurement_setup']['INPUT_ELL_MAX']))

    n_bandpowers = int(float(config['measurement_setup']['N_BANDPOWERS']))
    bandpower_spacing = str(config['measurement_setup']['BANDPOWER_SPACING'])
    accepted_bp_spacings = {'log', 'lin', 'custom'}

    if bandpower_spacing not in accepted_bp_spacings:
        print('Error! Bandpower Spacing Not Recognised - Exiting...')
        sys.exit()

    elif bandpower_spacing == 'log':
        bp_bin_edges = np.logspace(np.log10(output_lmin + 1e-5), np.log10(output_lmax + 1e-5), n_bandpowers + 1)

    elif bandpower_spacing == 'lin':
        bp_bin_edges = np.linspace(output_lmin + 1e-5, output_lmax + 1e-5, n_bandpowers + 1)

    elif bandpower_spacing == 'custom':
        # Placeholder at the moment!!!
        bp_bin_edges = np.empty(n_bandpowers)
        print('Need to figure out how to do something. Future JWong todo.')
        # I think it should probably be something like set_variables.ini points to a file in which the bandpower edges
        # are stored. Probably need to think about format for this. Also where to calculate this? Or just let user do
        # edges = np.loadtxt(dir + str(config['measure_3x2pt_cls']['BANDPOWER_EDGES']))
        # Then might need to transpose depending on the structure of supplied bandpower edges.
        # if os.path.isfile(dir + str(config['measure_3x2pt_cls']['BANDPOWER_EDGES'])) is False:
        #     print('WARNING - Must supply custom bandpower edges')
        #     sys.exit()

    else:
        # Bandpower type not recognised
        sys.exit()

    eff_bp_bin_edges = np.ceil(bp_bin_edges)
    eff_bp_bin_edges = eff_bp_bin_edges.astype(int)

    eff_bp_bin_edges_lo = eff_bp_bin_edges[:-1]
    eff_bp_bin_edges_hi = eff_bp_bin_edges[1:]

    bp_bins = nmt.NmtBin.from_edges(ell_ini=eff_bp_bin_edges_lo, ell_end=eff_bp_bin_edges_hi)
    ell_arr = bp_bins.get_effective_ells()

    pbl = gaussian_cl_likelihood.python.simulation.get_binning_matrix(
        n_bandpowers=n_bandpowers,
        output_lmin=output_lmin,
        output_lmax=output_lmax,
        bp_spacing=bandpower_spacing)

    # Prepare config dictionary
    config_dict = {
        'save_dir': save_dir,
        'nbins': nbins,
        'no_iter': no_iter,
        'mask_path': mask_path,
        'output_lmin': output_lmin,
        'output_lmax': output_lmax,
        'input_lmin': input_lmin,
        'input_lmax': input_lmax,
        'map_lmin': map_lmin,
        'map_lmax': map_lmax,
        'n_bandpowers': n_bandpowers,
        'bp_bins': bp_bins,
        'bandpower_spacing': bandpower_spacing,
        'ell_arr': ell_arr,
        'pbl': pbl
    }

    return config_dict


def measured_cls_to_obs_cls(measured_cls_dir, obs_cls_dir, bin_i, bin_j, lmin_out, lmax_out):

    measured_cl = np.loadtxt(measured_cls_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j))
    obs_cl = measured_cl[lmin_out:lmax_out+1]

    if not os.path.exists(obs_cls_dir):
        os.makedirs(obs_cls_dir)

    np.savetxt(obs_cls_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j), obs_cl)


def cl_to_bp(cl_dir, bp_dir, bin_i, bin_j, pbl):
    this_cl = np.loadtxt(cl_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j))
    this_bp = pbl @ this_cl

    if not os.path.exists(bp_dir):
        os.makedirs(bp_dir)

    np.savetxt(bp_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j), this_bp)


###############################################################################################
# Following bit of code is to make a prediction for the theory bandpowers based on the survey
###############################################################################################

def setup_theory_cls(cl_dir, spectra_type, bin_i, bin_j):
    """
    Opens some theory Cls based on field type and bin combination to prepare for decoupling into a binned Pseudo-Cl
    Parameters
    ----------
    cl_dir
    spectra_type
    bin_i
    bin_j

    Returns
    -------

    """

    accepted_spectra_types = {'TT', 'TE', 'TB', 'EE', 'EB', 'BE', 'BB', 'kk', 'gal_gal', 'gal_E', 'gal_B'}

    if spectra_type not in accepted_spectra_types:
        print(spectra_type)
        print('Error! Field Type Not Recognised - Exiting...')
        sys.exit()

    # Open correct cl folder - based on default CosmoSIS naming conventions

    cl_paths = {
        'TT': 'shear_cl/',
        'TE': 'null_spectra/',
        'TB': 'null_spectra/',
        'EE': 'shear_cl/',
        'EB': 'null_spectra/',
        'BE': 'null_spectra/',
        'BB': 'null_spectra/',
        'gal_gal': 'galaxy_cl/',
        'gal_E': 'galaxy_shear_cl/',
        'gal_B': 'null_spectra/'
    }

    theory_cls_dir = cl_dir + cl_paths[spectra_type]
    theory_cls = np.loadtxt(theory_cls_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j))
    return theory_cls, theory_cls_dir


def pad_cls(lmin, input_cls):
    """
    Convenience function - in case lmin!=0, we need to pad the theory Cls with zeros at 0<=l<lmin in order to
    combine with the mixing + binning matrices used in NaMaster to generate a binned theoretical Pseudo-Cl
    Parameters
    ----------
    lmin
    input_cls

    Returns
    -------

    """

    if lmin != 0:
        pad_cl_arr = np.zeros(lmin)
        input_cls = np.concatenate((pad_cl_arr, input_cls))
    else:
        pass
    return input_cls


def create_null_spectras(nbins, lmin, lmax, output_dir):
    # Create 'null spectra.' Pipeline assumes only non-zero components are TT, EE, gal_gal, gal_shear (gal_E).
    # So we need to create some 'zero' spectra of the same size as the non-zero components to do inference analysis,
    # comparison testing, systematics analysis etc.
    # This could probably be put somewhere else in the pipeline - in this way we rewrite the files every iteration.

    null_spectras_dir = output_dir + 'null_spectra/'
    if not os.path.exists(null_spectras_dir):
        os.makedirs(null_spectras_dir)

    null_ells = np.linspace(lmin, lmax + 1, (lmax - lmin) + 1)

    np.savetxt(null_spectras_dir + 'ell.txt',
               np.transpose(null_ells))

    for i in range(nbins):
        for j in range(nbins):
            null_cls = np.zeros(len(null_ells))
            np.savetxt(null_spectras_dir + 'bin_{}_{}.txt'.format(i + 1, j + 1),
                       np.transpose(null_cls))


def process_00_pcls(config_dict, theory_cl_dir, noise_cl_dir, spectra_type, bin_i, bin_j, obs_mask_path, bp_bins,
                    ell_arr, pbl):
    # Need to write some function to open the theory Cls to then decouple into a binned Pseudo-Cl calculated from each
    # of the measured ClKK, ClEE, ClBB and their cross correlations. Could save these matrices to a file and decouple
    # in a later module but it seems easier this way...

    output_lmin = config_dict['output_lmin']
    output_lmax = config_dict['output_lmax']
    input_lmin = config_dict['input_lmin']
    input_lmax = config_dict['input_lmax']
    map_lmin = config_dict['map_lmin']
    map_lmax = config_dict['map_lmax']

    accepted_spectra_types = {'TT', 'gal_gal'}
    if spectra_type not in accepted_spectra_types:
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    obs_mask = hp.read_map(obs_mask_path)

    # lmax_sht here should be input lmax
    f0 = nmt.NmtField(mask=obs_mask, maps=None, spin=0, lmax_sht=input_lmax, lite=True)

    w = nmt.NmtWorkspace()
    #w.compute_coupling_matrix(f0, f0, bins=bp_bins)
    w.compute_coupling_matrix(f0, f0, bins=nmt.NmtBin.from_lmax_linear(input_lmax, 1))

    # theory cls go from input lmin to input lmax
    theory_cls, bp_save_dir = setup_theory_cls(
        cl_dir=theory_cl_dir,
        spectra_type=spectra_type,
        bin_i=bin_i,
        bin_j=bin_j)

    noise_cls = setup_theory_cls(
        cl_dir=noise_cl_dir,
        spectra_type=spectra_type,
        bin_i=bin_i,
        bin_j=bin_j)[0]

    theory_cls = pad_cls(lmin=input_lmin, input_cls=theory_cls)
    noise_cls = pad_cls(lmin=output_lmin, input_cls=noise_cls)
    theory_cls = [theory_cls]
    theory_pcl = w.couple_cell(theory_cls)

    binned_theory_pcl = pbl @ (theory_pcl[0][output_lmin:output_lmax + 1] + noise_cls[output_lmin:output_lmax + 1])
    # binned_theory_pcl = pbl @ (theory_pcl[0][output_lmin:output_lmax + 1])

    np.savetxt(bp_save_dir + 'PCl_Bandpowers_{}_bin_{}_{}.txt'.format(spectra_type, bin_i, bin_j),
               np.transpose(binned_theory_pcl))

    if os.path.isfile(bp_save_dir + 'ell_measured.txt') is False:
        np.savetxt(bp_save_dir + 'ell_measured.txt',
                   np.transpose(ell_arr))


def process_02_pcls(config_dict, theory_cl_dir, noise_cl_dir, spectra_type, bin_i, bin_j, obs_mask_path, bp_bins,
                    ell_arr, pbl):

    output_lmin = config_dict['output_lmin']
    output_lmax = config_dict['output_lmax']
    input_lmin = config_dict['input_lmin']
    input_lmax = config_dict['input_lmax']
    map_lmin = config_dict['map_lmin']
    map_lmax = config_dict['map_lmax']

    accepted_spectra_types = {'TE', 'TB', 'gal_E', 'gal_B'}
    if spectra_type not in accepted_spectra_types:
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    obs_mask = hp.read_map(obs_mask_path)

    f0 = nmt.NmtField(mask=obs_mask, maps=None, spin=0, lmax_sht=input_lmax)
    f2 = nmt.NmtField(mask=obs_mask, maps=None, spin=2, lmax_sht=input_lmax)

    w = nmt.NmtWorkspace()
    #w.compute_coupling_matrix(f0, f2, bins=bp_bins)
    w.compute_coupling_matrix(f0, f2, bins=nmt.NmtBin.from_lmax_linear(input_lmax, 1))

    binned_theory_pcls = defaultdict(list)

    if spectra_type == 'TE' or 'TB':
        theory_TE = pad_cls(
            lmin=input_lmin,
            input_cls=setup_theory_cls(
                cl_dir=theory_cl_dir,
                spectra_type='TE',
                bin_i=bin_i,
                bin_j=bin_j)[0])

        theory_TB = pad_cls(
            lmin=input_lmin,
            input_cls=setup_theory_cls(
                cl_dir=theory_cl_dir,
                spectra_type='TB',
                bin_i=bin_i,
                bin_j=bin_j)[0])

        noise_cls = pad_cls(
            lmin=output_lmin,
            input_cls=setup_theory_cls(
                cl_dir=noise_cl_dir,
                spectra_type=spectra_type,
                bin_i=bin_i,
                bin_j=bin_j)[0])

        theory_cls = [theory_TE, theory_TB]
        theory_pcl = w.couple_cell(theory_cls)

        binned_theory_pcls['TE'] = pbl @ (theory_pcl[0][output_lmin:output_lmax + 1]
                                          + noise_cls[output_lmin:output_lmax + 1])
        binned_theory_pcls['TB'] = pbl @ (theory_pcl[1][output_lmin:output_lmax + 1]
                                          + noise_cls[output_lmin:output_lmax + 1])

    if spectra_type == 'gal_E' or 'gal_B':
        theory_gal_E = pad_cls(
            lmin=input_lmin,
            input_cls=setup_theory_cls(
                cl_dir=theory_cl_dir,
                spectra_type='gal_E',
                bin_i=bin_i,
                bin_j=bin_j)[0])

        theory_gal_B = pad_cls(
            lmin=input_lmin,
            input_cls=setup_theory_cls(
                cl_dir=theory_cl_dir,
                spectra_type='gal_B',
                bin_i=bin_i,
                bin_j=bin_j)[0])

        noise_cls = pad_cls(
            lmin=output_lmin,
            input_cls=setup_theory_cls(
                cl_dir=noise_cl_dir,
                spectra_type=spectra_type,
                bin_i=bin_i,
                bin_j=bin_j)[0])

        theory_cls = [theory_gal_E, theory_gal_B]
        theory_pcl = w.couple_cell(theory_cls)

        binned_theory_pcls['gal_E'] = pbl @ (theory_pcl[0][output_lmin:output_lmax + 1]
                                             + noise_cls[output_lmin:output_lmax + 1])
        #binned_theory_pcls['gal_E'] = pbl @ (theory_pcl[0][output_lmin:output_lmax + 1])
        binned_theory_pcls['gal_B'] = pbl @ (theory_pcl[1][output_lmin:output_lmax + 1]
                                             + noise_cls[output_lmin:output_lmax + 1])

    bp_save_dir = setup_theory_cls(cl_dir=theory_cl_dir, spectra_type=spectra_type, bin_i=bin_i, bin_j=bin_j)[1]

    np.savetxt(bp_save_dir + 'PCl_Bandpowers_{}_bin_{}_{}.txt'.format(spectra_type, bin_i, bin_j),
               np.transpose(binned_theory_pcls[spectra_type]))

    if os.path.isfile(bp_save_dir + 'ell_measured.txt') is False:
        np.savetxt(bp_save_dir + 'ell_measured.txt',
                   np.transpose(ell_arr))


def process_22_pcls(config_dict, theory_cl_dir, noise_cl_dir, spectra_type, bin_i, bin_j, obs_mask_path, bp_bins,
                    ell_arr, pbl):

    output_lmin = config_dict['output_lmin']
    output_lmax = config_dict['output_lmax']
    input_lmin = config_dict['input_lmin']
    input_lmax = config_dict['input_lmax']
    map_lmin = config_dict['map_lmin']
    map_lmax = config_dict['map_lmax']

    accepted_spectra_types = {'EE', 'EB', 'BE', 'BB'}
    if spectra_type not in accepted_spectra_types:
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    obs_mask = hp.read_map(obs_mask_path)

    f2 = nmt.NmtField(mask=obs_mask, maps=None, spin=2, lmax_sht=input_lmax)

    w = nmt.NmtWorkspace()
    #w.compute_coupling_matrix(f2, f2, bins=bp_bins)
    w.compute_coupling_matrix(f2, f2, bins=nmt.NmtBin.from_lmax_linear(input_lmax, 1))

    theory_EE, theory_EB, theory_BE, theory_BB = [pad_cls(
        lmin=input_lmin,
        input_cls=setup_theory_cls(
            cl_dir=theory_cl_dir,
            spectra_type=spec_type,
            bin_i=bin_i,
            bin_j=bin_j)[0]) for spec_type in ['EE', 'EB', 'BE', 'BB']]

    noise_cls = pad_cls(
        lmin=output_lmin,
        input_cls=setup_theory_cls(
            cl_dir=noise_cl_dir,
            spectra_type=spectra_type,
            bin_i=bin_i,
            bin_j=bin_j)[0])

    theory_cls = [theory_EE, theory_EB, theory_BE, theory_BB]
    theory_pcls = w.couple_cell(theory_cls)

    binned_theory_pcls = {
        'EE': pbl @ (theory_pcls[0][output_lmin:output_lmax + 1] + noise_cls[output_lmin:output_lmax + 1]),
        'EB': pbl @ (theory_pcls[1][output_lmin:output_lmax + 1] + noise_cls[output_lmin:output_lmax + 1]),
        'BE': pbl @ (theory_pcls[2][output_lmin:output_lmax + 1] + noise_cls[output_lmin:output_lmax + 1]),
        'BB': pbl @ (theory_pcls[3][output_lmin:output_lmax + 1] + noise_cls[output_lmin:output_lmax + 1])
    }

    #binned_theory_pcls = {
    #    'EE': pbl @ (theory_pcls[0][output_lmin:output_lmax + 1]),
    #    'EB': pbl @ (theory_pcls[1][output_lmin:output_lmax + 1]),
    #    'BE': pbl @ (theory_pcls[2][output_lmin:output_lmax + 1]),
    #    'BB': pbl @ (theory_pcls[3][output_lmin:output_lmax + 1])
    #}

    bp_save_dir = setup_theory_cls(cl_dir=theory_cl_dir, spectra_type=spectra_type, bin_i=bin_i, bin_j=bin_j)[1]

    np.savetxt(bp_save_dir + 'PCl_Bandpowers_{}_bin_{}_{}.txt'.format(spectra_type, bin_i, bin_j),
               np.transpose(binned_theory_pcls[spectra_type]))

    if os.path.isfile(bp_save_dir + 'ell_measured.txt') is False:
        np.savetxt(bp_save_dir + 'ell_measured.txt',
                   np.transpose(ell_arr))


def main():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']
    config_dict = measure_bps_config(pipeline_variables_path=pipeline_variables_path)

    save_dir = config_dict['save_dir']
    nbins = config_dict['nbins']
    no_iter = config_dict['no_iter']
    mask_path = config_dict['mask_path']

    bp_bins = config_dict['bp_bins']
    ell_arr = config_dict['ell_arr']
    pbl = config_dict['pbl']

    output_lmin = config_dict['output_lmin']
    output_lmax = config_dict['output_lmax']
    input_lmin = config_dict['input_lmin']
    input_lmax = config_dict['input_lmax']
    map_lmin = config_dict['map_lmin']
    map_lmax = config_dict['map_lmax']

    recov_cat_cls_dir = save_dir + 'raw_3x2pt_cls/'
    obs_cat_cls_dir = save_dir + 'measured_3x2pt_cls/'
    obs_cat_bps_dir = save_dir + 'measured_3x2pt_bps/'
    obs_noise_cls_dir = save_dir + 'measured_noise_cls/'

    gal_bps_dir = obs_cat_bps_dir + 'galaxy_bp/'
    shear_bps_dir = obs_cat_bps_dir + 'shear_bp/'
    gal_shear_bps_dir = obs_cat_bps_dir + 'galaxy_shear_bp/'

    k_bps_dir = shear_bps_dir + 'Cl_TT/'
    y1_bps_dir = shear_bps_dir + 'Cl_EE/'
    y1y2_bps_dir = shear_bps_dir + 'Cl_EB/'
    y2y1_bps_dir = shear_bps_dir + 'Cl_BE/'
    y2_bps_dir = shear_bps_dir + 'Cl_BB/'

    for folder in [gal_bps_dir, shear_bps_dir, gal_shear_bps_dir, k_bps_dir, y1_bps_dir, y1y2_bps_dir, y2y1_bps_dir,
                   y2_bps_dir]:
        if not os.path.exists(folder):
            os.makedirs(folder)
        np.savetxt(folder + 'ell_measured.txt',
                   np.transpose(ell_arr))

    gal_cl_dir = recov_cat_cls_dir + 'galaxy_cl/'
    shear_cl_dir = recov_cat_cls_dir + 'shear_cl/'
    gal_shear_cl_dir = recov_cat_cls_dir + 'galaxy_shear_cl/'

    obs_gal_cl_dir = obs_cat_cls_dir + 'galaxy_cl/'
    obs_shear_cl_dir = obs_cat_cls_dir + 'shear_cl/'
    obs_gal_shear_cl_dir = obs_cat_cls_dir + 'galaxy_shear_cl/'

    theory_cls_dir = save_dir + 'theory_cls/'
    noise_cls_dir = save_dir + 'raw_noise_cls/'

    create_null_spectras(nbins=nbins, lmin=input_lmin, lmax=input_lmax, output_dir=theory_cls_dir)
    create_null_spectras(nbins=nbins, lmin=map_lmin, lmax=map_lmax, output_dir=noise_cls_dir)

    shutil.copytree(save_dir + 'cosmosis/shear_cl', theory_cls_dir + 'shear_cl')
    shutil.copytree(save_dir + 'cosmosis/galaxy_cl', theory_cls_dir + 'galaxy_cl')
    shutil.copytree(save_dir + 'cosmosis/galaxy_shear_cl', theory_cls_dir + 'galaxy_shear_cl')

    # Get ready for some magic.. Lol.
    for i in range(nbins):
        for j in range(nbins):

            # Cut the raw PCls measured from Catalogues to a specific range in lmin, lmax
            measured_cls_to_obs_cls(
                measured_cls_dir=gal_shear_cl_dir,
                obs_cls_dir=obs_gal_shear_cl_dir,
                bin_i=i+1,
                bin_j=j+1,
                lmin_out=output_lmin,
                lmax_out=output_lmax)

            measured_cls_to_obs_cls(
                measured_cls_dir=noise_cls_dir + 'galaxy_shear_cl/',
                obs_cls_dir=obs_noise_cls_dir + 'galaxy_shear_cl/',
                bin_i=i+1,
                bin_j=j+1,
                lmin_out=output_lmin,
                lmax_out=output_lmax)

            # Convert Pseudo-Cl for measured GGL into bandpower
            cl_to_bp(cl_dir=obs_gal_shear_cl_dir, bp_dir=gal_shear_bps_dir, bin_i=i + 1, bin_j=j + 1, pbl=pbl)

            # Calculate the theoretical PCl bandpower for the fiducial full-sky GGL
            process_02_pcls(
                config_dict=config_dict,
                theory_cl_dir=theory_cls_dir,
                noise_cl_dir=obs_noise_cls_dir,
                spectra_type='gal_E',
                bin_i=i + 1,
                bin_j=j + 1,
                obs_mask_path=mask_path,
                bp_bins=bp_bins,
                ell_arr=ell_arr,
                pbl=pbl)

            if i >= j:

                measured_cls_to_obs_cls(
                    measured_cls_dir=gal_cl_dir,
                    obs_cls_dir=obs_gal_cl_dir,
                    bin_i=i+1, 
                    bin_j=j+1, 
                    lmin_out=output_lmin,
                    lmax_out=output_lmax)

                measured_cls_to_obs_cls(
                    measured_cls_dir=noise_cls_dir + 'galaxy_cl/',
                    obs_cls_dir=obs_noise_cls_dir + 'galaxy_cl/',
                    bin_i=i+1,
                    bin_j=j+1,
                    lmin_out=output_lmin,
                    lmax_out=output_lmax)

                cl_to_bp(cl_dir=obs_gal_cl_dir, bp_dir=gal_bps_dir, bin_i=i + 1, bin_j=j + 1, pbl=pbl)

                process_00_pcls(
                    config_dict=config_dict,
                    theory_cl_dir=theory_cls_dir,
                    noise_cl_dir=noise_cls_dir,
                    spectra_type='gal_gal',
                    bin_i=i + 1,
                    bin_j=j + 1,
                    obs_mask_path=mask_path,
                    bp_bins=bp_bins,
                    ell_arr=ell_arr,
                    pbl=pbl)

                measured_cls_to_obs_cls(
                    measured_cls_dir=noise_cls_dir + 'shear_cl/',
                    obs_cls_dir=obs_noise_cls_dir + 'shear_cl/',
                    bin_i=i + 1,
                    bin_j=j + 1,
                    lmin_out=output_lmin,
                    lmax_out=output_lmax)

                for shear_component_dir in ['Cl_TT/', 'Cl_EE/', 'Cl_EB/', 'Cl_BE/', 'Cl_BB/']:

                    measured_cls_to_obs_cls(
                        measured_cls_dir=shear_cl_dir + shear_component_dir,
                        obs_cls_dir=obs_shear_cl_dir + shear_component_dir,
                        bin_i=i+1, 
                        bin_j=j+1, 
                        lmin_out=output_lmin,
                        lmax_out=output_lmax)

                    cl_to_bp(
                        cl_dir=obs_shear_cl_dir + shear_component_dir,
                        bp_dir=shear_bps_dir + shear_component_dir,
                        bin_i=i + 1,
                        bin_j=j + 1,
                        pbl=pbl)

                for spin2_component in ['EE', 'BB']:
                    
                    process_22_pcls(
                        config_dict=config_dict,
                        theory_cl_dir=theory_cls_dir,
                        noise_cl_dir=noise_cls_dir,
                        spectra_type=spin2_component,
                        bin_i=i + 1,
                        bin_j=j + 1,
                        obs_mask_path=mask_path,
                        bp_bins=bp_bins,
                        ell_arr=ell_arr,
                        pbl=pbl)

    for it in range(no_iter):

        gal_cl_it_dir = gal_cl_dir + 'iter_{}/'.format(it + 1)
        gal_bp_it_dir = gal_bps_dir + 'iter_{}/'.format(it + 1)

        gal_shear_cl_it_dir = gal_shear_cl_dir + 'iter_{}/'.format(it + 1)
        gal_shear_bp_it_dir = gal_shear_bps_dir + 'iter_{}/'.format(it + 1)

        for i in range(nbins):
            for j in range(nbins):

                measured_cls_to_obs_cls(
                    measured_cls_dir=gal_shear_cl_it_dir,
                    obs_cls_dir=obs_gal_shear_cl_dir + 'iter_{}/'.format(it + 1),
                    bin_i=i + 1,
                    bin_j=j + 1,
                    lmin_out=output_lmin,
                    lmax_out=output_lmax)

                cl_to_bp(cl_dir=obs_gal_shear_cl_dir + 'iter_{}/'.format(it + 1),
                         bp_dir=gal_shear_bp_it_dir,
                         bin_i=i + 1,
                         bin_j=j + 1,
                         pbl=pbl)
                if i >= j:

                    measured_cls_to_obs_cls(
                        measured_cls_dir=gal_cl_it_dir,
                        obs_cls_dir=obs_gal_cl_dir + 'iter_{}/'.format(it + 1),
                        bin_i=i + 1,
                        bin_j=j + 1,
                        lmin_out=output_lmin,
                        lmax_out=output_lmax)

                    cl_to_bp(cl_dir=obs_gal_cl_dir + 'iter_{}/'.format(it + 1),
                             bp_dir=gal_bp_it_dir,
                             bin_i=i + 1,
                             bin_j=j + 1,
                             pbl=pbl)

                    for shear_component_dir in ['Cl_TT/', 'Cl_EE/', 'Cl_EB/', 'Cl_BE/', 'Cl_BB/']:
                        shear_cl_it_dir = shear_cl_dir + shear_component_dir + 'iter_{}/'.format(it + 1)
                        shear_bp_it_dir = shear_bps_dir + shear_component_dir + 'iter_{}/'.format(it + 1)

                        measured_cls_to_obs_cls(
                            measured_cls_dir=shear_cl_it_dir,
                            obs_cls_dir=obs_shear_cl_dir + shear_component_dir + 'iter_{}/'.format(it + 1),
                            bin_i=i + 1,
                            bin_j=j + 1,
                            lmin_out=output_lmin,
                            lmax_out=output_lmax)

                        cl_to_bp(cl_dir=obs_shear_cl_dir + shear_component_dir + 'iter_{}/'.format(it + 1),
                                 bp_dir=shear_bp_it_dir,
                                 bin_i=i + 1,
                                 bin_j=j + 1,
                                 pbl=pbl)


if __name__ == '__main__':
    main()
