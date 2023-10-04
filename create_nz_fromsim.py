import os
import configparser
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.colors import LogNorm


def nz_fromsim_config(pipeline_variables_path):
    """

    Parameters
    ----------
    pipeline_variables_path (str):  Path to location of set_variables.ini file

    Returns
    -------
    Dictionary of nz parameters
    """

    config = configparser.ConfigParser()
    config.read(pipeline_variables_path)

    zmin = float(config['create_nz']['ZMIN'])
    zmax = float(config['create_nz']['ZMAX'])

    # Precision/step-size of z-range that is sampled over.
    dz = float(config['create_nz']['DZ'])

    nbins = int(float(config['create_nz']['N_ZBIN']))
    bin_type = str(config['create_nz']['ZBIN_TYPE'])

    nz_table_filename = str(config['create_nz']['MEASURED_NZ_TABLE_NAME'])

    save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
    catalogue_dir = str(config['measurement_setup']['CATALOGUE_DIR'])
    realisations = int(float(config['measurement_setup']['REALISATIONS']))

    sigma_phot = str(config['noise_cls']['SIGMA_PHOT'])
    sigma_shear = str(config['noise_cls']['SIGMA_SHEAR'])

    # Prepare config dictionary
    config_dict = {
        'zmin': zmin,
        'zmax': zmax,
        'dz': dz,
        'nbins': nbins,
        'bin_type': bin_type,
        'nz_table_filename': nz_table_filename,
        'save_dir': save_dir,
        'catalogue_dir': catalogue_dir,
        'realisations': realisations

    }

    return config_dict


def create_zbin_boundaries(config_dict):
    zmin = config_dict['zmin']
    zmax = config_dict['zmax']
    dz = config_dict['dz']
    nbins = config_dict['nbins']
    bin_type = config_dict['bin_type']
    save_dir = config_dict['save_dir']
    catalogue_dir = config_dict['catalogue_dir']

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    z_boundaries_filename = 'z_boundaries.txt'

    if bin_type == 'EQUI_Z':
        z_boundaries_low = np.linspace(zmin, zmax, nbins + 1)
        z_boundaries_mid = z_boundaries_low + (((zmax - zmin) / nbins) / 2)
        z_boundaries_high = z_boundaries_mid + (((zmax - zmin) / nbins) / 2)

        z_boundaries = [z_boundaries_low, z_boundaries_mid, z_boundaries_high]
        np.savetxt(save_dir + z_boundaries_filename,
                   np.transpose(z_boundaries),
                   fmt=['%.2f', '%.2f', '%.2f'])

    if bin_type == 'EQUI_POP':

        # Need to generate a rnd_sample from the measured n(z), i.e. n(z)*z for each z

        mock_cat_filename = 'Raw_Galaxy_Sample.hdf5'
        mock_cat = catalogue_dir + mock_cat_filename
        with h5py.File(mock_cat, "r") as f:
            rnd_sample = f['Redshift_z'][()]

        # rnd_sample = np.load(save_dir)  # Placeholder for now!

        rnd_sample = np.round(rnd_sample, 2)
        sorted_sample = np.sort(rnd_sample)
        split_sorted_sample = np.array_split(sorted_sample, nbins)
        z_boundaries_low = [zmin]
        z_boundaries_high = []
        for i in range(nbins):
            z_boundaries_low.append(split_sorted_sample[i][-1])
            z_boundaries_high.append(split_sorted_sample[i][-1])
        z_boundaries_high.append(z_boundaries_high[-1] + dz)
        z_boundaries_mid = []
        for i in range(len(z_boundaries_low)):
            z_boundaries_mid.append(round(np.mean([z_boundaries_low[i], z_boundaries_high[i]]), 2))

        z_boundaries = [z_boundaries_low, z_boundaries_mid, z_boundaries_high]
        np.savetxt(save_dir + z_boundaries_filename,
                   np.transpose(z_boundaries),
                   fmt=['%.2f', '%.2f', '%.2f'])

    if bin_type == 'EQUI_D':
        # we need to go back to the directory of the simulation and into the cosmosis/distances file for the
        # comoving distance as a function of z. Then cut out the range that corresponds to the z_range of observation
        # then define equally spaced boundaries in d-space and take the corresponding z boundaries

        z_distances = np.loadtxt(catalogue_dir + 'cosmosis/distances/z.txt')
        d_m = np.loadtxt(catalogue_dir + 'cosmosis/distances/d_m.txt')
        zmin_id = np.where((z_distances == zmin))[0][0]
        zmax_id = np.where((z_distances == zmax))[0][0]
        #print(zmin_id, zmax_id)
        d_m_observed = d_m[zmin_id:zmax_id+1]
        z_observed = z_distances[zmin_id:zmax_id+1]
        d_m_range = d_m_observed[-1]-d_m_observed[0]
        d_m_separation = d_m_range/nbins
        z_boundaries_low = [zmin]
        z_boundaries_high = []

        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx

        for i in range(nbins):
            obs_id = find_nearest(d_m_observed, d_m_observed[0] + (d_m_separation*(i+1)))
            #print(obs_id)
            z_boundaries_low.append(z_observed[obs_id])
            z_boundaries_high.append(z_observed[obs_id])
        z_boundaries_high.append(z_boundaries_high[-1] + dz)

        z_boundaries_mid = []

        for i in range(len(z_boundaries_low)):
            z_boundaries_mid.append(round(np.mean([z_boundaries_low[i], z_boundaries_high[i]]), 2))

        z_boundaries = [z_boundaries_low, z_boundaries_mid, z_boundaries_high]

        np.savetxt(save_dir + z_boundaries_filename,
                   np.transpose(z_boundaries),
                   fmt=['%.2f', '%.2f', '%.2f'])

    return np.asarray(z_boundaries)


def create_nz_fromsim(config_dict, z_boundaries, redshift_type):
    zmin = config_dict['zmin']
    zmax = config_dict['zmax']
    dz = config_dict['dz']
    nbins = config_dict['nbins']
    bin_type = config_dict['bin_type']
    save_dir = config_dict['save_dir']
    catalogue_dir = config_dict['catalogue_dir']
    realisations = config_dict['realisations']
    nz_table_filename = config_dict['nz_table_filename']
    sigma_phot = config_dict['sigma_phot']
    sigma_shear = config_dict['sigma_shear']

    z_boundaries_low = z_boundaries[0]
    z_boundaries_mid = z_boundaries[1]
    z_boundaries_high = z_boundaries[2]

    z_boundaries_low = z_boundaries_low.round(decimals=2)
    z_boundaries_mid = z_boundaries_mid.round(decimals=2)
    z_boundaries_high = z_boundaries_high.round(decimals=2)

    sub_hist_bins = np.linspace(
        zmin,
        zmax + dz,
        (round((zmax + dz - zmin) / dz)) + 1
    )

    hists = defaultdict(list)
    for b in range(nbins):
        hists["BIN_{}".format(b + 1)] = []

    for i in range(realisations):

        f = h5py.File(catalogue_dir + 'cat_products/master_cats/master_cat_poisson_sampled_{}.hdf5'.format(i + 1), 'r')
        true_z = np.array(f.get('True_Redshift_z'))

        if redshift_type == 'Observed':
            # obs_z = np.array(f.get('Redshift_z'))
            obs_z = np.array(f.get('Gaussian_Redshift_z')) # when including Catastrophic Photo-z errors
        elif redshift_type == 'True':
            obs_z = true_z
        else:
            print('Unknown Redshift Type - Please specify \'Observed\' or \'True\' for n(z) generation')
            sys.exit()

        if dz == 0.1:
            obs_z = np.around(obs_z, decimals=1)
        elif dz == 0.01:
            obs_z = np.around(obs_z, decimals=2)

        f.close()
        for b in range(nbins):
            bin_pop = true_z[np.where((obs_z >= z_boundaries_low[b]) & (obs_z < z_boundaries_high[b]))[0]]
            bin_hist = np.histogram(bin_pop, bins=int(np.rint((zmax + dz - zmin) / dz)), range=(zmin, zmax))[0]
            hists["BIN_{}".format(b + 1)].append(bin_hist)

    nz = []
    nz.append(sub_hist_bins[0:-1])

    for b in range(nbins):
        iter_hist_sample = hists["BIN_{}".format(b + 1)]
        nz.append(np.mean(np.asarray(iter_hist_sample), axis=0))

    final_cat_tab = np.asarray(nz)
    '''
    final_cat_tab = np.transpose(final_cat_tab)
    zmax_pad = np.concatenate((np.array([zmax+dz]), np.zeros(nbins)))

    final_cat_tab = np.vstack((final_cat_tab, zmax_pad))
    final_cat_tab = np.transpose(final_cat_tab)
    '''

    if sigma_phot == 0:
        if zmin != 0:
            final_cat_tab = np.transpose(final_cat_tab)
            pad_vals = int((zmin-0)/dz)
            for i in range(pad_vals):
                z_pad = np.array([zmin-((i+1)*dz)])
                pad_arr = np.concatenate((z_pad, np.zeros(nbins)))
                final_cat_tab = np.vstack((pad_arr, final_cat_tab))

            final_cat_tab = np.transpose(final_cat_tab)

    return final_cat_tab


def main():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']
    # Create 'Observed Redshift'
    config_dict = nz_fromsim_config(pipeline_variables_path=pipeline_variables_path)
    z_boundaries = create_zbin_boundaries(config_dict=config_dict)
    observed_nz = create_nz_fromsim(config_dict=config_dict, z_boundaries=z_boundaries, redshift_type='Observed')
    true_nz = create_nz_fromsim(config_dict=config_dict, z_boundaries=z_boundaries, redshift_type='True')
    save_dir = config_dict['save_dir']
    nz_table_filename = config_dict['nz_table_filename']

    np.savetxt(
        save_dir + nz_table_filename,
        np.transpose(observed_nz)
    )

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax1 = plt.subplots(figsize=(7.5, 6))

    for i in range(len(observed_nz) - 1):
        zs = observed_nz[0]

        true_nzs = true_nz[i + 1]
        true_nzs = true_nzs.astype(np.float64)
        true_nzs[true_nzs == 0] = np.nan

        obs_nzs = observed_nz[i + 1]
        obs_nzs = obs_nzs.astype(np.float64)
        obs_nzs[obs_nzs == 0] = np.nan

        ax1.plot(zs, true_nzs, marker=None, linestyle=':', markersize=5, color=colors[i])
        ax1.plot(zs, obs_nzs, 'o', markersize=5, color=colors[i])

    ax1.set_xlabel('Redshift ' r'$z$', fontsize=15, labelpad=10)
    ax1.set_ylabel(r'$n(z)$' ' [No. Galaxies/' r'$dz=0.1$' ']', fontsize=15, labelpad=10)
    ax1.tick_params(axis="both", direction="in")

    ax1.tick_params(right=True, top=True, labelright=False, labeltop=False)
    ax1.tick_params(axis='both', which='major', labelsize=13.5)

    plt.setp(ax1.xaxis.get_majorticklabels(), ha="center")

    plt.savefig(save_dir + 'nz.png')


if __name__ == '__main__':
    main()
