import os
import linecache
import statistics
import configparser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem


def av_cls_config(pipeline_variables_path):
    config = configparser.ConfigParser()
    config.read(pipeline_variables_path)

    save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
    nside = int(config['measurement_setup']['NSIDE'])
    realisations = int(config['measurement_setup']['realisations'])
    pcl_lmin_out = 0
    pcl_lmax_out = (3*nside)-1

    nbins = int(config['create_nz']['N_ZBIN'])

    nz_table_filename = str(config['create_nz']['MEASURED_NZ_TABLE_NAME'])

    # Prepare config dictionary
    config_dict = {
        'nside': nside,
        'nbins': nbins,
        'pcl_lmin_out': pcl_lmin_out,
        'pcl_lmax_out': pcl_lmax_out,
        'save_dir': save_dir,
        'realisations': realisations,
        'nz_table_filename': nz_table_filename
    }

    return config_dict


def calc_av_cls(cl_dir, ell_min, ell_max, bin_i, bin_j, realisations):
    cls = []
    ell = np.arange(ell_min, ell_max + 1)
    for x in range(len(ell)):
        cl_av = []
        for y in range(realisations):
            cl_file = cl_dir + 'iter_{}/bin_{}_{}.txt'.format(y + 1, bin_i, bin_j)
            a = linecache.getline(
                cl_file,
                x + 1).split()
            cl_av.append(float(a[0]))

        cls.append(statistics.mean(cl_av))

    np.savetxt(cl_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j),
               np.transpose(cls))


def calc_stdem_cls(cl_dir, ell_min, ell_max, bin_i, bin_j, realisations):
    cls_err = []
    ell = np.arange(ell_min, ell_max + 1)

    for x in range(len(ell)):

        cl_av = []

        for y in range(realisations):
            cl_file = cl_dir + 'iter_{}/bin_{}_{}.txt'.format(y + 1, bin_i, bin_j)
            a = linecache.getline(
                cl_file,
                x + 1).split()
            cl_av.append(float(a[0]))

        cls_err.append(sem(cl_av))

    np.savetxt(cl_dir + 'bin_%s_%s_err.txt' % (bin_i, bin_j),
               np.transpose(cls_err))


def calc_av_nz(nz_tables_dir, realisations):

    """
    Function to generate the n(z) for a given tomographic analysis as measured from the simulated catalogues

    Parameters
    ----------
    nz_tables_dir
    realisations

    Returns
    -------

    """

    nz_dat = []
    for i in range(realisations):

        nz_dat.append(np.transpose(np.loadtxt(nz_tables_dir + 'nz_iter{}.txt'.format(i+1))))

    return np.mean(np.asarray(nz_dat), axis=0)


def main():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']
    config_dict = av_cls_config(pipeline_variables_path=pipeline_variables_path)

    save_dir = config_dict['save_dir']
    nbins = config_dict['nbins']

    pcl_lmin_out = config_dict['pcl_lmin_out']
    pcl_lmax_out = config_dict['pcl_lmax_out']

    realisations = config_dict['realisations']

    nz_table_filename = config_dict['nz_table_filename']

    noise_cls_dir = save_dir + 'raw_noise_cls/'
    measured_cls_dir = save_dir + 'raw_3x2pt_cls/'

    final_nz_table = calc_av_nz(nz_tables_dir=save_dir+'nz_tables/', realisations=realisations)
    np.savetxt(save_dir + nz_table_filename, np.transpose(final_nz_table))

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax1 = plt.subplots(figsize=(7.5, 6))

    for i in range(len(final_nz_table) - 1):
        zs = final_nz_table[0]

        obs_nzs = final_nz_table[i + 1]
        obs_nzs = obs_nzs.astype(np.float64)
        obs_nzs[obs_nzs == 0] = np.nan

        ax1.plot(zs, obs_nzs, 'o', linestyle='--',markersize=5, color=colors[i])

    ax1.set_xlabel('Redshift ' r'$z$', fontsize=15, labelpad=10)
    ax1.set_ylabel(r'$n(z)$' ' [No. Galaxies/' r'$dz=0.1$' ']', fontsize=15, labelpad=10)
    ax1.tick_params(axis="both", direction="in")

    ax1.tick_params(right=True, top=True, labelright=False, labeltop=False)
    ax1.tick_params(axis='both', which='major', labelsize=13.5)

    plt.setp(ax1.xaxis.get_majorticklabels(), ha="center")

    plt.savefig(save_dir + 'nz.png')

    for i in range(nbins):
        for j in range(nbins):

            calc_av_cls(cl_dir=noise_cls_dir + 'galaxy_shear_cl/',
                        ell_min=pcl_lmin_out,
                        ell_max=pcl_lmax_out,
                        bin_i=i + 1,
                        bin_j=j + 1,
                        realisations=realisations)

            calc_av_cls(cl_dir=measured_cls_dir + 'galaxy_shear_cl/',
                        ell_min=pcl_lmin_out,
                        ell_max=pcl_lmax_out,
                        bin_i=i + 1,
                        bin_j=j + 1,
                        realisations=realisations)

            calc_stdem_cls(cl_dir=measured_cls_dir + 'galaxy_shear_cl/',
                           ell_min=pcl_lmin_out,
                           ell_max=pcl_lmax_out,
                           bin_i=i + 1,
                           bin_j=j + 1,
                           realisations=realisations)
            if i >= j:
                calc_av_cls(noise_cls_dir + 'galaxy_cl/',
                            ell_min=pcl_lmin_out,
                            ell_max=pcl_lmax_out,
                            bin_i=i + 1,
                            bin_j=j + 1,
                            realisations=realisations)

                calc_av_cls(measured_cls_dir + 'galaxy_cl/',
                            ell_min=pcl_lmin_out,
                            ell_max=pcl_lmax_out,
                            bin_i=i + 1,
                            bin_j=j + 1,
                            realisations=realisations)

                calc_stdem_cls(measured_cls_dir + 'galaxy_cl/',
                               ell_min=pcl_lmin_out,
                               ell_max=pcl_lmax_out,
                               bin_i=i + 1,
                               bin_j=j + 1,
                               realisations=realisations)

                calc_av_cls(noise_cls_dir + 'shear_cl/',
                            ell_min=pcl_lmin_out,
                            ell_max=pcl_lmax_out,
                            bin_i=i + 1,
                            bin_j=j + 1,
                            realisations=realisations)

    cl_shear_types = ['Cl_TT', 'Cl_EE', 'Cl_EB', 'Cl_BE', 'Cl_BB']

    for shear_type in cl_shear_types:
        for i in range(nbins):
            for j in range(nbins):
                if i >= j:
                    calc_av_cls(measured_cls_dir + 'shear_cl/' + shear_type + '/',
                                ell_min=pcl_lmin_out,
                                ell_max=pcl_lmax_out,
                                bin_i=i + 1,
                                bin_j=j + 1,
                                realisations=realisations)


if __name__ == '__main__':
    main()
