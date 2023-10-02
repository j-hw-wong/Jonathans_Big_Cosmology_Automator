import os
import linecache
import statistics
import configparser
import numpy as np
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

    # Prepare config dictionary
    config_dict = {
        'nside': nside,
        'nbins': nbins,
        'pcl_lmin_out': pcl_lmin_out,
        'pcl_lmax_out': pcl_lmax_out,
        'save_dir': save_dir,
        'realisations': realisations,
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


def main():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']
    config_dict = av_cls_config(pipeline_variables_path=pipeline_variables_path)

    save_dir = config_dict['save_dir']
    nbins = config_dict['nbins']

    pcl_lmin_out = config_dict['pcl_lmin_out']
    pcl_lmax_out = config_dict['pcl_lmax_out']

    realisations = config_dict['realisations']

    noise_cls_dir = save_dir + 'raw_noise_cls/'
    measured_cls_dir = save_dir + 'raw_3x2pt_cls/'

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
