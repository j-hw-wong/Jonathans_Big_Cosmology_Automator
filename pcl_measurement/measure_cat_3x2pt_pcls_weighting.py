import os
import sys
import h5py
import configparser
import numpy as np
import healpy as hp
import pymaster as nmt
import matplotlib.pyplot as plt
from collections import defaultdict


def measure_pcls_config(pipeline_variables_path):
    config = configparser.ConfigParser()
    config.read(pipeline_variables_path)

    save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
    cats_dir = str(config['measurement_setup']['CATALOGUE_DIR'])
    realisations = int(config['measurement_setup']['REALISATIONS'])

    zmin = float(config['create_nz']['ZMIN'])
    zmax = float(config['create_nz']['ZMAX'])

    # Precision/step-size of z-range that is sampled over.
    dz = float(config['create_nz']['DZ'])

    nside = int(config['measurement_setup']['NSIDE'])
    nbins = int(config['create_nz']['N_ZBIN'])
    npix = hp.nside2npix(nside)

    # These lmin, lmax out could be changed - lrange that is measured out from Healpix maps
    raw_pcl_lmin_out = 0
    raw_pcl_lmax_out = (3*nside)-1

    sigma_phot = float(config['noise_cls']['SIGMA_PHOT'])
    sigma_e = float(config['noise_cls']['SIGMA_SHEAR'])

    mask_path = str(config['measurement_setup']['PATH_TO_MASK'])

    # Prepare config dictionary
    config_dict = {
        'nside': nside,
        'nbins': nbins,
        'npix': npix,
        'zmin': zmin,
        'zmax': zmax,
        'dz': dz,
        'raw_pcl_lmin_out': raw_pcl_lmin_out,
        'raw_pcl_lmax_out': raw_pcl_lmax_out,
        'save_dir': save_dir,
        'cats_dir': cats_dir,
        'realisations': realisations,
        'sigma_phot': sigma_phot,
        'sigma_e': sigma_e,
        'mask_path': mask_path
    }

    return config_dict


def maps_from_cats(config_dict, iter_no):
    nside = config_dict['nside']
    nbins = config_dict['nbins']
    npix = config_dict['npix']
    zmin = config_dict['zmin']
    zmax = config_dict['zmax']
    dz = config_dict['dz']
    lmin_out = config_dict['raw_pcl_lmin_out']
    lmax_out = config_dict['raw_pcl_lmax_out']
    save_dir = config_dict['save_dir']
    cats_dir = config_dict['cats_dir']
    sigma_phot = config_dict['sigma_phot']
    sigma_e = config_dict['sigma_e']
    mask_path = config_dict['mask_path']

    master_cat_file = cats_dir + 'cat_products/master_cats/master_cat_poisson_sampled_{}.hdf5'.format(iter_no)

    f = h5py.File(master_cat_file, 'r')
    ras = np.array(f.get('RA'))
    decs = np.array(f.get('Dec'))
    zs = np.array(f.get('Redshift_z'))
    conv = np.array(f.get('Convergence'))
    shear1 = np.array(f.get('Shear_y1'))
    shear2 = np.array(f.get('Shear_y2'))
    indices = np.array(f.get('Healpix_Index_(Position)'))
    weightings = np.array(f.get('Galaxy_Weighting'))

    # Extract all the redshift information to make the n(z)
    true_z = np.array(f.get('True_Redshift_z'))
    obs_gaussian_zs = np.array(f.get('Gaussian_Redshift_z'))

    f.close()

    ras = np.asarray(ras)
    decs = np.asarray(decs)
    zs = np.asarray(zs)
    conv = np.asarray(conv)
    shear1 = np.asarray(shear1)
    shear2 = np.asarray(shear2)
    indices = np.asarray(indices)
    indices = indices.astype(int)

    zs = np.around(zs, decimals=1)

    z_boundaries_filename = 'z_boundaries.txt'
    z_boundaries = np.loadtxt(save_dir + z_boundaries_filename)
    z_boundary_columns = np.transpose(z_boundaries)
    # could import and run function (try?)

    z_boundaries_low = z_boundary_columns[0][0:-1]
    z_boundaries_mid = z_boundary_columns[1][0:-1]
    z_boundaries_high = z_boundary_columns[2][0:-1]

    z_boundaries_low = np.round(z_boundaries_low, 2)
    z_boundaries_mid = np.round(z_boundaries_mid, 2)
    z_boundaries_high = np.round(z_boundaries_high, 2)

    # Let's make the n(z) here - useful for making theoretical predictions and inference analysis

    sub_hist_bins = np.linspace(
        zmin,
        zmax + dz,
        (round((zmax + dz - zmin) / dz)) + 1
    )

    hists = defaultdict(list)
    for b in range(nbins):
        hists["BIN_{}".format(b + 1)] = []

    if dz == 0.1:
        obs_gaussian_zs = np.around(obs_gaussian_zs, decimals=1)
    elif dz == 0.01:
        obs_gaussian_zs = np.around(obs_gaussian_zs, decimals=2)
    '''
    for b in range(nbins):
        bin_pop = true_z[np.where((obs_gaussian_zs >= z_boundaries_low[b]) & (obs_gaussian_zs < z_boundaries_high[b]))[0]]
        bin_hist = np.histogram(bin_pop, bins=int(np.rint((zmax + dz - zmin) / dz)), range=(zmin, zmax))[0]
        hists["BIN_{}".format(b + 1)].append(bin_hist)

    nz = []
    nz.append(sub_hist_bins[0:-1])

    for b in range(nbins):
        iter_hist_sample = hists["BIN_{}".format(b + 1)][0]
        nz.append(iter_hist_sample)

    final_cat_tab = np.asarray(nz)

    if sigma_phot == 0:
        if zmin != 0:
            final_cat_tab = np.transpose(final_cat_tab)
            pad_vals = int((zmin-0)/dz)
            for i in range(pad_vals):
                z_pad = np.array([zmin-((i+1)*dz)])
                pad_arr = np.concatenate((z_pad, np.zeros(nbins)))
                final_cat_tab = np.vstack((pad_arr, final_cat_tab))

            final_cat_tab = np.transpose(final_cat_tab)

    nz_save_dir = save_dir + 'nz_tables/'
    if not os.path.exists(nz_save_dir):
        os.makedirs(nz_save_dir)

    np.savetxt(nz_save_dir + 'nz_iter{}.txt'.format(iter_no),
               np.transpose(final_cat_tab))
    '''
    # Let's make the 3x2pt maps now

    ell_arr = np.arange(lmin_out, lmax_out + 1, 1)

    recov_cat_cls_dir = save_dir + 'raw_3x2pt_cls/'

    gal_dir = recov_cat_cls_dir + 'galaxy_cl/iter_{}/'.format(iter_no)
    if not os.path.exists(gal_dir):
        os.makedirs(gal_dir)

    shear_dir = recov_cat_cls_dir + 'shear_cl/'.format(iter_no)
    if not os.path.exists(shear_dir):
        os.makedirs(shear_dir)

    gal_shear_dir = recov_cat_cls_dir + 'galaxy_shear_cl/iter_{}/'.format(iter_no)
    if not os.path.exists(gal_shear_dir):
        os.makedirs(gal_shear_dir)

    k_dir = shear_dir + 'Cl_TT/iter_{}/'.format(iter_no)
    y1_dir = shear_dir + 'Cl_EE/iter_{}/'.format(iter_no)
    y1y2_dir = shear_dir + 'Cl_EB/iter_{}/'.format(iter_no)
    y2y1_dir = shear_dir + 'Cl_BE/iter_{}/'.format(iter_no)
    y2_dir = shear_dir + 'Cl_BB/iter_{}/'.format(iter_no)

    for folder in [gal_dir, gal_shear_dir, k_dir, y1_dir, y1y2_dir, y2y1_dir, y2_dir]:
        if not os.path.exists(folder):
            os.makedirs(folder)
        np.savetxt(folder + 'ell.txt',
                   np.transpose(ell_arr))

    poisson_cls_theory_dir = save_dir + 'raw_noise_cls/galaxy_cl/iter_{}/'.format(iter_no)
    if not os.path.exists(poisson_cls_theory_dir):
        os.makedirs(poisson_cls_theory_dir)

    shape_noise_cls_dir = save_dir + 'raw_noise_cls/shear_cl/iter_{}/'.format(iter_no)
    if not os.path.exists(shape_noise_cls_dir):
        os.makedirs(shape_noise_cls_dir)

    print('Dividing mock catalogue into %s bins and generating weak lensing Healpix Maps' % nbins)

    cut_maps_dic = defaultdict(list)

    # obs_mask = hp.read_map(mask_path + "Euclid_DR1_{}.fits".format(nside))
    obs_mask = hp.read_map(mask_path)
    w_survey = np.count_nonzero(obs_mask == 1) / hp.nside2npix(nside)

    unobserved_pixels = np.count_nonzero(obs_mask == 0)
    observed_pixels_ids = np.where((obs_mask==1))[0]
    unobserved_pixels_ids = np.where((obs_mask==0))[0]

    nz = []
    nz.append(sub_hist_bins)

    # create nz
    
    for i in range(len(z_boundaries_low)):
        nz_inbin = np.zeros(len(sub_hist_bins))

        # this is for all stuff within a bin. We then need to do the number counts map for each redshift slice within the bin
        z_slices_inbin = np.arange(round(z_boundaries_low[i], 2), round(z_boundaries_high[i], 2), dz)

        for z_slice in z_slices_inbin:
            inds_slice = np.where((obs_gaussian_zs == round(z_slice, 2)))
            inds_slice = np.asarray(inds_slice[0])

            weightings_slice = weightings[inds_slice]

            true_zs_at_slice = true_z[inds_slice]
            for z_count, z_i in enumerate(sub_hist_bins):
                ids_z_i = np.where((true_zs_at_slice == z_i))
                nz_i = sum(weightings_slice[ids_z_i])
                nz_inbin[z_count] = nz_inbin[z_count] + nz_i

        nz.append(nz_inbin)
    
    final_cat_tab = np.asarray(nz)

    if zmin != 0:
        final_cat_tab = np.transpose(final_cat_tab)
        pad_vals = int((zmin - 0) / dz)
        for i in range(pad_vals):
            z_pad = np.array([zmin - ((i + 1) * dz)])
            pad_arr = np.concatenate((z_pad, np.zeros(nbins)))
            final_cat_tab = np.vstack((pad_arr, final_cat_tab))

        final_cat_tab = np.transpose(final_cat_tab)
    
    nz_save_dir = save_dir + 'nz_tables/'
    if not os.path.exists(nz_save_dir):
        os.makedirs(nz_save_dir)

    np.savetxt(nz_save_dir + 'nz_iter{}.txt'.format(iter_no),
               np.transpose(final_cat_tab))


    for i in range(len(z_boundaries_low)):
        # Calculate indices in catalog array of gals. that occur in each successive bin
        inds = np.where((zs >= round(z_boundaries_low[i], 2)) & (zs < z_boundaries_high[i]))
        inds = np.asarray(inds[0])

        # this is for all stuff within a bin. We then need to do the number counts map for each redshift slice within the bin
        z_slices_inbin = np.arange(round(z_boundaries_low[i], 2), round(z_boundaries_high[i], 2), dz)

        ngals_inbin = []
        k_maps_inbin = []
        y1_maps_inbin = []
        y2_maps_inbin = []
        numb_dens_maps_inbin = []

        for z_slice in z_slices_inbin:

            inds_slice = np.where((obs_gaussian_zs == round(z_slice, 2)))
            inds_slice = np.asarray(inds_slice[0])

            hp_indices_slice = indices[inds_slice]
            weightings_slice = weightings[inds_slice]

            numb_dens_map_slice = np.zeros(npix)

            positive_weightings_ids = np.where((weightings_slice>0))
            negative_weightings_ids = np.where((weightings_slice<0))
            #np.add.at(numb_dens_map_slice, hp_indices_slice, hp_indices_slice + 1)

            np.add.at(numb_dens_map_slice, hp_indices_slice[positive_weightings_ids], np.around(hp_indices_slice[positive_weightings_ids]*weightings_slice[positive_weightings_ids]) + 1)
            np.add.at(numb_dens_map_slice, hp_indices_slice[negative_weightings_ids], np.around(hp_indices_slice[negative_weightings_ids]*weightings_slice[negative_weightings_ids]) - 1)

            numb_dens_map_slice = numb_dens_map_slice / np.arange(1, npix + 1)

            ngal_at_slice = sum(weightings_slice)
            ngals_inbin.append(ngal_at_slice)

            k_map_slice = np.zeros(npix)
            y1_map_slice = np.zeros(npix)
            y2_map_slice = np.zeros(npix)

            np.add.at(k_map_slice, hp_indices_slice, conv[inds_slice])
            np.add.at(y1_map_slice, hp_indices_slice, shear1[inds_slice])
            np.add.at(y2_map_slice, hp_indices_slice, shear2[inds_slice])

            #np.add.at(k_map_slice, hp_indices_slice[positive_weightings_ids], conv[inds_slice][positive_weightings_ids])
            #np.add.at(y1_map_slice, hp_indices_slice[positive_weightings_ids], shear1[inds_slice][positive_weightings_ids])
            #np.add.at(y2_map_slice, hp_indices_slice[positive_weightings_ids], shear2[inds_slice][positive_weightings_ids])

            #np.add.at(k_map_slice, hp_indices_slice, conv[inds_slice]*weightings_slice)
            #np.add.at(y1_map_slice, hp_indices_slice, shear1[inds_slice]*weightings_slice)
            #np.add.at(y2_map_slice, hp_indices_slice, shear2[inds_slice]*weightings_slice)

            k_maps_inbin.append(k_map_slice)
            y1_maps_inbin.append(y1_map_slice)
            y2_maps_inbin.append(y2_map_slice)
            numb_dens_maps_inbin.append(numb_dens_map_slice)

        k_map = sum(k_maps_inbin)
        y1_map = sum(y1_maps_inbin)
        y2_map = sum(y2_maps_inbin)

        dens_map = sum(numb_dens_maps_inbin)

        #print([n[observed_pixels_ids] for n in numb_dens_maps_inbin])
        #print(dens_map[observed_pixels_ids])

        # Prepare dictionary entries for pCl estimates
        cut_maps_dic["BIN_{}".format(i + 1)] = []

        searchval = 0
        ii = np.where(dens_map == searchval)[0]
        if 0 in dens_map:
            if len(ii) != unobserved_pixels:

                with open(recov_cat_cls_dir + 'README.txt', 'a+') as f:
                    f.write('\nIter {}: WARNING! BIN {} HAS {} HEALPIX PIXELS NOT FILLED: INDICES {}'.format(
                        iter_no, (i + 1), len(ii), ii))
            else:
                pass

        n_gals = sum(ngals_inbin)
        #print(n_gals)
        #n_gals = len(inds)
        #print(n_gals)
        count_map = dens_map
        # Calculate average occupancy of each pixel
        delta_bar = n_gals / (npix - len(ii))
        # The clustering density field is then just the count in each pixel normalised to the average occupancy
        dens_map = (dens_map - delta_bar) / delta_bar

        # Calculate the theoretical noise cls based on observed number density
        pi = np.pi
        sky_coverage = hp.nside2pixarea(nside, degrees=True) * (npix - len(ii)) * ((pi / 180) ** 2)
        poisson_cls_theory = np.zeros(len(np.arange(lmin_out, lmax_out + 1)))
        poisson_cls_theory = poisson_cls_theory + (w_survey * sky_coverage / n_gals)

        np.savetxt(poisson_cls_theory_dir + 'bin_{}_{}.txt'.format(i + 1, i + 1),
                   np.transpose(poisson_cls_theory))

        np.savetxt(poisson_cls_theory_dir + 'ell.txt', np.transpose(ell_arr))

        # Need to calculate the Shear shape noise - dependent on number density in bin. See Upham 2021 Eq.8
        she_nl = ((((sigma_e / np.sqrt(2)) ** 2) / (n_gals / sky_coverage)) * w_survey) + np.zeros(
            len(np.arange(lmin_out, lmax_out + 1)))

        np.savetxt(shape_noise_cls_dir + 'bin_{}_{}.txt'.format(i + 1, i + 1),
                   np.transpose(she_nl))

        np.savetxt(shape_noise_cls_dir + 'ell.txt', np.transpose(ell_arr))

        # Now we need to normalise each of the weak lensing observables maps - take the average value in
        # each pixel, i.e. divide by the count in each pixel

        # Temporary mask so we are not dividing by zero
        count_map[ii] = 1

        # Now normalise each weak lensing map
        k_map = k_map / count_map
        y1_map = y1_map / count_map
        y2_map = y2_map / count_map

        # Change back temporary mask into Healpix unseen quantities so they can be plotted as null values
        k_map[ii] = hp.UNSEEN
        y1_map[ii] = hp.UNSEEN
        y2_map[ii] = hp.UNSEEN
        dens_map[ii] = hp.UNSEEN

        mask = np.zeros(npix)
        mask = mask + 1
        mask[ii] = 0

        # Now append each cut-sky map into an array in dictionary entry - format is e.g.
        cut_maps_dic["BIN_{}".format(i + 1)].append(
            nmt.NmtField(mask=mask, maps=[k_map], spin=0, lmax_sht=lmax_out))
        cut_maps_dic["BIN_{}".format(i + 1)].append(
            nmt.NmtField(mask=mask, maps=[y1_map, y2_map], spin=2, lmax_sht=lmax_out))
        cut_maps_dic["BIN_{}".format(i + 1)].append(
            nmt.NmtField(mask=mask, maps=[dens_map], spin=0, lmax_sht=lmax_out))

        # Plot each of the weak lensing maps
        plt.figure(figsize=(15, 3.5))
        hp.mollview(k_map,
                    sub=(131),
                    title=r'$\kappa$' + " Convergence ['TEMPERATURE']")
        hp.mollview(y1_map,
                    sub=(132),
                    title="Redshift {} < z < {}\n\n".format(z_boundaries_low[i], round(z_boundaries_high[i], 2))
                          + r'$\gamma$' + "1 Shear ['Q_POLARISATION']",
                    cmap='ocean')
        hp.mollview(y2_map,
                    sub=(133),
                    title=r'$\gamma$' + "2 Shear ['U_POLARISATION']",
                    cmap='Purples')
        plt.tight_layout()

        plt.savefig(shear_dir + 'map_zbin_%s.png' % ((i + 1)))

        # Plot density field
        plt.figure(figsize=(5, 3.5))
        hp.mollview(dens_map,
                    title="Galaxy Clustering (Matter Distribution) %.2f < z < %.2f" % (
                        z_boundaries_low[i], round(z_boundaries_high[i], 2))
                    )
        plt.savefig(gal_dir + 'map_zbin_%s.png' % (i + 1))
        plt.close()


    return cut_maps_dic


def measure_00_pcls(lmin_out, lmax_out, cut_maps_dic, spectra_type, bin_i, bin_j, measured_pcl_save_dir):

    accepted_spectra_types = {'TT', 'gal_gal'}
    if spectra_type not in accepted_spectra_types:
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    nmt_fields = {
        'TT': [cut_maps_dic["BIN_{}".format(bin_i)][0],
               cut_maps_dic["BIN_{}".format(bin_j)][0],
               ],

        'gal_gal': [cut_maps_dic["BIN_{}".format(bin_i)][-1],
                    cut_maps_dic["BIN_{}".format(bin_j)][-1],
                    ],
    }

    pcl_coupled = nmt.compute_coupled_cell(nmt_fields[spectra_type][0], nmt_fields[spectra_type][1])

    measured_pcl = pcl_coupled[0][lmin_out:lmax_out + 1]

    np.savetxt(measured_pcl_save_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j),
               np.transpose(measured_pcl))

    ells = np.arange(lmin_out, lmax_out+1, 1)

    #np.savetxt(measured_pcl_save_dir + 'ell.txt', np.transpose(ells))


def measure_02_pcls(lmin_out, lmax_out, cut_maps_dic, spectra_type, bin_i, bin_j, measured_pcl_save_dir):

    accepted_spectra_types = {'TE', 'TB', 'gal_E', 'gal_B'}
    if spectra_type not in accepted_spectra_types:
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    nmt_fields = {
        'TE': [cut_maps_dic["BIN_{}".format(bin_i)][0],
               cut_maps_dic["BIN_{}".format(bin_j)][1],
               ],

        'TB': [cut_maps_dic["BIN_{}".format(bin_i)][0],
               cut_maps_dic["BIN_{}".format(bin_j)][1],
               ],

        'gal_E': [cut_maps_dic["BIN_{}".format(bin_i)][-1],
                  cut_maps_dic["BIN_{}".format(bin_j)][1],
                  ],

        'gal_B': [cut_maps_dic["BIN_{}".format(bin_i)][-1],
                  cut_maps_dic["BIN_{}".format(bin_j)][1],
                  ]
    }

    pcl_coupled = nmt.compute_coupled_cell(nmt_fields[spectra_type][0], nmt_fields[spectra_type][1])
    measured_pcls = defaultdict(list)

    if spectra_type == 'TE' or 'TB':
        measured_pcls['TE'] = (pcl_coupled[0][lmin_out:lmax_out + 1])
        measured_pcls['TB'] = (pcl_coupled[1][lmin_out:lmax_out + 1])

    if spectra_type == 'gal_E' or 'gal_B':
        measured_pcls['gal_E'] = (pcl_coupled[0][lmin_out:lmax_out + 1])
        measured_pcls['gal_B'] = (pcl_coupled[1][lmin_out:lmax_out + 1])

    np.savetxt(measured_pcl_save_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j),
               np.transpose(measured_pcls[spectra_type]))

    ells = np.arange(lmin_out, lmax_out+1, 1)

    #np.savetxt(measured_pcl_save_dir + 'ell.txt', np.transpose(ells))


def measure_22_pcls(lmin_out, lmax_out, cut_maps_dic, spectra_type, bin_i, bin_j, measured_pcl_save_dir):

    accepted_spectra_types = {'EE', 'EB', 'BE', 'BB'}
    if spectra_type not in accepted_spectra_types:
        # print(spectra_type)
        print('Warning! Field Type Not Recognised - Exiting...')
        sys.exit()

    pcl_coupled = nmt.compute_coupled_cell(cut_maps_dic["BIN_{}".format(bin_i)][1],
                                           cut_maps_dic["BIN_{}".format(bin_j)][1])

    measured_pcl_components = {
        'EE': (pcl_coupled[0][lmin_out:lmax_out + 1]),
        'EB': (pcl_coupled[1][lmin_out:lmax_out + 1]),
        'BE': (pcl_coupled[2][lmin_out:lmax_out + 1]),
        'BB': (pcl_coupled[3][lmin_out:lmax_out + 1])
    }

    np.savetxt(measured_pcl_save_dir + 'bin_{}_{}.txt'.format(bin_i, bin_j),
               np.transpose(measured_pcl_components[spectra_type]))

    ells = np.arange(lmin_out, lmax_out+1, 1)

    #np.savetxt(measured_pcl_save_dir + 'ell.txt', np.transpose(ells))


def execute_pcl_measurement():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']
    config_dict = measure_pcls_config(pipeline_variables_path=pipeline_variables_path)

    realisation = os.environ['ITER_NO']

    cut_sky_map_dicts = maps_from_cats(config_dict=config_dict, iter_no=realisation)

    nbins = config_dict['nbins']
    raw_pcl_lmin_out = config_dict['raw_pcl_lmin_out']
    raw_pcl_lmax_out = config_dict['raw_pcl_lmax_out']
    save_dir = config_dict['save_dir']

    recov_cat_cls_dir = save_dir + 'raw_3x2pt_cls/'

    gal_dir = recov_cat_cls_dir + 'galaxy_cl/iter_{}/'.format(realisation)
    gal_shear_dir = recov_cat_cls_dir + 'galaxy_shear_cl/iter_{}/'.format(realisation)

    shear_dir = recov_cat_cls_dir + 'shear_cl/'
    k_dir = shear_dir + 'Cl_TT/iter_{}/'.format(realisation)
    y1_dir = shear_dir + 'Cl_EE/iter_{}/'.format(realisation)
    y1y2_dir = shear_dir + 'Cl_EB/iter_{}/'.format(realisation)
    y2y1_dir = shear_dir + 'Cl_BE/iter_{}/'.format(realisation)
    y2_dir = shear_dir + 'Cl_BB/iter_{}/'.format(realisation)

    # Let's measure some Pseudo-Cls
    for bin_i in range(nbins):
        for bin_j in range(nbins):

            # Galaxy-Galaxy Lensing
            measure_02_pcls(
                lmin_out=raw_pcl_lmin_out,
                lmax_out=raw_pcl_lmax_out,
                cut_maps_dic=cut_sky_map_dicts,
                spectra_type='gal_E',
                bin_i=bin_i + 1,
                bin_j=bin_j + 1,
                measured_pcl_save_dir=gal_shear_dir
            )

            if bin_i >= bin_j:
                # Weak Lensing Convergence
                measure_00_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='TT',
                    bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=k_dir
                )

                # Galaxy Clustering
                measure_00_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='gal_gal',
                    bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=gal_dir)

                # Weak Lensing E-Mode
                measure_22_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='EE',
                    bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=y1_dir
                )

                # Weak Lensing EB Cross
                measure_22_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='EB', bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=y1y2_dir
                )

                # Weak Lensing BE Cross
                measure_22_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='BE',
                    bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=y2y1_dir
                )

                # Weak Lensing B-Mode
                measure_22_pcls(
                    lmin_out=raw_pcl_lmin_out,
                    lmax_out=raw_pcl_lmax_out,
                    cut_maps_dic=cut_sky_map_dicts,
                    spectra_type='BB',
                    bin_i=bin_i + 1,
                    bin_j=bin_j + 1,
                    measured_pcl_save_dir=y2_dir
                )

    # Fill in the remaining theory Cls - off diagonal spectra will be zero:

    gal_shear_noise_cls_dir = save_dir + 'raw_noise_cls/galaxy_shear_cl/iter_%s/' % realisation
    if not os.path.exists(gal_shear_noise_cls_dir):
        os.makedirs(gal_shear_noise_cls_dir)

    ell_arr = np.arange(raw_pcl_lmin_out, raw_pcl_lmax_out + 1, 1)

    np.savetxt(save_dir + 'raw_noise_cls/galaxy_cl/ell.txt',
               np.transpose(ell_arr))

    np.savetxt(save_dir + 'raw_noise_cls/shear_cl/ell.txt',
               np.transpose(ell_arr))

    np.savetxt(save_dir + 'raw_noise_cls/galaxy_shear_cl/ell.txt',
               np.transpose(ell_arr))

    for i in range(nbins):
        for j in range(nbins):
            null_noise_cls = np.zeros(len(ell_arr))
            np.savetxt(
                gal_shear_noise_cls_dir + 'bin_%s_%s.txt' % (i + 1, j + 1),
                np.transpose(null_noise_cls)
            )

            np.savetxt(gal_shear_noise_cls_dir + 'ell.txt', np.transpose(ell_arr))

            if i > j:
                poisson_cls_theory_dir = save_dir + 'raw_noise_cls/galaxy_cl/iter_{}/'.format(realisation)
                np.savetxt(
                    poisson_cls_theory_dir + 'bin_%s_%s.txt' % (i + 1, j + 1),
                    np.transpose(null_noise_cls)
                )

                shape_noise_cls_dir = save_dir + 'raw_noise_cls/shear_cl/iter_{}/'.format(realisation)
                np.savetxt(
                    shape_noise_cls_dir + 'bin_%s_%s.txt' % (i + 1, j + 1),
                    np.transpose(null_noise_cls)
                )


if __name__ == '__main__':
    execute_pcl_measurement()
