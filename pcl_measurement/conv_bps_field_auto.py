import os
import sys
import configparser
import numpy as np

'''
Collect all 3x2pt observed spectra (pseudo-Cl bandpowers) into one file for input into inference code
'''


def mysplit(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return head, tail


def conv_3x2pt_bps(n_zbin, n_bp, save_dir, recov_cat_bps_path, obs_type='obs'):
    n_field = 2 * n_zbin
    n_spec = n_field * (n_field + 1) // 2

    # Form list of power spectra
    fields = [f'{f}{z}' for z in range(1, n_zbin + 1) for f in ['N', 'E']]
    print(fields)
    assert len(fields) == n_field

    spectra = [fields[row] + fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]
    spec_1 = [fields[row] for diag in range(n_field) for row in range(n_field - diag)]
    spec_2 = [fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]

    field_1 = [mysplit(spec_1_id)[0] for spec_1_id in spec_1]
    zbin_1 = [mysplit(spec_1_id)[1] for spec_1_id in spec_1]
    field_2 = [mysplit(spec_2_id)[0] for spec_2_id in spec_2]
    zbin_2 = [mysplit(spec_2_id)[1] for spec_2_id in spec_2]

    obs_bp = np.full((n_spec, n_bp), np.nan)

    for spec_idx in range(len(spectra)):
        f1 = field_1[spec_idx]
        f2 = field_2[spec_idx]
        z1 = zbin_1[spec_idx]
        z2 = zbin_2[spec_idx]

        if obs_type == 'obs':

            if f1 == 'N' and f2 == 'N':
                measured_bp_file = recov_cat_bps_path + 'galaxy_bp/bin_{}_{}.txt'.format(max(z1, z2), min(z1, z2))
                #measured_bp_file = save_dir + 'theory_cls/galaxy_cl/PCl_Bandpowers_gal_gal_bin_{}_{}.txt'.\
                #    format(max(z1, z2), min(z1, z2))

            elif f1 == 'E' and f2 == 'E':
                measured_bp_file = recov_cat_bps_path + 'shear_bp/Cl_EE/bin_{}_{}.txt'.format(max(z1, z2), min(z1, z2))
                #measured_bp_file = save_dir + 'theory_cls/shear_cl/PCl_Bandpowers_EE_bin_{}_{}.txt'. \
                #    format(max(z1, z2), min(z1, z2))

            elif f1 == 'N' and f2 == 'E':
                measured_bp_file = recov_cat_bps_path + 'galaxy_shear_bp/bin_{}_{}.txt'.format(z1, z2)
                #measured_bp_file = save_dir + 'theory_cls/galaxy_shear_cl/PCl_Bandpowers_gal_E_bin_{}_{}.txt'. \
                #    format(z1, z2)

            elif f1 == 'E' and f2 == 'N':
                # switch around, i.e. open f2z2f1z1
                measured_bp_file = recov_cat_bps_path + 'galaxy_shear_bp/bin_{}_{}.txt'.format(z2, z1)
                #measured_bp_file = save_dir + 'theory_cls/galaxy_shear_cl/PCl_Bandpowers_gal_E_bin_{}_{}.txt'. \
                #    format(z2, z1)

            else:
                print('Unexpected Spectra Found - Please check inference pipeline')
                sys.exit()

            measured_bp = np.loadtxt(measured_bp_file)
            obs_bp[spec_idx, :] = measured_bp

        else:

            assert obs_type == 'fid'

            if f1 == 'N' and f2 == 'N':
                measured_bp_file = save_dir + 'theory_cls/galaxy_cl/PCl_Bandpowers_gal_gal_bin_{}_{}.txt'. \
                    format(max(z1, z2), min(z1, z2))

            elif f1 == 'E' and f2 == 'E':
                measured_bp_file = save_dir + 'theory_cls/shear_cl/PCl_Bandpowers_EE_bin_{}_{}.txt'. \
                    format(max(z1, z2), min(z1, z2))

            elif f1 == 'N' and f2 == 'E':
                measured_bp_file = save_dir + 'theory_cls/galaxy_shear_cl/PCl_Bandpowers_gal_E_bin_{}_{}.txt'. \
                    format(z1, z2)

            elif f1 == 'E' and f2 == 'N':
                # switch around, i.e. open f2z2f1z1
                measured_bp_file = save_dir + 'theory_cls/galaxy_shear_cl/PCl_Bandpowers_gal_E_bin_{}_{}.txt'. \
                    format(z2, z1)

            else:
                print('Unexpected Spectra Found - Please check inference pipeline')
                sys.exit()

            measured_bp = np.loadtxt(measured_bp_file)
            obs_bp[spec_idx, :] = measured_bp

    obs_bp_path = os.path.join(recov_cat_bps_path, f'obs_{n_bp}bp.npz')
    obs_bp_header = (f'Observed bandpowers for 3x2pt simulation')
    np.savez_compressed(obs_bp_path, obs_bp=obs_bp, header=obs_bp_header)


def conv_1x2pt_bps(n_zbin, n_bp, save_dir, recov_cat_bps_path, obs_type='obs', field='E'):
    # Assuming a shear-only analysis

    n_field = n_zbin
    n_spec = n_field * (n_field + 1) // 2

    # Form list of power spectra
    if field == 'E':
        fields = [f'E{z}' for z in range(1, n_zbin + 1)]
    else:
        assert field == 'N'
        fields = [f'N{z}' for z in range(1, n_zbin + 1)]

    assert len(fields) == n_field
    spectra = [fields[diag] + fields[diag] for diag in range(n_field)]
    print(spectra)
    spec_1 = [fields[diag] for diag in range(n_field)]
    spec_2 = [fields[diag] for diag in range(n_field)]
    #print(spec_1)
    #print(spec_2)
    field_1 = [mysplit(spec_1_id)[0] for spec_1_id in spec_1]
    zbin_1 = [mysplit(spec_1_id)[1] for spec_1_id in spec_1]
    field_2 = [mysplit(spec_2_id)[0] for spec_2_id in spec_2]
    zbin_2 = [mysplit(spec_2_id)[1] for spec_2_id in spec_2]

    obs_bp = np.full((n_spec, n_bp), np.nan)

    for spec_idx in range(len(spectra)):
        f1 = field_1[spec_idx]
        f2 = field_2[spec_idx]
        z1 = zbin_1[spec_idx]
        z2 = zbin_2[spec_idx]
        #assert f1 == 'E'
        #assert f2 == 'E'

        if obs_type == 'obs':

            if field == 'E':
                measured_bp_file = recov_cat_bps_path + 'shear_bp/Cl_EE/bin_{}_{}.txt'.format(max(z1, z2), min(z1, z2))
                measured_bp = np.loadtxt(measured_bp_file)
                obs_bp[spec_idx, :] = measured_bp

            else:
                assert field == 'N'
                measured_bp_file = recov_cat_bps_path + 'galaxy_bp/bin_{}_{}.txt'.format(max(z1, z2), min(z1, z2))
                measured_bp = np.loadtxt(measured_bp_file)
                obs_bp[spec_idx, :] = measured_bp
        else:
            assert obs_type == 'fid'
            if field == 'E':
                measured_bp_file = save_dir + 'theory_cls/shear_cl/PCl_Bandpowers_EE_bin_{}_{}.txt'. \
                    format(max(z1, z2), min(z1, z2))
                measured_bp = np.loadtxt(measured_bp_file)
                obs_bp[spec_idx, :] = measured_bp

            else:
                assert field == 'N'
                measured_bp_file = save_dir + 'theory_cls/galaxy_cl/PCl_Bandpowers_gal_gal_bin_{}_{}.txt'. \
                    format(max(z1, z2), min(z1, z2))
                measured_bp = np.loadtxt(measured_bp_file)
                obs_bp[spec_idx, :] = measured_bp

    obs_bp_path = os.path.join(recov_cat_bps_path, f'obs_{n_bp}bp.npz')
    obs_bp_header = (f'Observed bandpowers for 1x2pt simulation')
    np.savez_compressed(obs_bp_path, obs_bp=obs_bp, header=obs_bp_header)


def main():
    pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']

    config = configparser.ConfigParser()
    config.read(pipeline_variables_path)

    save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
    recov_cat_bps_path = save_dir + 'measured_3x2pt_bps/'

    n_zbin = int(float(config['create_nz']['N_ZBIN']))
    n_bp = int(float(config['measurement_setup']['N_BANDPOWERS']))

    obs_type = str(config['measurement_setup']['OBS_TYPE'])
    obs_field = str(config['measurement_setup']['FIELD'])

    if obs_type == '3X2PT':
        conv_3x2pt_bps(
            n_zbin=n_zbin,
            n_bp=n_bp,
            save_dir=save_dir,
            recov_cat_bps_path=recov_cat_bps_path,
            obs_type='obs')

    elif obs_type == '1X2PT':
        conv_1x2pt_bps(
            n_zbin=n_zbin,
            n_bp=n_bp,
            save_dir=save_dir,
            recov_cat_bps_path=recov_cat_bps_path,
            obs_type='obs',
            field=obs_field)


if __name__ == '__main__':
    main()
