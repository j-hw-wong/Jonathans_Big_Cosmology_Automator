import os
import sys
import configparser
import numpy as np


def mysplit(s):
	head = s.rstrip('0123456789')
	tail = s[len(head):]
	return head, tail


def open_spectrum(id_a, id_b, measured_cls_dir, no_iter):
	# id_a is the first signal i.e. N1
	# id_b is the second signal i.e. E2
	# so the power spectrum of interest is N1E2

	spec_1_field = mysplit(id_a)[0]
	spec_1_zbin = mysplit(id_a)[1]

	spec_2_field = mysplit(id_b)[0]
	spec_2_zbin = mysplit(id_b)[1]

	if spec_1_field == 'N' and spec_2_field == 'N':
		spec_dir = measured_cls_dir + 'galaxy_cl/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'iter_{}/bin_{}_{}.txt'.format(it + 1, spec_2_zbin, spec_1_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'E' and spec_2_field == 'E':
		spec_dir = measured_cls_dir + 'shear_cl/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'Cl_EE/iter_{}/bin_{}_{}.txt'.format(it + 1, spec_2_zbin, spec_1_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'N' and spec_2_field == 'E':
		spec_dir = measured_cls_dir + 'galaxy_shear_cl/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'iter_{}/bin_{}_{}.txt'.format(it + 1, spec_1_zbin, spec_2_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'E' and spec_2_field == 'N':
		spec_dir = measured_cls_dir + 'galaxy_shear_cl/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'iter_{}/bin_{}_{}.txt'.format(it + 1, spec_2_zbin, spec_1_zbin))

			cov_dat.append(dat)
	else:
		sys.exit()

	return cov_dat


def main():
	pipeline_variables_path = os.environ['PIPELINE_VARIABLES_PATH']

	config = configparser.ConfigParser()
	config.read(pipeline_variables_path)

	save_dir = str(config['measurement_setup']['MEASUREMENT_SAVE_DIR'])
	measured_cls_dir = save_dir + 'measured_3x2pt_cls/'
	no_iter = int(config['measurement_setup']['REALISATIONS'])
	n_zbin = int(config['create_nz']['N_ZBIN'])

	obs_type = str(config['measurement_setup']['OBS_TYPE'])
	field_type = str(config['measurement_setup']['FIELD'])

	if obs_type == '3X2PT':
		n_field = 2 * n_zbin
		fields = [f'{f}{z}' for z in range(1, n_zbin + 1) for f in ['N', 'E']]

	else:
		assert obs_type == '1X2PT'
		n_field = n_zbin
		if field_type == 'E':
			fields = [f'{f}{z}' for z in range(1, n_zbin + 1) for f in ['E']]

		else:
			assert field_type == 'N'
			fields = [f'{f}{z}' for z in range(1, n_zbin + 1) for f in ['N']]

	spectra = [fields[row] + fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]
	spec_1 = [fields[row] for diag in range(n_field) for row in range(n_field - diag)]
	spec_2 = [fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]

	for i in range(len(spectra)):
		for j in range(len(spectra)):
			if i >= j:
				spectrum_i_a = spec_1[i]
				spectrum_i_b = spec_2[i]

				spectrum_j_a = spec_1[j]
				spectrum_j_b = spec_2[j]

				spectrum_i_dat = open_spectrum(
					spectrum_i_a,
					spectrum_i_b,
					measured_cls_dir=measured_cls_dir,
					no_iter=no_iter)
				spectrum_i_dat_av = np.mean(np.array(spectrum_i_dat), axis=0)

				spectrum_j_dat = open_spectrum(
					spectrum_j_a,
					spectrum_j_b,
					measured_cls_dir=measured_cls_dir,
					no_iter=no_iter)
				spectrum_j_dat_av = np.mean(np.array(spectrum_j_dat), axis=0)

				total_cov_dat = []
				for c in range(len(spectrum_i_dat)):
					cov_iter = np.matmul(
						np.transpose([np.asarray(spectrum_i_dat[c] - spectrum_i_dat_av)]),
						[np.asarray(spectrum_j_dat[c] - spectrum_j_dat_av)])
					total_cov_dat.append(cov_iter)

				cov = np.mean(np.array(total_cov_dat), axis=0)
				# cov = np.abs(cov)
				# cov = np.diag(np.diag(cov))
				save_sim_cov_dir = save_dir + 'cov_fromsim/'

				if not os.path.exists(save_sim_cov_dir):
					os.makedirs(save_sim_cov_dir)

				np.savez(
					save_sim_cov_dir + 'cov_spec1_{}_spec2_{}.npz'.format(i, j),
					cov_block=cov,
					spec1_idx=i,
					spec2_idx=j)


if __name__ == '__main__':
	main()
