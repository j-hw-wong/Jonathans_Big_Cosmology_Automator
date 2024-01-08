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
		spec_dir = measured_cls_dir + 'galaxy_bp/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'iter_{}/bin_{}_{}.txt'.format(it + 1, spec_2_zbin, spec_1_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'E' and spec_2_field == 'E':
		spec_dir = measured_cls_dir + 'shear_bp/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'Cl_EE/iter_{}/bin_{}_{}.txt'.format(it + 1, spec_2_zbin, spec_1_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'N' and spec_2_field == 'E':
		spec_dir = measured_cls_dir + 'galaxy_shear_bp/'
		cov_dat = []
		for it in range(no_iter):
			dat = np.loadtxt(
				spec_dir + 'iter_{}/bin_{}_{}.txt'.format(it + 1, spec_1_zbin, spec_2_zbin))
			cov_dat.append(dat)

	elif spec_1_field == 'E' and spec_2_field == 'N':
		spec_dir = measured_cls_dir + 'galaxy_shear_bp/'
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
	measured_bps_dir = save_dir + 'measured_3x2pt_bps/'
	no_iter = int(config['measurement_setup']['REALISATIONS'])
	n_zbin = int(config['create_nz']['N_ZBIN'])
	n_bp = int(config['measurement_setup']['N_BANDPOWERS'])
	obs_type = str(config['measurement_setup']['OBS_TYPE'])
	field_type = str(config['measurement_setup']['FIELD'])

	if obs_type == '3X2PT':
		n_field = 2 * n_zbin
		fields = [f'{f}{z}' for z in range(1, n_zbin + 1) for f in ['N', 'E']]

	else:
		assert obs_type == '1X2PT'
		n_field = n_zbin
		if field_type == 'E':
			fields = [f'E{z}' for z in range(1, n_zbin + 1)]
		else:
			assert field_type == 'N'
			fields = [f'N{z}' for z in range(1, n_zbin + 1)]

	spectra = [fields[row] + fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]

	spec_1 = [fields[row] for diag in range(n_field) for row in range(n_field - diag)]
	spec_2 = [fields[row + diag] for diag in range(n_field) for row in range(n_field - diag)]

	full_cov_dat = []
	cov_data_av = np.array([])

	for i in range(len(spectra)):
		spectrum_a = spec_1[i]
		spectrum_b = spec_2[i]

		spectrum_i_dat = open_spectrum(
			spectrum_a,
			spectrum_b,
			measured_cls_dir=measured_bps_dir,
			no_iter=no_iter)
		full_cov_dat.append(spectrum_i_dat)
		spectrum_i_dat_av = np.mean(np.array(spectrum_i_dat), axis=0)
		cov_data_av = np.concatenate((cov_data_av, spectrum_i_dat_av), axis=0)

	cov_total = []
	for n in range(no_iter):
		this_cov_dat = np.array([])
		for it in range(len(full_cov_dat)):
			this_cov_dat = np.concatenate((this_cov_dat, full_cov_dat[it][n]), axis=0)

		cov_iter = np.zeros([len(this_cov_dat), len(this_cov_dat)])

		for x in range(len(this_cov_dat)):
			for y in range(len(this_cov_dat)):
				cov_iter[x][y] = (this_cov_dat[x] - cov_data_av[x])*(this_cov_dat[y] - cov_data_av[y])
		cov_total.append(cov_iter)

	cov = np.mean(np.array(cov_total), axis=0)
	# cov = np.abs(cov)
	# cov = np.diag(np.diag(cov))
	save_sim_cov_dir = save_dir + 'cov_fromsim/'

	if not os.path.exists(save_sim_cov_dir):
		os.makedirs(save_sim_cov_dir)

	np.savez(
		save_sim_cov_dir + 'cov_{}bp.npz'.format(n_bp),
		cov=cov)


if __name__ == '__main__':
	main()
