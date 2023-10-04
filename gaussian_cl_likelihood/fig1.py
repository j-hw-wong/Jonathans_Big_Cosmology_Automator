from python import cosmosis_utils
from python import simulation

'''
#Step a) 1
params = {
    'cosmological_parameters--w': {
        'min': -1.06,
        'max': -0.94,
        'steps': 99
    }
}
n_chains = 11 # for a 12-CPU machine, to leave one free
output_dir = './fig1/'

cosmosis_utils.generate_chain_input(params, n_chains, output_dir)
'''

'''
#Step a) 3
input_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/chains/'
chain_subdir_mask = 'chain{i}/'
filemask = '_{n}.tgz'
clean = False

cosmosis_utils.combine_chain_output(input_dir, chain_subdir_mask, filemask, clean)
'''

'''
#Step a) 4
input_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/chains/'
filemask = '_*.tgz'
params_filename = 'cosmological_parameters/values.txt'
nodelete = False # or True if you want to retain the .tgz files in case you want to extract anything else later
nbin_3x2 = 5

cosmosis_utils.extract_data(input_dir, filemask=filemask, params_filename=params_filename, nodelete=nodelete, nbin_3x2=nbin_3x2)
'''

'''
#Step a) 5
n_zbin = 5
gals_per_sq_arcmin = 30 / n_zbin
sigma_e = 0.3
lmin = 2
lmax = 2000
save_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/noise_cls'


simulation.noise_cls(gals_per_sq_arcmin, sigma_e, lmin, lmax, save_dir)
'''

'''
#Step b)
n_zbin = 5
fiducial_cls_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/chains/_0' # This can be one of the sub-directories in the grid produced in step (a).
pos_pos_in_dir = fiducial_cls_dir + '/galaxy_cl'
she_she_in_dir = fiducial_cls_dir + '/shear_cl'
pos_she_in_dir = fiducial_cls_dir + '/galaxy_shear_cl'
lmax = 2000
lmin_in = 2
pos_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/noise_cls/pos_nl.txt'
she_nl_path = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1/noise_cls/she_nl.txt'
lmin_out = 2
n_loop = 27000
batch_size = 1000
save_dir = '/raid/scratch/wongj/mywork/gaussian_cl_likelihood_jw/fig1'

simulation.sim_cls_fullsky(n_zbin, pos_pos_in_dir, she_she_in_dir, pos_she_in_dir, lmax, lmin_in, pos_nl_path, she_nl_path,
                    lmin_out, n_loop, batch_size, save_dir)
'''
