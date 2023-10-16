import os
import configparser
import sys

angular_binning_path = os.environ['ANGULAR_BINNING_PATH']
gaussian_cl_likelihood_path = os.environ['GAUSSIAN_CL_LIKELIHOOD_PATH']
pipeline_dir = os.environ['PIPELINE_DIR']

sys.path.insert(1, pipeline_dir)
sys.path.insert(1, angular_binning_path)
sys.path.insert(1, gaussian_cl_likelihood_path)

from gaussian_cl_likelihood.python import cosmosis_utils

config = configparser.ConfigParser()
config.read(pipeline_variables_path)

save_dir = str(config['inference_analysis_params']['MEASUREMENT_SAVE_DIR'])

params = {
    'cosmological_parameters--w': {
        'min': -1.05,
        'max': -0.95,
        'steps': 25
    },
    'cosmological_parameters--wa': {
        'min': -0.1,
        'max': 0.1,
        'steps': 25
    }
}

n_chains = 11
output_dir = save_dir + 'inference_analysis/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

cosmosis_utils.generate_chain_input(params, n_chains, output_dir)
