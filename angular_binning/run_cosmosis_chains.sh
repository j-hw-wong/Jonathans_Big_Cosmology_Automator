#!/bin/bash

start=$SECONDS

PIPELINE_VARIABLES_PATH="/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim/set_variables.ini"
GAUSSIAN_CL_LIKELIHOOD_PATH="/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood/"
export PIPELINE_VARIABLES_PATH

#source <(grep = $PIPELINE_VARIABLES_PATH)
#export SAVE_DIR
#export INFERENCE_PIPELINE_DIR
#export INFERENCE_OUTPUT_DIR
#export COSMOSIS_ROOT_DIR
#export ELL_MIN
#export ELL_MAX

SAVE_DIR="/raid/scratch/wongj/mywork/3x2pt/ang_bin_fig5/"
INFERENCE_PIPELINE_DIR="/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood/"
INFERENCE_OUTPUT_DIR="chains/"
COSMOSIS_ROOT_DIR="/raid/scratch/wongj/cosmosis/"
ELL_MIN=2.0
ELL_MAX=75.0

export SAVE_DIR
export NZ_TABLE_FILENAME
export INFERENCE_PIPELINE_DIR
export INFERENCE_OUTPUT_DIR
export COSMOSIS_ROOT_DIR
export ELL_MIN
export ELL_MAX

#echo ${INFERENCE_PIPELINE_DIR}cosmosis-standard-library/structure/projection/project_2d.py
NZ_TABLE_PATH=${SAVE_DIR}${NZ_TABLE_FILENAME}
export NZ_TABLE_PATH

COSMOSIS_N_ELL=$(echo "($ELL_MAX - $ELL_MIN + 1)/1" | bc)
export COSMOSIS_N_ELL

cd ${GAUSSIAN_CL_LIKELIHOOD_PATH}bash/

bash run_cosmosis_chains.sh