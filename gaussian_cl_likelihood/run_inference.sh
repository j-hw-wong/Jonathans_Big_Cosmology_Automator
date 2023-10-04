#!/bin/bash

start=$SECONDS

PIPELINE_VARIABLES_PATH="/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim/set_variables.ini"
export PIPELINE_VARIABLES_PATH

source <(grep = $PIPELINE_VARIABLES_PATH)
export SAVE_DIR
export INFERENCE_PIPELINE_DIR
export INFERENCE_OUTPUT_DIR
export COSMOSIS_ROOT_DIR
export ELL_MIN
export ELL_MAX
export NZ_TABLE_FILENAME

#COSMOSIS_N_ELL="$("$ELL_MAX-$ELL_MIN" |bc)"
#COSMOSIS_N_ELL=$ELL_MAX - $ELL_MIN + 1| bc
COSMOSIS_N_ELL=$(echo "($ELL_MAX - $ELL_MIN + 1)/1" | bc)
#echo $COSMOSIS_N_ELL
export COSMOSIS_N_ELL

#python fig3.py

#cd ./bash/
#bash run_cosmosis_chains.sh