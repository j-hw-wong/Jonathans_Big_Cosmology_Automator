#!/bin/bash

COSMOSIS_DIR=$COSMOSIS_ROOT_DIR

cd $COSMOSIS_DIR
source config/setup-cosmosis

make

echo $PIPELINE_DIR
COSMOSIS_PIPELINE_PATH=${PIPELINE_DIR}cosmosis_config.ini
COSMOSIS_VARIABLES_PATH=${PIPELINE_DIR}cosmosis_params.ini
export COSMOSIS_PIPELINE_PATH
export COSMOSIS_VARIABLES_PATH

COSMOSIS_OUT_DIR=${SAVE_DIR}cosmosis/
mkdir -p ${COSMOSIS_OUT_DIR}
export COSMOSIS_OUT_DIR

COSMOSIS_ELL_MAX=$(echo $ELL_MAX+1.0|bc)
COSMOSIS_N_ELL=$(echo "($ELL_MAX - $ELL_MIN + 1)/1" | bc)
export COSMOSIS_ELL_MAX
export COSMOSIS_N_ELL

NZ_TABLE_PATH="$SAVE_DIR$NZ_TABLE_FILENAME"
export NZ_TABLE_PATH

cosmosis $COSMOSIS_PIPELINE_PATH
