#!/bin/bash

# Bash script to submit N_CHAINS cosmosis chains using nohup,
# with output from each chain to a separate log file.
# Pipeline must be set up to use the CHAIN_NO environment variable.

#PIPELINE_VARIABLES_PATH="/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim/set_variables.ini"
#export PIPELINE_VARIABLES_PATH

#source <(grep = $PIPELINE_VARIABLES_PATH)
#export SAVE_DIR
#export INFERENCE_PIPELINE_DIR
#export INFERENCE_OUTPUT_DIR
#export COSMOSIS_ROOT_DIR
#export ELL_MIN
#export ELL_MAX


# Parameters to set before using
N_CHAINS=11
COSMOSIS_DIR=/raid/scratch/wongj/cosmosis
#echo $INFERENCE_PIPELINE_DIR
#COSMOSIS_DIR=$COSMOSIS_ROOT_DIR
PIPELINE_PATH=/raid/scratch/wongj/mywork/3x2pt/gaussian_cl_likelihood/ini/tomo_3x2_pipeline_0.ini
#PIPELINE_PATH=${INFERENCE_PIPELINE_DIR}ini/tomo_3x2_pipeline_obs.ini
#echo $PIPELINE_PATH
#echo $INFERENCE_PIPELINE_DIR
#NZ_TABLE_PATH=${SAVE_DIR}${NZ_TABLE_FILENAME}
#echo $NZ_TABLE_PATH
#export NZ_TABLE_PATH

cd $COSMOSIS_DIR

source config/setup-cosmosis

# Loop
for ((CHAIN_NO = 0; CHAIN_NO < $N_CHAINS; CHAIN_NO++))
do
   CHAIN_LOG="chain$CHAIN_NO.out"
   export CHAIN_NO
   #echo $PIPELINE_PATH
   #cosmosis $PIPELINE_PATH
   nohup cosmosis $PIPELINE_PATH > $CHAIN_LOG 2>&1 &
   echo "Chain $CHAIN_NO running with PID $!"
done

# Tidy up
unset CHAIN_NO
echo "Done"

wait
