#!/bin/bash

cd /raid/scratch/wongj/mywork/3x2pt/Jonathans_Big_Cosmology_Automator/catalogue_sim/

bash run_cat_sim.sh

cd /raid/scratch/wongj/mywork/3x2pt/Jonathans_Big_Cosmology_Automator/pcl_measurement/

bash run_3x2pt_tomo_measurement.sh

cd /raid/scratch/wongj/mywork/3x2pt/Jonathans_Big_Cosmology_Automator/inference_analysis/

bash run_inference.sh