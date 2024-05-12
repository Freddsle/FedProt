#!/bin/bash

base_path=${1:-"/home/yuliya/repos/cosybio/FedProt/tmt_simulation_test_gradually"}
mode=${2:-"balanced"} 
number_of_files=${3:-"50"}

# Navigate to the directory where the files and scripts are located
cd ${base_path}/${mode}/meta

# Loop through each set of files from 1 to 50
for j in $(seq 1 $number_of_files)
do
    # Run the R scripts for meta-analysis on each triple of files
    Rscript /home/yuliya/repos/cosybio/FedProt/evaluation_utils/meta_code/run_MetaDE.R ./ ${j}_lab1 ${j}_lab2 ${j}_lab3
    Rscript /home/yuliya/repos/cosybio/FedProt/evaluation_utils/meta_code/run_MetaVolcanoR.R ./ ${j}_lab1 ${j}_lab2 ${j}_lab3
    Rscript /home/yuliya/repos/cosybio/FedProt/evaluation_utils/meta_code/run_RankProd.R ./ ${j}_lab1 ${j}_lab2 ${j}_lab3

    # Copy the MA_* files to a new location and rename them by adding the prefix corresponding to j
    for file in MA_*
    do
        cp "$file" "${base_path}${mode}/results/${j}_$file"
    done
done

echo "Meta-analysis and file copying completed."
