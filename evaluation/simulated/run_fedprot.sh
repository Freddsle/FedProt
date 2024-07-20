#!/bin/bash

data_path=${1:-"/home/yuliya/repos/cosybio/FedProt/tmt_simulation_test_gradually"}
base_path=${2:-"/home/yuliya/repos/cosybio/FedProt/tmt_simulation_test_gradually"}
mode=${3:-"balanced"} 
number_of_files=${4:-"50"}

# Navigate to the directory where the files and scripts are located
cd "${data_path}/${mode}/data"


# Loop through each set of files from 1 to 50
for j in $(seq 1 $number_of_files)
do
    # copy lab1 lab2 lab3 files to the corresponding folders
    cp ./${j}_lab_1_intensities_data.tsv ../lab1/intensities.tsv
    cp ./${j}_lab_1_design.tsv ../lab1/design.tsv

    cp ./${j}_lab_2_intensities_data.tsv ../lab2/intensities.tsv
    cp ./${j}_lab_2_design.tsv ../lab2/design.tsv

    cp ./${j}_lab_3_intensities_data.tsv ../lab3/intensities.tsv
    cp ./${j}_lab_3_design.tsv ../lab3/design.tsv

    # Run the R scripts for meta-analysis on each triple of files
    python /home/yuliya/repos/cosybio/FedProt/evaluation_utils/fedprot_prototype/fedprot_script.py \
    ${data_path}/ \
    ${mode} \
    lab1,lab2,lab3 \
    ${base_path}/

    # rename adj.P.Val and P.Value columns in the DPE.tsv file to sca.adj.pval and sca.P.Value
    # do not rename other columns
    sed -i 's/adj.P.Val/sca.adj.pval/g' ${base_path}/${mode}/results/DPE.csv
    sed -i 's/P.Value/sca.P.Value/g' ${base_path}/${mode}/results/DPE.csv

    # Copy the DPE.tsv file with the j prefix to a new location
    mv "${base_path}/${mode}/results/DPE.csv" "${base_path}/${mode}/results/${j}_DPE.tsv"
    
done

echo "FedProt run and file copying completed."
