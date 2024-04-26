# Run DIA-NN

1. Pull docker image  
    `docker pull biocontainers/diann:1.8.1_cv2`
2. Using escherichia_coli_MG1655.fasta
3. For lab B - unpack:  
    `nohup bash -c 'for zipfile in ./raw_zip/*.zip; do dir="./raw/${zipfile##*/}"; dir="${dir%.zip}"; mkdir -p "$dir" && unzip -d "$dir" "$zipfile"; done' &`

4. DIA-NN run for labs: 
   1. lab A - script for lab A.  
        New run (updated DIA-NN version). Copy mzML file to new dir. Run:
        ```
        nohup docker run --rm \
        -v /cosybio/project/burankova/FedProt_lab_storage/labA_results/:/data biocontainers/diann:1.8.1_cv2 diann \
        --dir "/data/raw_diann/" \
        --lib "" \
        --threads 10 \
        --verbose 1 \
        --out "/data/2024-04-20_21-25-12_diann_run/report.tsv" \
        --qvalue 0.01 --matrices \
        --out-lib "/data/2024-04-20_21-25-12_diann_run/report-lib.tsv" \
        --gen-spec-lib --predictor \
        --fasta "/data/escherichia_coli_MG1655.fasta" \
        --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 360 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --smart-profiling  > diann_output.log 2>&1 &
        ```

   2. lab B - "When running on Linux (native builds, not Wine), only .d, .mzML, and .dia data are supported." (https://github.com/vdemichev/DiaNN). Run (`chmod -R 777 /cosybio/project/burankova/FedProt_lab_storage/labB_results/2024-04-20_21-56-12_diann_run/`):
        ```
        nohup docker run --rm \
        -v /cosybio/project/burankova/FedProt_lab_storage/labB_results/:/data biocontainers/diann:1.8.1_cv2 diann \
        --f /data/raw/Clinspect_E_coli_A_QC1_Slot1-6_1_8640.d --f /data/raw/Clinspect_E_coli_A_S34_Slot1-18_1_8666.d \
        --f /data/raw/Clinspect_E_coli_B_66_Slot1-13_1_8647.d --f /data/raw/Clinspect_E_coli_B_S78_Slot1-17_1_8665.d \
        --f /data/raw/Clinspect_E_coli_A_QC2_Slot1-15_1_8649.d --f /data/raw/Clinspect_E_coli_A_S39_Slot1-26_1_8675.d \
        --f /data/raw/Clinspect_E_coli_B_QC3_Slot1-23_1_8672.d --f /data/raw/Clinspect_E_coli_B_S85_Slot1-12_1_8646.d \
        --f /data/raw/Clinspect_E_coli_A_S14_Slot1-11_1_8645.d --f /data/raw/Clinspect_E_coli_A_S42_Slot1-21_1_8670.d \
        --f /data/raw/Clinspect_E_coli_B_QC4_Slot1-10_1_8644.d --f /data/raw/Clinspect_E_coli_B_S86_Slot1-4_1_8638.d \
        --f /data/raw/Clinspect_E_coli_A_S20_Slot1-20_1_8669.d --f /data/raw/Clinspect_E_coli_A_S47_Slot1-14_1_8648.d  \
        --f /data/raw/Clinspect_E_coli_B_S53_Slot1-22_1_8671.d --f /data/raw/Clinspect_E_coli_B_S91_Slot1-25_1_8674.d \
        --f /data/raw/Clinspect_E_coli_A_S23_Slot1-8_1_8642.d --f /data/raw/Clinspect_E_coli_A_S5_Slot1-7_1_8641.d \
        --f /data/raw/Clinspect_E_coli_B_S58_Slot1-5_1_8639.d --f /data/raw/Clinspect_E_coli_B_S98_Slot1-9_1_8643.d \
        --f /data/raw/Clinspect_E_coli_A_S29_Slot1-19_1_8667.d --f /data/raw/Clinspect_E_coli_A_S7_Slot1-24_1_8673.d \
        --f /data/raw/Clinspect_E_coli_B_S62_Slot1-16_1_8664.d \
        --lib "" \
        --threads 40 \
        --verbose 1 \
        --out "/data/2024-04-25_15-05-38_diann_run/report.tsv" \
        --out-lib "/data/2024-04-25_15-05-38_diann_run/lib.tsv" \
        --qvalue 0.01 --matrices \
        --fasta "/data/escherichia_coli_MG1655.fasta" \
        --predictor --gen-spec-lib  \
        --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 360 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --smart-profiling > diann_output.log 2>&1 &
        ```
        Note: raw data folder should have writable for docker - to save quant files.

   3. lab C - script for lab C.  
        New run (updated DIA-NN version):  
        ```
        nohup docker run  --rm \
        -v /cosybio/project/burankova/FedProt_lab_storage/labC_results/:/data biocontainers/diann:1.8.1_cv2 diann \
        --dir "/data/raw_diann/" \
        --lib "" \
        --threads 3 \
        --verbose 1 \
        --out "/data/2024-04-20_21-41-34_diann_run/report.tsv" \
        --qvalue 0.01 --matrices \
        --out-lib "/data/2024-04-20_21-41-34_diann_run/report-lib.tsv" \
        --gen-spec-lib --predictor \
        --fasta "/data/escherichia_coli_MG1655.fasta" \
        --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 360 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --smart-profiling > diann_output.log 2>&1 &
        ```

   4. lab D - same as for lab C and lab A.  
        New run (updated DIA-NN):
        ```
        nohup docker run --rm \
        -v /cosybio/project/burankova/FedProt_lab_storage/labD_results/:/data biocontainers/diann:1.8.1_cv2 diann \
        --dir "/data/raw_diann/" \
        --lib "" \
        --threads 2 \
        --verbose 1 \
        --out "/data/2024-04-20_21-45-56_diann_run/report.tsv" \
        --qvalue 0.01 --matrices \
        --out-lib "/data/2024-04-20_21-45-56_diann_run/report-lib.tsv" \
        --gen-spec-lib --predictor \
        --fasta "/data/escherichia_coli_MG1655.fasta" \
        --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 360 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --smart-profiling > diann_output.log 2>&1 &
        ```

   5. lab E - same as for lab C and lab A.  
        New run (updated DIA-NN):
        ```
        nohup docker run --rm \
        -v /cosybio/project/burankova/FedProt_lab_storage/labE_results/:/data biocontainers/diann:1.8.1_cv2 diann \
        --dir "/data/raw_diann/" \
        --lib "" \
        --threads 2 \
        --verbose 1 \
        --out "/data/2024-04-20_21-51-06_diann_run/report.tsv" \
        --qvalue 0.01 --matrices \
        --out-lib "/data/2024-04-20_21-51-06_diann_run/report-lib.tsv" \
        --gen-spec-lib --predictor \
        --fasta "/data/escherichia_coli_MG1655.fasta" \
        --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 360 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --smart-profiling > diann_output.log 2>&1 &
        ```



# OLD scripts (used to get mzML files):

1. Script for Lab A:
```
#!/bin/bash

# Check for required arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_dir> <threads> <fasta_file>"
  exit 1
fi

# Set shared variables
input_dir="$1"
threads="$2"
fasta_file="$3"
newdir="$(date +%Y-%m-%d_%H-%M-%S)_diann_run"

# Create folders on the host system
mkdir -m 777 "$input_dir/$newdir"
mkdir -m 777 "$input_dir/$newdir/mzXML"

# Common options for both containers
docker_opts="--rm --workdir /data -e WINEDEBUG=-all -v $input_dir:/data"

# First container
container1=$(docker run -d $docker_opts chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest \
bash -c "cd ./raw; find . -maxdepth 1 -type f -name \"*.raw\" | xargs -I {} -P $threads wine msconvert {} -o \"/data/$newdir/mzXML\"")

# Wait for the first container to finish
docker wait $container1

# Second container
docker run -d $docker_opts biocontainers/diann:v1.8.1_cv1 diann \
--dir "/data/$newdir/mzXML/" \
--lib "" \
--threads $threads \
--verbose 1 \
--out "/data/$newdir/report.tsv" \
--qvalue 0.01 \
--matrices \
--out-lib "/data/$newdir/report-lib.tsv" \
--gen-spec-lib \
--predictor \
--fasta "/data/$fasta_file" \
--fasta-search \
--min-fr-mz 200 \
--max-fr-mz 1800 \
--met-excision \
--cut K*,R* \
--missed-cleavages 2 \
--min-pep-len 7 \
--max-pep-len 30 \
--min-pr-mz 360 \
--max-pr-mz 1800 \
--min-pr-charge 1 \
--max-pr-charge 4 \
--unimod4 \
--var-mods 1 \
--var-mod UniMod:35,15.994915,M \
--var-mod UniMod:1,42.010565,*n \
--monitor-mod UniMod:1 \
--reanalyse \
--smart-profiling 
```


2. Script for lab C:
```
#!/bin/bash

# Check for required arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_dir> <threads> <fasta_file>"
  exit 1
fi

# Set shared variables
input_dir="$1"
threads="$2"
fasta_file="$3"
newdir="$(date +%Y-%m-%d_%H-%M-%S)_diann_run"

# Create folders on the host system
mkdir -m 777 "$input_dir/$newdir"
mkdir -m 777 "$input_dir/$newdir/mzXML"

# Common options for both containers
docker_opts="--rm --workdir /data -e WINEDEBUG=-all -v $input_dir:/data"

# First container
container1=$(docker run -d $docker_opts chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest \
bash -c "cd ./raw; find . -maxdepth 1 -type f -name \"*.raw\" | xargs -I {} -P $threads wine msconvert {} -o \"/data/$newdir/mzXML\"")

# Wait for the first container to finish
docker wait $container1

# Second container
docker run -d  $docker_opts biocontainers/diann:v1.8.0_cv1 diann \
--dir "/data/$newdir/mzXML/" \
--lib "/data/labC_results/DIANN/lib.tsv.speclib" \
--threads $threads \
--verbose 1 \
--out "/data/$newdir/report.tsv" \
--qvalue 0.01 \
--matrices \
--fasta "/data/$fasta_file" \
--met-excision \
--cut K*,R* \
--unimod4 \
--var-mods 2 \
--var-mod UniMod:35,15.994915,M \
--var-mod UniMod:1,42.010565,*n \
--monitor-mod UniMod:1 \
--reanalyse \
--smart-profiling
```







