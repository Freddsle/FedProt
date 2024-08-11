import os
import time
import zipfile
import logging

import featurecloud_test as fc_utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# save time of the start
start_time = time.time()

### SETTINGS
# The directory that contains ALL data needed in ANY experiment.

## CHANGE THIS TO THE DIRECTORY WHERE YOU HAVE THE REPO CLONED
data_dir = os.path.abspath("/home/yuliya/repos/cosybio/FedProt/")

### SETTING THE EXPERIMENTS
experiments = list()
result_file_names = list()

## Baterial dataset
bacterial_experiment = fc_utils.Experiment(
        clients=[os.path.join(data_dir, "data", "bacterial_data", "balanced", "lab_A"),
                 os.path.join(data_dir, "data", "bacterial_data", "balanced", "lab_B"),
                 os.path.join(data_dir, "data", "bacterial_data", "balanced", "lab_C"),
                 os.path.join(data_dir, "data", "bacterial_data", "balanced", "lab_D"),
                 os.path.join(data_dir, "data", "bacterial_data", "balanced", "lab_E"),
        ]
)
result_file_names.append(os.path.join(data_dir, "evaluation", "fc_results", "DEP_bacterial_FEDERATED.tsv"))
experiments.append(bacterial_experiment)

# # human serum
serum_experiment = fc_utils.Experiment(
        clients=[os.path.join(data_dir, "data", "TMT_data", "01_smaller_lib_balanced_PG_MajorPG", "Center1"),
                 os.path.join(data_dir, "data", "TMT_data", "01_smaller_lib_balanced_PG_MajorPG", "Center2"),
                 os.path.join(data_dir, "data", "TMT_data", "01_smaller_lib_balanced_PG_MajorPG", "Center3")
        ]
)
result_file_names.append(os.path.join(data_dir, "evaluation", "fc_results", "DEP_TMT_FEDERATED.tsv"))
experiments.append(serum_experiment)

# simulated balanced
simbal_experiment = fc_utils.Experiment(
        clients=[os.path.join(data_dir, "data", "simulated_data", "balanced", "lab1"),
                 os.path.join(data_dir, "data", "simulated_data", "balanced", "lab2"),
                 os.path.join(data_dir, "data", "simulated_data", "balanced", "lab3")
        ]
)
result_file_names.append(os.path.join(data_dir, "evaluation", "fc_results", "DEP_simulBal_FEDERATED.tsv"))
experiments.append(simbal_experiment)

# simulated mil imbalanced
simMimbal_experiment = fc_utils.Experiment(
        clients=[os.path.join(data_dir, "data", "simulated_data", "mild_imbalanced", "lab1"),
                 os.path.join(data_dir, "data", "simulated_data", "mild_imbalanced", "lab2"),
                 os.path.join(data_dir, "data", "simulated_data", "mild_imbalanced", "lab3")
        ]
)
result_file_names.append(os.path.join(data_dir, "evaluation", "fc_results", "DEP_simulMildImbal_FEDERATED.tsv"))
experiments.append(simMimbal_experiment)

# simulated balanced 10
simImbal_experiment = fc_utils.Experiment(
        clients=[os.path.join(data_dir, "data", "simulated_data", "imbalanced", "lab1"),
                 os.path.join(data_dir, "data", "simulated_data", "imbalanced", "lab2"),
                 os.path.join(data_dir, "data", "simulated_data", "imbalanced", "lab3")
        ]
)
result_file_names.append(os.path.join(data_dir, "evaluation", "fc_results", "DEP_simuImb_FEDERATED.tsv"))
experiments.append(simImbal_experiment)

### MAIN
# Starts the FeatureCloud controler and runs the experiments.
# Start the controller
try:
    fc_utils.start_controller(data_dir)
except Exception as e:
    raise RuntimeError(f"Experiment could not be started! Error: \n{e}") from e

# Log how long it took to start the controller
logging.info(f"Controller started in {int(time.time() - start_time)} seconds.")

# Run the experiments
for exp, result_filename in zip(experiments, result_file_names):
    start_time = time.time()

    ### RUN THE FEATURECLOUD EXPERIMENT AND EXTRACT INDIVIDUAL RESULTS
    logging.info(f"Starting experiment:\n{exp}")
    try:
        result_files_zipped = fc_utils.run_test(exp, data_dir)
    except Exception as e:
        logging.error(f"Experiment could not be started or aborted too many times! Error: \n{e}")
        logging.error("_______________¡FAILED_EXPERIMENT!_______________")
        continue
    logging.info(f"Experiment finished successfully! Result files: {result_files_zipped}")

    ### Postprocessing
    time.sleep(10)
    logging.info("Starting postprocessing...")
    result_folder = os.path.dirname(result_filename)
    logging.info(f"Results will be saved in: {result_folder}")

    # Extract one (first) file from the results
    idx = 0
    zipfilename = result_files_zipped[idx]
    if not os.path.exists(os.path.join(data_dir, "tests", "tests", zipfilename)):
        logging.info(f"Waiting for file {zipfilename} to be available...")
        time.sleep(5) 
    
    with zipfile.ZipFile(os.path.join(data_dir, "tests", "tests", zipfilename), 'r') as zip_ref:
        zip_ref.extract(f"DPE.csv", os.path.join(result_folder))
        # rename the just extracted file so it doesn't get overwritten
        os.rename(os.path.join(result_folder, "DPE.csv"), os.path.join(result_folder, result_filename))

    logging.info(f"Experiment finished in {int((time.time() - start_time) // 60)} min {int((time.time() - start_time) % 60)} sec.")
    logging.info("_______________¡SUCCESSFUL_EXPERIMENT!_______________")