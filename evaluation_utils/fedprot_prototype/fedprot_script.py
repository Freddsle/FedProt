import logging
import sys

import pandas as pd
import numpy as np
import yaml

sys.path.append('/home/yuliya/repos/cosybio/FedProt/FedProt')
from client import Client
import utils

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)


try:
    if len(sys.argv) > 1:
        data_dir = str(sys.argv[1])
        MODE = str(sys.argv[2])
        cohorts = str(sys.argv[3]).split(",")
        output_path = str(sys.argv[4])
    else:
        # test data
        data_dir = "/home/yuliya/repos/cosybio/FedProt/data/bacterial_data/"
        MODE = "balanced"
        cohorts = ["lab_A", "lab_B", "lab_C", "lab_D", "lab_E"]
        output_path = "/home/yuliya/repos/cosybio/FedProt/evaluation/bacterial/"
except ValueError:
    print("Invalid input.")
    sys.exit(1)

data_dir = f"{data_dir}/{MODE}"  # path to data folder
output_path = f"{output_path}/{MODE}/results"  # path to output folder

# read yml file with config from first cohort
# config = pd.read_csv(f"{data_dir}/{cohorts[0]}/config.yml", sep=":", header=None)

with open(f"{data_dir}/{cohorts[0]}/config.yml", 'r') as file:
    config = yaml.safe_load(file)

# parse config
intensity_path = config["fedprot"]["intensities"]
design_path = config["fedprot"]["design"]
use_counts = config["fedprot"]["use_counts"]
if use_counts:
    counts_path = config["fedprot"]["counts"]
result_name = config["fedprot"]["result_table"]

max_na_rate = config["fedprot"]["max_na_rate"]
log_transformed = config["fedprot"]["log_transformed"]

experiment_type = config["fedprot"]["experiment_type"]
remove_single_pep_protein = config["fedprot"]["remove_single_pep_protein"]
target_classes = config["fedprot"]["target_classes"]
covariates = config["fedprot"]["covariates"]


#############################################################
# initialize the client
store_clients = {}

prot_names = []

# clinets are joining
for cohort_name in cohorts:
    intensity_file_path = f"{data_dir}/{cohort_name}/{intensity_path}"
    design_file_path = f"{data_dir}/{cohort_name}/{design_path}"
    if use_counts:
        count_file_path = f"{data_dir}/{cohort_name}/{counts_path}"
    else:
        count_file_path = None
    
    # Initialize the client
    client = Client(
        cohort_name,
        intensity_file_path,
        count_file_path,
        design_file_path,
        experiment_type,
        log_transformed = log_transformed,
    )
    store_clients[client.cohort_name] = client

    # add list (as list of lists) of protein names
    prot_names.append(client.prot_names)

# SERVER SIDE
prot_names = utils.get_common_proteins(prot_names)
logging.info(f"SERVER: Common proteins: {len(prot_names)}")
variables = target_classes + covariates

# CLIENT SIDE
# validate inputs and add cohort effects
# apply filters

list_of_na_counts_tuples = []

for c in cohorts:
    client = store_clients[c]
    client.validate_inputs(prot_names, variables)

    # add cohort effect columns to each design matrix
    # add 1 column less than the number of cohorts
    client.add_cohort_effects_to_design(cohorts[1:])

    logging.info(f"Samples in {client.cohort_name} data: {len(client.sample_names)}")        
    logging.info(f"Protein groups in {client.cohort_name} data:  {len(client.prot_names)}")
    logging.info(f"Design {client.cohort_name} has beed updated with cohort effects")

    na_count_in_variable, samples_per_class = client.apply_filters(
            min_f=max_na_rate, 
            remove_single_peptide_prots=remove_single_pep_protein
        )
    # add both as a tuple to the list
    list_of_na_counts_tuples.append((
        na_count_in_variable.to_dict(orient='index'), samples_per_class))
    
# SERVER SIDE
keep_proteins = utils.filter_features_na_rate(list_of_na_counts_tuples, max_na_rate)
logging.info(f"SERVER: Proteins after filtering: {len(keep_proteins)}")

# CLIENT SIDE
XtX_XtY_list = []

for c in cohorts:
    client = store_clients[c]
    client.update_prot_names(keep_proteins)

    logging.info("Start computing XtX and XtY")
    client.prepare_for_limma(keep_proteins)
    XtX, XtY = client.compute_XtX_XtY()
    XtX_XtY_list.append((XtX, XtY))

    logging.info(f"XtX and XtY have been computed for {client.cohort_name}")
    logging.info(f'Design colnames: {client.design.columns.values}')

variables = client.design.columns.values

# SERVER SIDE
k = len(variables)
n = len(keep_proteins)
logging.info(f"SERVER: Number of proteins: {n}, number of variables: {k}")

XtX_glob, XtY_glob = utils.aggregate_XtX_XtY(XtX_XtY_list, n, k, used_SMPC=False)
logging.info("SERVER: XtX and XtY have been aggregated")
logging.info('SERVER: Computing beta and stdev')

beta, stdev_unscaled = utils.compute_beta_and_stdev(XtX_glob, XtY_glob, n, k)

# CLIENT SIDE
list_of_sse_cov_coef = []

for c in cohorts:
    client = store_clients[c]
    logging.info(f"Start computation of SSE and cov_coef...")
    SSE, cov_coef = client.compute_SSE_and_cov_coef(beta)
    intensities_sum = np.array(client.sum_intensities())
    n_measurements = np.array(client.get_not_na())
    list_of_sse_cov_coef.append((SSE, cov_coef, intensities_sum, n_measurements))
    logging.info(f"SSE and cov_coef have been computed for {client.cohort_name}")

# SERVER SIDE
logging.info("SERVER: Aggregating SSE and cov_coef...")
SSE, cov_coef, Amean, n_measurements = utils.aggregate_SSE_and_cov_coef(
        list_of_sse_cov_coef, n, k, False, len(store_clients))
logging.info('Aggregation of SSE is done, start computing global parameters...')
sigma, cov_coef, df_residual, Amean, var = utils.compute_SSE_and_cov_coef_global(
            cov_coef, SSE, Amean, n_measurements, n, k
        )

logging.info('Making contrasts...')
contrasts = utils.make_contrasts(target_classes, variables)

logging.info('Fit contrasts...')
contrast_matrix = contrasts.values
ncoef = cov_coef.shape[1]

beta, stdev_unscaled, cov_coef = utils.fit_contrasts(
    beta, 
    contrast_matrix,
    cov_coef, 
    ncoef, 
    stdev_unscaled
)

logging.info('Start eBays stage...')
results, df_total = utils.moderatedT(
    var, 
    df_residual,
    beta, 
    stdev_unscaled
)
results["AveExpr"] = Amean
logging.info('"Result table is pre-computed...')

logging.info("Computing B statistic...")
results = utils.Bstat(
    df_total, 
    stdev_unscaled, 
    results,
    stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)
logging.info("B statistic has been computed...")

logging.info("Computing p-values...")
results = utils.topTableT( 
    results,
    keep_proteins, 
    beta, 
    stdev_unscaled, 
    df_total,
    adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
logging.info("P-values have been computed...")

# CLIENT SIDE
list_of_counts = []

for c in cohorts:
    client = store_clients[c]
    min_counts = client.get_min_count()
    list_of_counts.append(min_counts.to_dict())
logging.info("Min counts have been computed...")

# SERVER SIDE
if use_counts:
    min_counts = list()
    for local_counts in list_of_counts:
        min_counts.append(pd.DataFrame.from_dict(local_counts, orient='index'))
    global_min_counts = pd.concat(min_counts, axis=1)
    logging.info("Min counts have been aggregated...")
    global_min_counts = global_min_counts.min(axis=1).loc[keep_proteins] + 1
    logging.info(f"Size of global min counts: {global_min_counts.shape}")

    logging.info("Start spectral count eBayes stage...")
    results = utils.spectral_count_ebayes(
        results, 
        global_min_counts, 
        keep_proteins,
        beta, 
        stdev_unscaled,
        df_residual, 
        sigma,
        return_sorted=False, fit_method="loess")
    logging.info("Spectral count eBayes stage is done...")

logging.info("Writing results...")
output_file = f"{output_path}/{result_name}"
results.to_csv(output_file, sep="\t")
logging.info(f"Results have been saved to {output_file}")