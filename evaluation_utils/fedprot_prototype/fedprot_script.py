import logging
import sys
# for json files
import json

import pandas as pd
import numpy as np
import yaml

sys.path.append('/home/yuliya/repos/cosybio/FedProt/FedProt')
from client import Client
import utils

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

# Test mode - if True, turn off "less then 2 non-NA samples" filter - keeping all values, even if there is only one sample
# If False, filter out proteins with less then 2 non-NA samples - but the result will be different from the centralized version. 
TEST_MODE = True

try:
    sys_arguments = sys.argv
    number_of_arguments = len(sys_arguments)
    logging.info(f"Number of arguments: {number_of_arguments}")
    
    if number_of_arguments > 1:
        data_dir = str(sys_arguments[1])
        MODE = str(sys_arguments[2])
        cohorts = str(sys_arguments[3]).split(",")
        output_path = str(sys_arguments[4])
    
    # ONLY FOR EVALUATION
    if number_of_arguments > 5:
        keep_proteins_path = str(sys_arguments[5])
        CHECK_NUMBER = True
    else:
        CHECK_NUMBER = False

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
else:
    counts_path = None
result_name = config["fedprot"]["result_table"]

max_na_rate = config["fedprot"]["max_na_rate"]
log_transformed = config["fedprot"]["log_transformed"]

experiment_type = config["fedprot"]["experiment_type"]

if experiment_type == "TMT":
    ref_type = config["fedprot"]["ref_type"]
    plex_covariate = config["fedprot"]["plex_covariate"]
    plex_column = config["fedprot"]["plex_column"]

    use_median = config["fedprot"]["use_median"]
    use_irs = config["fedprot"]["use_irs"]
else:
    ref_type = None
    plex_column = None

remove_single_pep_protein = config["fedprot"]["remove_single_pep_protein"]
target_classes = config["fedprot"]["target_classes"]
covariates = config["fedprot"]["covariates"]

only_shared_proteins = config["fedprot"]["only_shared_proteins"]

logging.info(f"Data directory: {data_dir}")

#############################################################
# initialize the client
store_clients = {}

prot_names = []
if experiment_type == "TMT":
    if plex_covariate:
        plex_covariates_list = []

if use_counts:
    list_of_counts = []

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
        ref_type = ref_type,
        plex_column = plex_column,
        target_classes = target_classes,
        TEST_MODE = TEST_MODE
    )
    store_clients[client.cohort_name] = client

    # add plex covariates to the list
    if experiment_type == "TMT":
        if plex_covariate:
            plex_covariates_list += client.tmt_names

    if use_counts:
        min_counts = client.get_min_count()
        list_of_counts.append(min_counts.to_dict())
        logging.info("Min counts have been computed...")
    
    # add list (as list of lists) of protein names
    prot_names.append(client.prot_names)

if experiment_type == "TMT":
    if plex_covariate:
        plex_covariates_list = sorted(list(set(plex_covariates_list)))

# SERVER SIDE
prot_names = utils.get_analyzed_proteins(prot_names, only_shared_proteins)
logging.info(f"SERVER: Common proteins: {len(prot_names)}")
variables = target_classes + covariates

if use_counts:
    min_counts = list()
    for local_counts in list_of_counts:
        min_counts.append(pd.DataFrame.from_dict(local_counts, orient='index'))
    global_min_counts = pd.concat(min_counts, axis=1)
    global_min_counts = global_min_counts.min(axis=1)
    global_min_counts = global_min_counts[prot_names]
    logging.info("Min counts have been aggregated...")
    logging.info(f"Size of global min counts: {global_min_counts.shape}")

# CLIENT SIDE
# validate inputs and add cohort effects
# apply filters
list_of_na_counts_tuples = []

for c in cohorts:
    client = store_clients[c]
    if use_counts:
        client.counts = global_min_counts.copy()

    client.validate_inputs(prot_names, variables)

    # add cohort effect columns to each design matrix
    # if plex_covariate exists, use this column as a cohort effect
    if experiment_type == "TMT":
        if plex_covariate:
            client.add_cohort_effects_to_design(plex_covariates_list, plex_covariate)
    else:
        client.add_cohort_effects_to_design(cohorts)

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

logging.info("CLIENT: Filters have been applied...")
logging.info("Min counts have been computed...")


# SERVER SIDE
keep_proteins = utils.filter_features_na_rate(list_of_na_counts_tuples, max_na_rate)
logging.info(f"SERVER: Proteins after filtering: {len(keep_proteins)}")

# ONLY FOR EVALUATION
if CHECK_NUMBER:
    # read keep_proteins from json file 
    with open(keep_proteins_path, 'r') as file:
        list_from_prefilter = json.load(file)

    # add keep proteins to the json file
    list_from_prefilter[MODE]["FedProt"] = keep_proteins
    # write to the json file
    with open(keep_proteins_path, 'w') as file:
        json.dump(list_from_prefilter, file)
    logging.info(f"SERVER: Proteins after filtering FOR EVALUATION: {len(keep_proteins)}")

# CLIENT SIDE
for c in cohorts:
    client = store_clients[c]
    client.update_prot_names(keep_proteins)

# Normalize intensities
# CLIENT SIDE
# if TMT and use_median, compute medians
if experiment_type == "TMT" and use_median:
    mead_samples_tuples = []
    for c in cohorts:
        client = store_clients[c]
        avg_median = client.compute_medians()
        total_samples = len(client.sample_names)
        mead_samples_tuples.append((avg_median, total_samples))
        
# SERVER SIDE
# if TMT and use_median, compute medians
if experiment_type == "TMT" and use_median:
    global_median_mean = utils.aggregate_medians(mead_samples_tuples)

# CLIENT SIDE
# if TMT and use_median, compute medians
if experiment_type == "TMT" and use_median:
    for c in cohorts:
        client = store_clients[c]
        client.mean_median_centering(global_median_mean)

        if experiment_type == "TMT" and use_irs:
            client.irsNorm_in_silico_single_center()

#############################################################
# Mask preparation
# CLIENT SIDE
masks_list = []
for c in cohorts:
    client = store_clients[c]
    client.prepare_for_limma(keep_proteins)

    masks_list.append(client.get_mask())

# SERVER SIDE
variables = client.design.columns.values
k = len(variables)
n = len(keep_proteins)
logging.info(f"SERVER: Number of proteins: {n}, number of variables: {k}")

# aggregate mask
global_mask = utils.aggregate_masks(masks_list, n, k, used_SMPC=False)

# CLIENT SIDE
list_of_masks = []
for c in cohorts:
    client = store_clients[c]
    updated_masks = client.updated_mask(global_mask)
    list_of_masks.append(updated_masks)

#############################################################

# CLIENT SIDE
XtX_XtY_list = []

for c in cohorts:
    client = store_clients[c]

    logging.info("Start computing XtX and XtY")
    XtX, XtY = client.compute_XtX_XtY()
    XtX_XtY_list.append((XtX, XtY))

    logging.info(f"XtX and XtY have been computed for {client.cohort_name}")
    logging.info(f'Design colnames: {client.design.columns.values}')


# SERVER SIDE
mask_glob = utils.aggregate_masks(list_of_masks, n, k, second_round=True, used_SMPC=False)

XtX_glob, XtY_glob = utils.aggregate_XtX_XtY(XtX_XtY_list, n, k, used_SMPC=False)
logging.info("SERVER: XtX and XtY have been aggregated")
logging.info('SERVER: Computing beta and stdev')

beta, stdev_unscaled = utils.compute_beta_and_stdev(XtX_glob, XtY_glob, n, k, mask_glob)

# CLIENT SIDE
list_of_sse_cov_coef = []

for c in cohorts:
    client = store_clients[c]
    logging.info(f"Start computation of SSE and cov_coef...")
    SSE, cov_coef = client.compute_SSE_and_cov_coef(beta, mask_glob)
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
            cov_coef, SSE, Amean, n_measurements, n, mask_glob
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

# if not using counts - stop here and save the results
if not use_counts:
    logging.info("Writing results to file...")
    logging.info("File: " + f"{output_path}/{result_name}")
    output_file = f"{output_path}/{result_name}"
    results.to_csv(output_file, sep="\t")
    logging.info(f"Results have been saved to {output_file}")
    sys.exit(0)

# SERVER SIDE
if use_counts:
    global_min_counts = global_min_counts[keep_proteins] + 1
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
