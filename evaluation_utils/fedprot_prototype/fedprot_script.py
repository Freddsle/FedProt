import logging
import sys
# for json files
import json

import pandas as pd
import numpy as np
import yaml

# CHANGE THIS TO THE DIRECTORY WHERE YOU HAVE THE REPO CLONED
sys.path.append('/home/yuliya/repos/cosybio/FedProt/FedProt')

from client import Client
import utils

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

# Test mode - if True, turn off "less then 2 non-NA samples" filter - keeping all values, even if there is only one sample
# If False, filter out proteins with less then 2 non-NA samples - but the result will be different from the centralized version. 

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
with open(f"{data_dir}/{cohorts[0]}/config.yml", 'r') as file:
    config = yaml.safe_load(file)

# parse config
intensity_path = config["fedprot"]["intensities"]
design_path = config["fedprot"]["design"]
use_counts = config["fedprot"]["use_counts"]


if "TEST_MODE" in config["fedprot"]:
    TEST_MODE = config["fedprot"]["TEST_MODE"]
else:
    TEST_MODE = False

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
remove_single_value_design = config["fedprot"]['remove_single_value_design']
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
        TEST_MODE = TEST_MODE,
        remove_single_value_design = remove_single_value_design
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
    logging.info(f"Type check after aggregation: global_min_counts type: {type(global_min_counts)}")

# CLIENT SIDE
# validate inputs and add cohort effects
# apply filters
list_of_na_counts_tuples = []

for c in cohorts:
    client = store_clients[c]
    if use_counts:
        client.counts = global_min_counts.copy()
        logging.info(f"Type check after updating counts: client.counts type: {type(client.counts)}")

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
    logging.debug(f"Type check before adding to list: na_count_in_variable type: {type(na_count_in_variable)}, samples_per_class type: {type(samples_per_class)}")
    # add both as a tuple to the list
    list_of_na_counts_tuples.append((
        na_count_in_variable.to_dict(orient='index'), samples_per_class))

logging.info("CLIENT: Filters have been applied...")
logging.debug(f"Type check after filtering: list_of_na_counts_tuples type: {type(list_of_na_counts_tuples)}")

# SERVER SIDE
logging.debug(f"NA counts are received, length: {len(list_of_na_counts_tuples)}")
keep_proteins = utils.filter_features_na_rate(list_of_na_counts_tuples, max_na_rate)
logging.info(f"SERVER: Proteins after filtering: {len(keep_proteins)}")
logging.debug(f"Type check after filtering: keep_proteins type: {type(keep_proteins)}")

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
    logging.debug(f"SERVER: Proteins after filtering FOR EVALUATION: {len(keep_proteins)}")

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
        logging.debug(f"Type check before adding to tuple: avg_median type: {type(avg_median)}, total_samples type: {type(total_samples)}")
        mead_samples_tuples.append((avg_median, total_samples))
        
# SERVER SIDE
# if TMT and use_median, compute medians
if experiment_type == "TMT" and use_median:
    global_median_mean = utils.aggregate_medians(mead_samples_tuples)
    logging.debug(f"Type check after aggregation: global_median_mean type: {type(global_median_mean)}")

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
for i, c in enumerate(cohorts):
    client = store_clients[c]
    client.prepare_for_limma(keep_proteins)

    mask = client.get_mask()
    logging.info(f"Type check before adding to list: mask type: {type(mask)}")

    masks_list.append(mask)

# SERVER SIDE
variables = client.design.columns.values
k = len(variables)
n = len(keep_proteins)
logging.info(f"SERVER: Number of proteins: {n}, number of variables: {k}")


# aggregate mask
global_mask = utils.aggregate_masks(masks_list, n, k, used_SMPC=False)
logging.debug(f"Type check after aggregation: global_mask type: {type(global_mask)}")

# # save mask to the file
# output_file = f"{output_path}/{result_name}_mask1.tsv"
# pd.DataFrame(global_mask).to_csv(output_file, sep="\t")

global_mask = utils.process_mask(global_mask,
                                 client_number=len(store_clients))
logging.debug(f"Type check after processing: global_mask type: {type(global_mask)}")

# output_file = f"{output_path}/{result_name}_mask1_1.tsv"
# pd.DataFrame(global_mask).to_csv(output_file, sep="\t")


# CLIENT SIDE
list_of_masks = []
for c in cohorts:
    updated_masks = None
    client = store_clients[c]    
    updated_masks = client.updated_mask(global_mask.copy())
    logging.info(f"Type check after updating: updated_masks type: {type(updated_masks)}")
    logging.debug(f"Updated mask mean: {np.mean(updated_masks)})")
    list_of_masks.append(updated_masks)

#############################################################

# CLIENT SIDE
XtX_XtY_list = []

for c in cohorts:
    client = store_clients[c]
    XtX, XtY = None, None
    logging.info("Start computing XtX and XtY")
    XtX, XtY = client.compute_XtX_XtY()
    logging.info(f"Type check before adding to list: XtX type: {type(XtX)}, XtY type: {type(XtY)}")
    XtX_XtY_list.append((XtX, XtY))

    logging.info(f"XtX and XtY have been computed for {client.cohort_name}")
    logging.info(f'Design colnames: {client.design.columns.values}')

# SERVER SIDE
mask_glob = utils.aggregate_masks(list_of_masks, n, k, used_SMPC=False)
logging.debug(f"Type check after aggregation: mask_glob type: {type(mask_glob)}")

# # save mask to the file
# output_file = f"{output_path}/{result_name}_mask2.tsv"
# pd.DataFrame(mask_glob).to_csv(output_file, sep="\t")

mask_glob = utils.process_mask(mask_glob, second_round=True)

# # save mask to the file
# output_file = f"{output_path}/{result_name}_mask2_1.tsv"
# pd.DataFrame(mask_glob).to_csv(output_file, sep="\t")

XtX_glob, XtY_glob = utils.aggregate_XtX_XtY(XtX_XtY_list, n, k, used_SMPC=False)
logging.info("SERVER: XtX and XtY have been aggregated")
logging.debug(f"Type check after aggregation: XtX_glob type: {type(XtX_glob)}, XtY_glob type: {type(XtY_glob)}")
logging.info('SERVER: Computing beta and stdev')

beta, stdev_unscaled = utils.compute_beta_and_stdev(XtX_glob, XtY_glob, n, k, mask_glob)
logging.info(f"Type check after computing: beta type: {type(beta)}, stdev_unscaled type: {type(stdev_unscaled)}")

# # write beta and stdev_unscaled to the file
# output_file = f"{output_path}/{result_name}_beta.tsv"
# pd.DataFrame(beta).to_csv(output_file, sep="\t")
# logging.info(f"Beta has been saved to {output_file}")
# output_file = f"{output_path}/{result_name}_stdev_unscaled.tsv"
# pd.DataFrame(stdev_unscaled).to_csv(output_file, sep="\t")


# CLIENT SIDE
list_of_sse_cov_coef = []

for c in cohorts:
    client = store_clients[c]
    SSE, cov_coef = None, None
    intensities_sum, n_measurements = None, None
    
    logging.info(f"Start computation of SSE and cov_coef...")
    SSE, cov_coef = client.compute_SSE_and_cov_coef(beta, mask_glob)
    intensities_sum = np.array(client.sum_intensities())
    n_measurements = np.array(client.get_not_na())
    logging.info(f"Type check before adding to list: SSE type: {type(SSE)}, cov_coef type: {type(cov_coef)}, intensities_sum type: {type(intensities_sum)}, n_measurements type: {type(n_measurements)}")
    list_of_sse_cov_coef.append((SSE, cov_coef, intensities_sum, n_measurements))
    logging.info(f"SSE and cov_coef have been computed for {client.cohort_name}")

# SERVER SIDE
SSE, cov_coef, Amean, n_measurements = None, None, None, None
logging.info("SERVER: Aggregating SSE and cov_coef...")
SSE, cov_coef, Amean, n_measurements = utils.aggregate_SSE_and_cov_coef(
        list_of_sse_cov_coef, n, k, False, len(store_clients))
logging.info(f"Type check after aggregation: SSE type: {type(SSE)}, cov_coef type: {type(cov_coef)}, Amean type: {type(Amean)}, n_measurements type: {type(n_measurements)}")

# # save SSE, cov_coef, Amean, n_measurements to the file
# output_file = f"{output_path}/{result_name}_SSE.tsv"
# pd.DataFrame(SSE).to_csv(output_file, sep="\t")
# logging.info(f"SSE has been saved to {output_file}")
# output_file = f"{output_path}/{result_name}_cov_coef.tsv"
# pd.DataFrame(cov_coef).to_csv(output_file, sep="\t")
# logging.info(f"cov_coef has been saved to {output_file}")
# output_file = f"{output_path}/{result_name}_Amean.tsv"
# pd.DataFrame(Amean).to_csv(output_file, sep="\t")
# logging.info(f"Amean has been saved to {output_file}")
# output_file = f"{output_path}/{result_name}_n_measurements.tsv"
# pd.DataFrame(n_measurements).to_csv(output_file, sep="\t")


logging.info('Aggregation of SSE is done, start computing global parameters...')
sigma, cov_coef, df_residual, Amean, var = utils.compute_SSE_and_cov_coef_global(
            cov_coef, SSE, Amean, n_measurements, n, mask_glob)
logging.info(f"Type check after computing global parameters: sigma type: {type(sigma)}, cov_coef type: {type(cov_coef)}, df_residual type: {type(df_residual)}, Amean type: {type(Amean)}, var type: {type(var)}")

logging.info('Making contrasts...')
contrasts = utils.make_contrasts(target_classes, variables)
logging.info(f"Type check after making contrasts: contrasts type: {type(contrasts)}")

logging.info('Fit contrasts...')
contrast_matrix = contrasts.values
ncoef = cov_coef.shape[1]
logging.info(f"Type check before fitting contrasts: contrast_matrix type: {type(contrast_matrix)}, ncoef type: {type(ncoef)}")

beta, stdev_unscaled, cov_coef = utils.fit_contrasts(
    beta, 
    contrast_matrix,
    cov_coef, 
    ncoef, 
    stdev_unscaled
)
logging.info(f"Type check after fitting contrasts: beta type: {type(beta)}, stdev_unscaled type: {type(stdev_unscaled)}, cov_coef type: {type(cov_coef)}")

logging.info('Start eBays stage...')
results, df_total = utils.moderatedT(
    var, 
    df_residual,
    beta, 
    stdev_unscaled
)
results["AveExpr"] = Amean
logging.info(f"Type check after computing moderated t-statistics: results type: {type(results)}, df_total type: {type(df_total)}")
logging.info('"Result table is pre-computed...')

logging.info("Computing B statistic...")
results = utils.Bstat(
    df_total, 
    stdev_unscaled, 
    results,
    stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)
logging.info(f"Type check after computing B statistic: results type: {type(results)}")

logging.info("Computing p-values...")
results = utils.topTableT( 
    results,
    keep_proteins, 
    beta, 
    stdev_unscaled, 
    df_total,
    adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
logging.info(f"Type check after computing p-values: results type: {type(results)}")

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
    logging.info(f"Type check after updating global min counts: global_min_counts type: {type(global_min_counts)}")
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
    logging.info(f"Type check after spectral count eBayes: results type: {type(results)}")
    logging.info("Spectral count eBayes stage is done...")

logging.info("Writing results...")
output_file = f"{output_path}/{result_name}"
results.to_csv(output_file, sep="\t")
logging.info(f"Results have been saved to {output_file}")