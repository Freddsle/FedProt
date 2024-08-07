from __future__ import annotations
from FeatureCloud.app.engine.app import AppState, app_state, Role
import time
import os
import bios

import pandas as pd
import numpy as np

from client import Client
import utils

CONFIG_FILE = "config.yml"

APP_DIR = "/app"

MOUNT_DIR = "/mnt/input"
OUTPUT_DIR = "/mnt/output"


@app_state(name='initial', role=Role.BOTH)
class InitialState(AppState):
    def register(self):
        self.register_transition('get_counts', Role.COORDINATOR)
        self.register_transition('validation', Role.PARTICIPANT)

    def run(self):
        # read config
        self.read_config()

        if self.load('use_smpc'):
            self.configure_smpc()
        # defining the client
        cohort_name = self.id

        if self.load('use_counts'):
            counts_path = os.path.join(MOUNT_DIR, self.load('counts_name'))
        else:
            counts_path = None
        
        # # initializing the client includes loading and preprocessing of the data
        client = Client(
            cohort_name,
            intensities_file_path = os.path.join(MOUNT_DIR, self.load('intensities_name')),
            count_file_path = counts_path,
            annotation_file_path = os.path.join(MOUNT_DIR, self.load('design_name')),
            experiment_type = self.load('experiment_type'),
            log_transformed = self.load('log_transformed'),
            ref_type = self.load('ref_type'),
            plex_column = self.load('plex_column'),
            target_classes = self.load('target_classes')
        ) 
        self.log(f"Data read! {client.intensities.shape[1]} samples, {client.intensities.shape[0]} proteins")

        # store client
        self.store(key='client', value=client)

        # get client's counts
        if self.load('use_counts'):
            counts = client.get_min_count()
            self.log("Counts are computed...")

        # send list of protein names (genes) to coordinator
        if self.load('use_counts'):
            dict_counts = counts.to_dict()
        else:
            dict_counts = None

        if self.load('experiment_type') == "TMT":
                self.log("[initial] Sending the counts, protein names list and plexes to the coordinator")
                self.send_data_to_coordinator((dict_counts, client.prot_names, client.tmt_names),
                                            send_to_self=True, use_smpc=False)
        else:
                self.log("[initial] Sending the counts and protein names list to the coordinator")
                self.send_data_to_coordinator((dict_counts, client.prot_names),
                                            send_to_self=True, use_smpc=False)

        # initialize app
        if self.is_coordinator:
            self.log("Transition to the next step...")
            return 'get_counts'
        self.log("Transition to validation...")
        return 'validation'

    def read_config(self):
        self.log("Reading config...")
        config = bios.read(os.path.join(MOUNT_DIR, CONFIG_FILE))
        config = config["fedprot"]

        self.store(key='intensities_name', value=config['intensities'])
        self.store(key='design_name', value=config['design'])
        self.store(key='use_counts', value=config['use_counts'])
        if config['use_counts']:
            self.store(key='counts_name', value=config['counts'])

        self.store(key='sep', value=config['sep'])
        self.store(key='result_table', value=config['result_table'])

        self.store(key='max_na_rate', value=config['max_na_rate'])
        self.store(key='log_transformed', value=config['log_transformed'])

        self.store(key='experiment_type', value=config['experiment_type'])

        if config['experiment_type'] == "TMT":
            self.store(key='ref_type', value=config['ref_type'])
            self.store(key='plex_covariate', value=config['plex_covariate'])
            self.store(key='plex_column', value=config['plex_column'])

            self.store(key='use_median', value=config['use_median'])
            self.store(key='use_irs', value=config['use_irs'])
        else:
            self.store(key='ref_type', value=None)
            self.store(key='plex_column', value=None)
            self.store(key='use_median', value=False)
            self.store(key='use_irs', value=False)

        self.store(key='remove_single_pep_protein', value=config['remove_single_pep_protein'])
        self.store(key='target_classes', value=config['target_classes'])
        self.store(key='covariates', value=config['covariates'])

        self.store(key='only_shared_proteins', value=config['only_shared_proteins'])

        self.store(key='use_smpc', value=config['use_smpc'])
        self.log("Config read!")
        self.log(f"Experiment type: {self.load('experiment_type')}; Shared proteins only: {self.load('only_shared_proteins')}")


@app_state(name='get_counts')
class GetCountsProteinsPlexesState(AppState):
    """
    Collects counts and protein lists from all clients.
    Additionally, if experiment type is "TMT" - collect plex names.
    """
    
    def register(self):
        self.register_transition('validation', role=Role.COORDINATOR)
    
    def run(self):
        self.log("[get_counts] Start collecting counts and protein lists...")
        experiment_type = self.load('experiment_type')
        plex_covariate = self.load('plex_covariate')
        
        if experiment_type == "TMT" and plex_covariate:
            list_of_counts, list_of_proteins, list_of_plexes = self.collect_tmt_data()
        else:
            list_of_counts, list_of_proteins = self.collect_non_tmt_data()
        
        self.log("[get_counts] Gathering done!")
        client_list = self.update_client_list()
        self.store(key='client_list', value=client_list)
        
        if experiment_type == "TMT" and plex_covariate:
            plex_covariates_list = self.aggregate_plex_covariates(list_of_plexes)
            self.store(key='plex_covariates_list', value=plex_covariates_list)            

        prot_names = self.aggregate_protein_groups(list_of_proteins)
        self.store(key='protein_groups', value=prot_names)

        if self.load('use_counts'):
            global_min_counts = self.aggregate_counts(list_of_counts, prot_names)
            self.store(key='min_counts', value=global_min_counts)
            counts_index = global_min_counts.index.tolist()
            counts_values = global_min_counts.values.tolist()
        else:
            counts_index = None
            counts_values = None
        
        variables = self.load_and_store_variables()
        self.store(key='variables', value=variables)

        self.log("Transition to validation ...")
        if experiment_type == "TMT" and plex_covariate:
            self.broadcast_data((client_list, 
                                 prot_names, 
                                 counts_index, counts_values, 
                                 variables, 
                                 plex_covariates_list), 
                                send_to_self=True, memo="ProteinsCounts")
        else:
            self.broadcast_data((client_list, 
                                 prot_names, 
                                 counts_index, counts_values, 
                                 variables), 
                                send_to_self=True, memo="ProteinsCounts")
        return 'validation'
    
    def collect_tmt_data(self):
        list_of_counts_proteins_plexes = self.gather_data(is_json=False, use_smpc=False)
        list_of_counts = [counts for counts, _, _ in list_of_counts_proteins_plexes]
        list_of_proteins = [proteins for _, proteins, _ in list_of_counts_proteins_plexes]
        list_of_plexes = [plexes for _, _, plexes in list_of_counts_proteins_plexes]
        return list_of_counts, list_of_proteins, list_of_plexes
    
    def collect_non_tmt_data(self):
        list_of_counts_proteins = self.gather_data(is_json=False, use_smpc=False)
        list_of_counts = [counts for counts, _ in list_of_counts_proteins]
        list_of_proteins = [proteins for _, proteins in list_of_counts_proteins]
        return list_of_counts, list_of_proteins
    
    def update_client_list(self):
        # move the coodinator cohort name to the 0 place in the list of self.clients
        cohort_name = self.id
        client_list = [cohort_name] + [client for client in self.clients if client != cohort_name]
        self.log(f'List of clients: {client_list}')
        return client_list
    
    def aggregate_plex_covariates(self, list_of_plexes):
        self.log("Start aggregation of plex covariates...")
        plex_covariates_list = []
        for plexes in list_of_plexes:
            plex_covariates_list.extend(plexes)
        plex_covariates_list = sorted(list(set(plex_covariates_list)))
        self.log(f"Size of plex covariates list: {len(plex_covariates_list)}")
        return plex_covariates_list
    
    def aggregate_counts(self, list_of_counts, prot_names):
        self.log("Start aggregation of counts...")
        min_counts = [pd.DataFrame.from_dict(local_counts, orient='index') for local_counts in list_of_counts]
        global_min_counts = pd.concat(min_counts, axis=1).min(axis=1)
        global_min_counts = global_min_counts.loc[prot_names]
        # if min in global_min_counts is 0, add 1 to all counts
        if (global_min_counts == 0).any():
            global_min_counts = global_min_counts + 1
        self.log(f"Size of global min counts: {global_min_counts.shape}")
        return global_min_counts
    
    def aggregate_protein_groups(self, list_of_proteins):
        self.log("Start aggregation of protein groups...")
        prot_names = utils.get_analyzed_proteins(list_of_proteins, self.load('only_shared_proteins'))
        self.log(f"Protein group list for analysis was created. Number of proteins: {len(prot_names)}")
        return prot_names
    
    def load_and_store_variables(self):
        self.log("Loading target classes and covariates...")
        target_classes = self.load("target_classes")
        covariates = self.load("covariates")
        variables = target_classes + covariates
        return variables

############################################################################################################
# Validation and filtering states
############################################################################################################
@app_state(name='validation')
class ValidationState(AppState):
    def register(self):
        self.register_transition('prot_na_counting', role=Role.BOTH)

    def run(self):
        self.log("Start validation...")
        # get list of common proteins from coordinator
        if self.load('experiment_type') == "TMT" and self.load('plex_covariate'):
            client_list, stored_features, global_min_counts_in, \
                global_min_counts_val, variables, plex_covariates = \
                    self.await_data(n=1, is_json=False,
                                    memo="ProteinsCounts")
        else:
            client_list, stored_features, global_min_counts_in, \
                global_min_counts_val, variables = self.await_data(n=1, is_json=False, 
                                                                      memo="ProteinsCounts")
        # load client data
        client = self.load('client')

        # update counts if use_counts is True
        if self.load('use_counts'):
            global_min_counts = pd.DataFrame({"count" : global_min_counts_val}, index = global_min_counts_in)
            client.counts = global_min_counts.copy()
            self.log("Counts are updated...")

        self.log(f"Data validation starts...")
        client.validate_inputs(stored_features, variables)

        # add cohort effect columns to each design matrix
        # if plex_covariate exists, use this column as a cohort effect
        if self.load('experiment_type') == "TMT" and self.load('plex_covariate'):
            client.add_cohort_effects_to_design(plex_covariates, self.load('plex_covariate'))
        else:
            client.add_cohort_effects_to_design(client_list)

        self.log(f"Samples in {client.cohort_name} data: {len(client.sample_names)}")        
        self.log(f"Protein groups in {client.cohort_name} data:  {len(client.prot_names)}")
        self.log(f"Design {client.cohort_name} has beed updated with cohort effects")

        self.log("Transition to filtering based on NA values...")
        return 'prot_na_counting'
    

@app_state(name='prot_na_counting')
class NACountState(AppState):
    def register(self):
        self.register_transition('prot_na_filtering', role=Role.COORDINATOR)
        self.register_transition('do_normalization', role=Role.PARTICIPANT)
    
    def run(self):
        self.log("Start NA counting...")
        client = self.load('client')
        na_count_in_variable, samples_per_class = client.apply_filters(
            min_f=self.load('max_na_rate'), 
            remove_single_peptide_prots=self.load('remove_single_pep_protein')
        )
        self.send_data_to_coordinator(
            [na_count_in_variable.to_dict(orient='index'), samples_per_class], 
            send_to_self=True, 
            use_smpc=self.load('use_smpc')
        )

        if self.is_coordinator:
            self.log(f"Transition to creation of list of proteins passed NA filter ({round(self.load('max_na_rate') * 100, 0)}%)...")
            return 'prot_na_filtering'
        
        self.log("Transition to normalization state...")
        return 'do_normalization'
    

@app_state(name='prot_na_filtering')
class NAFilterState(AppState):
    def register(self):
        self.register_transition('do_normalization', role=Role.COORDINATOR)

    def run(self):
        self.log("Start NA filtering...")
        self.log(f"Number of proteins before NA filtering: {len(self.load('protein_groups'))}")
        
        list_of_na_counts_tuples = self.gather_data(is_json=False, 
                                                    use_smpc=self.load('use_smpc'))
                                             
        # sort keep_proteins
        keep_proteins = utils.filter_features_na_rate(list_of_na_counts_tuples, self.load('max_na_rate'))

        self.log(f"Number of proteins after NA filtering: {len(keep_proteins)}")
        self.store(key='stored_features', value=keep_proteins)
        
        self.broadcast_data(keep_proteins, send_to_self=True, memo="KeepProteins")
        self.log("Transition to normalization state...")
        return 'do_normalization'

############################################################################################################
# Normalization states
############################################################################################################
@app_state(name='do_normalization')
class NormalizationState(AppState):
    def register(self):
        self.register_transition('median_calculation', role=Role.BOTH)
        self.register_transition('irs_normalization', role=Role.BOTH)
        self.register_transition('mask_preparation', role=Role.BOTH)

    def run(self):
        keep_proteins = self.await_data(n=1, is_json=False, memo="KeepProteins")
        client = self.load('client')
        client.update_prot_names(keep_proteins)
        # save keep_proteins
        self.store(key='stored_features', value=keep_proteins)

        if not self.load("use_median") and not self.load("use_irs"):
            self.log("No normalization is needed...")
            self.log("Transition to computation state...")
            return 'mask_preparation'  

        self.log("Start normalization...")        
        if self.load("use_median"):
            self.log("Transition to median calculation...")
            return "median_calculation"
        if self.load("use_irs"):
            self.log("Transition to IRS normalization...")
            return "irs_normalization"


@app_state(name='median_calculation')
class MedianCalculationState(AppState):
    def register(self):
        self.register_transition('global_median_calculation', role=Role.COORDINATOR)
        self.register_transition('median_centering', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start median normalization...")
        client = self.load('client')
        self.log("Computing local medians...")
        avg_medians = client.compute_medians()
        number_of_samples = len(client.sample_names)

        self.log("Medians are computed, sending to coordinator...")
        self.send_data_to_coordinator((avg_medians, number_of_samples),
                                      send_to_self=True, use_smpc=False)

        if self.is_coordinator:
            self.log("Transition to global medians computation...")
            return 'global_median_calculation'
        self.log("Transition to median centering...")
        return 'median_centering'


@app_state(name='global_median_calculation')
class GlobalMedianNormalizationState(AppState):
    def register(self):
        self.register_transition('median_centering', role=Role.COORDINATOR)

    def run(self):
        self.log("Start global median normalization...")
        mead_samples_tuples = self.gather_data(use_smpc=False)
        global_median = utils.aggregate_medians(mead_samples_tuples)
        self.log("Global median is computed, sending to participants...")

        self.broadcast_data(global_median, send_to_self=True, memo="GlobalMedian")
        self.log("Transition to median centering...")
        return 'median_centering'


@app_state(name='median_centering')
class MedianCenteringState(AppState):
    def register(self):
        self.register_transition('mask_preparation', role=Role.BOTH)
        self.register_transition('irs_normalization', role=Role.BOTH)

    def run(self):
        self.log("Start median centering...")
        global_median_mean = self.await_data(n=1, memo="GlobalMedian")
        client = self.load('client')
        client.mean_median_centering(global_median_mean)
        self.log("Median normalization is done...")

        if self.load("use_irs"):
            self.log("Transition to IRS normalization...")
            return 'irs_normalization'
        
        self.log("Transition to computation state...")
        return 'mask_preparation'


@app_state(name='irs_normalization')
class IRSNormalizationState(AppState):
    def register(self):
        self.register_transition('mask_preparation', role=Role.BOTH)

    def run(self):
        self.log("Start IRS normalization...")
        client = self.load('client')
        client.irsNorm_in_silico_single_center()
        self.log("IRS normalization is done...")

        self.log("Transition to computation state...")
        return 'mask_preparation'


############################################################################################################
# Mask prepation states
############################################################################################################
@app_state(name='mask_preparation')
class MaskPreparationState(AppState):
    def register(self):
        self.register_transition('mask_update', role=Role.PARTICIPANT)
        self.register_transition('mask_aggregation', role=Role.COORDINATOR)

    def run(self):
        client = self.load('client')

        # prepare for limma
        self.log("Preparing for limma...")
        client.prepare_for_limma(self.load('stored_features'))        
        self.log(f'Design colnames: {client.design.columns.values}')

        self.log("Start mask preparation...")
        mask = client.get_mask()
        self.log(f"Mask is prepared. Size: {mask.shape}")

        self.send_data_to_coordinator(mask, send_to_self=True, use_smpc=self.load('use_smpc'))
        
        if self.is_coordinator:
            self.store(key='variables', value=client.design.columns.values)
            self.log("Transition to mask aggregation...")
            return 'mask_aggregation'
        self.log("Transition to mask update...")
        return 'mask_update'


@app_state(name='mask_aggregation')
class MaskAggregationState(AppState):
    def register(self):
        self.register_transition('mask_update', role=Role.COORDINATOR)

    def run(self):
        self.log("Start aggregating mask...")
        mask = self.gather_data(use_smpc=self.load('use_smpc'), is_json=False)

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))        
        client_number = len(self.load('client_list'))

        global_mask = utils.aggregate_masks(mask, n, k, used_SMPC=self.load('use_smpc'),
                                          client_number=client_number)
        
        self.log(f"Mask is aggregated, size: {global_mask.shape}")
        
        self.broadcast_data(global_mask, send_to_self=True, memo="Mask")
        self.log("Transition to mask update...")
        return 'mask_update'


@app_state(name='mask_update')
class MaskUpdateState(AppState):
    def register(self):
        self.register_transition('compute_XtY_XTX', role=Role.PARTICIPANT)
        self.register_transition('save_mask', role=Role.COORDINATOR)

    def run(self):
        self.log("Start updating mask...")
        mask = self.await_data(n=1, is_json=False, memo="Mask")
        mask = np.array(mask)
        self.log(f"Mask is received..., size: {mask.shape}, type: {type(mask)}")

        client = self.load('client')
        updated_masks = client.updated_mask(mask)
        
        self.log(f"Mask is updated, size: {updated_masks.shape}")
        self.send_data_to_coordinator(updated_masks, send_to_self=True, use_smpc=self.load('use_smpc'))

        if self.is_coordinator:
            self.log("Transition to saving mask...")
            return 'save_mask'

        self.log("Transition to computation state...")
        return 'compute_XtY_XTX'


@app_state(name='save_mask')
class SaveMaskState(AppState):
    def register(self):
        self.register_transition('compute_XtY_XTX', role=Role.COORDINATOR)

    def run(self):
        self.log("Start saving mask...")
        masks = self.gather_data(use_smpc=self.load('use_smpc'), is_json=False)
        
        k = len(self.load('variables'))
        n = len(self.load('stored_features'))
        mask_glob = utils.aggregate_masks(masks, n, k, second_round=True, used_SMPC=self.load('use_smpc'))
        self.log("Final global mask is saved...")

        self.log("Broadcasting mask...")
        self.broadcast_data(mask_glob, send_to_self=True, memo="MaskGlobal")
        self.log("Transition to computation state...")
        return 'compute_XtY_XTX'


############################################################################################################
# Linear model states
############################################################################################################
@app_state(name='compute_XtY_XTX')
class ComputeXtState(AppState):
    def register(self):
        self.register_transition('compute_SSE', role=Role.PARTICIPANT,
                                 label='Compute SSE and cov_coef')
        self.register_transition('compute_beta', role=Role.COORDINATOR,
                                 label='Compute beta and beta stdev')

    def run(self):
        client = self.load('client')

        global_mask = self.await_data(n=1, memo="MaskGlobal")
        global_mask = np.array(global_mask)
        self.store(key='MaskGlobal', value=global_mask)

        self.log("Start computation of XtY and XTX...")
        XtX, XtY = client.compute_XtX_XtY()

        self.log(f"k for XTX computation: {client.design.shape[1]}")
        self.log("XtX and XtY are computed, sending to coordinator...")

        self.send_data_to_coordinator([XtX, XtY], send_to_self=True, use_smpc=self.load('use_smpc'))

        if self.is_coordinator:
            self.log("Transition to computing beta...")
            return 'compute_beta'
        self.log("Transition to computing SSE...")
        return 'compute_SSE'


@app_state(name='compute_beta')
class ComputeBetaState(AppState):
    def register(self):
        self.register_transition('compute_SSE', role=Role.COORDINATOR,
                                 label='Compute SSE and cov_coef')
        
    def run(self):
        self.log("Start computation of beta and beta stdev...")
        list_of_xt_lists = self.gather_data(use_smpc=self.load('use_smpc'))

        self.log("XtX and XtY are received, computing beta and beta stdev...")

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))
        self.log(f"Size for computing global beta and beta stdev, k = {k}, n = {n}...")
    
        XtX_glob, XtY_glob = utils.aggregate_XtX_XtY(list_of_xt_lists, n, k, self.load('use_smpc'))
        self.log("Computing beta and beta stdev...")        
    
        beta, stdev_unscaled = utils.compute_beta_and_stdev(np.array(XtX_glob), np.array(XtY_glob), 
                                                            n, k, 
                                                            self.load('MaskGlobal'))
        self.store(key='beta', value=beta)
        self.store(key='stdev_unscaled', value=stdev_unscaled)

        self.log("Beta and beta stdev are computed, sending to participants...")
        self.broadcast_data(beta, 
                            send_to_self=True, 
                            memo="BetaStdev")
        return 'compute_SSE'
    

@app_state(name='compute_SSE')
class ComputeSSEState(AppState):
    def register(self):
        self.register_transition('aggregate_SSE', role=Role.COORDINATOR)
        self.register_transition('write_results', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start computation of SSE and cov_coef...")
        beta = self.await_data(n=1, memo="BetaStdev")
        beta = np.array(beta)

        client = self.load('client')
        SSE, cov_coef = client.compute_SSE_and_cov_coef(beta, self.load('MaskGlobal'))
        intensities_sum = np.array(client.sum_intensities())
        n_measurements = np.array(client.get_not_na())

        self.log("SSE and cov_coef are computed, sending to coordinator...")
        self.send_data_to_coordinator(
            [SSE, cov_coef, intensities_sum, n_measurements],
            send_to_self=True,
            use_smpc=self.load('use_smpc')
        )
        
        if self.is_coordinator:
            self.log("Transition to aggregation of SSE...")
            return 'aggregate_SSE'
        self.log("Transition to writing results...")
        return 'write_results'


@app_state(name='aggregate_SSE')
class AggregateSSEState(AppState):
    def register(self):
        self.register_transition('make_contrasts', role=Role.COORDINATOR)

    def run(self):
        self.log("Start aggregation of SSE...")
        list_of_sse_cov_coef = self.gather_data(use_smpc=self.load('use_smpc'))

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))

        SSE, cov_coef, n_measurements, Amean = utils.aggregate_SSE_and_cov_coef(
            list_of_sse_cov_coef, n, k, self.load('use_smpc'), len(self.load('client_list'))
        )
        
        self.log("Aggregation of SSE is done, start computing global parameters...")
        sigma, cov_coef, df_residual, Amean, var = utils.compute_SSE_and_cov_coef_global(
            cov_coef, SSE, Amean, n_measurements, n, self.load('MaskGlobal')
        )

        self.store(key='cov_coef', value=cov_coef)
        self.store(key='var', value=var)
        self.store(key='sigma', value=sigma)
        self.store(key='df_residual', value=df_residual)
        self.store(key='Amean', value=Amean)

        self.log("Global SSE and cov_coef are computed...")
        self.log("Transition to making contrasts...")

        return 'make_contrasts'


############################################################################################################
# Apply contrasts states
############################################################################################################
@app_state(name='make_contrasts')
class MakeContrastsState(AppState):
    """
    Creates contrast matrix given deisgn matrix and pairs or columns to compare.
    For example:
    contrasts = [([A],[B]),([A,B],[C,D])] defines two contrasts:
    A-B and (A and B) - (C and D).
    """
    def register(self):
        self.register_transition('fit_contasts', role=Role.COORDINATOR)

    def run(self):
        self.log("Start making contrasts...")
        target_classes = self.load('target_classes')        
        contrasts_df = utils.make_contrasts(target_classes, self.load('variables'))
        self.store(key='contrasts', value=contrasts_df)

        self.log("Contrasts are computed...")
        self.log("Transition to fitting contrasts...")
        return 'fit_contasts'


@app_state(name='fit_contasts')
class FitContrastsState(AppState):
    def register(self):
        self.register_transition('ebayes', role=Role.COORDINATOR)

    def run(self):
        self.log("Start fitting contrasts...")
        contrast_matrix = self.load('contrasts').values
        ncoef = self.load('cov_coef').shape[1]
        
        beta, stdev_unscaled, cov_coef = utils.fit_contrasts(
            self.load('beta'), 
            contrast_matrix,
            self.load('cov_coef'), 
            ncoef,
            self.load('stdev_unscaled')
        )
        self.store(key='beta', value=beta)
        self.store(key='cov_coef', value=cov_coef)
        self.store(key='stdev_unscaled', value=stdev_unscaled)
        
        self.log("Contrasts are fitted...")
        self.log("Transition to eBayes stage...")
        return 'ebayes'


############################################################################################################
# eBayes state
############################################################################################################
@app_state(name='ebayes')
class eBayesState(AppState):
    def register(self):
        self.register_transition('spectral_count_ebayes', role=Role.COORDINATOR)
        self.register_transition('write_results', role=Role.COORDINATOR)

    def run(self):
        self.log("Start eBayes stage...")

        self.log("Calculating moderated t-statistics...")
        results, df_total = utils.moderatedT(
            self.load('var'), 
            self.load('df_residual'),
            self.load('beta'), 
            self.load('stdev_unscaled')
        )
        results["AveExpr"] = self.load('Amean')
        self.log("Result table is pre-computed...")

        self.log("Computing B statistic...")
        results = utils.Bstat(
            df_total, 
            self.load('stdev_unscaled'), 
            results,
            stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)
        self.log("B statistic is computed...")

        self.log("Computing p-values...")        
        results = utils.topTableT(results,
                                  self.load('stored_features'), 
                                  self.load('beta'), 
                                  self.load('stdev_unscaled'), 
                                  df_total,
                                  adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
        self.log("P-values are computed...")
        self.store(key='results', value=results)
        
        if self.load('use_counts'):
            self.log("Transition to count adjustment...")
            return 'spectral_count_ebayes'
        else:
            self.log("Broadcasting results...")
            self.broadcast_data(results.to_dict(orient = 'index'), send_to_self=True, memo="Results")
            self.log("Transition to writing results...")
            return 'write_results'


############################################################################################################
# Spectral count eBayes state
############################################################################################################
@app_state(name='spectral_count_ebayes')
class SpectralCounteBayesState(AppState):
    def register(self):
        self.register_transition('write_results', role=Role.COORDINATOR)

    def run(self):
        self.log("Start spectral count eBayes stage...")

        global_min_counts = self.load('min_counts') + 1
        global_min_counts = global_min_counts.loc[self.load('stored_features')]
        self.log(f"Size of global min counts: {global_min_counts.shape}")

        results = self.load('results')

        results = utils.spectral_count_ebayes(
            results, 
            global_min_counts, 
            self.load('stored_features'),
            self.load('beta'), 
            self.load('stdev_unscaled'),
            self.load('df_residual'), 
            self.load('sigma'),
            return_sorted=False, fit_method="loess")
        
        self.store(key='results', value=results)
        self.log("Spectral count eBayes stage is done...")
        
        self.log("Broadcasting results...")
        self.broadcast_data(results.to_dict(orient = 'index'), send_to_self=True, memo="Results")

        self.log("Transition to writing results...")
        return 'write_results'


############################################################################################################
# Writing results state
############################################################################################################
@app_state(name='write_results')
class WriteResults(AppState):
    def register(self):
        self.register_transition('terminal', label='Done!')

    def run(self):
        self.log("Awaiting results from coordinator...")
        results = self.await_data(n=1, memo="Results")
        results = pd.DataFrame.from_dict(results, orient='index')

        self.log("Writing results...")
        output_file = os.path.join(os.getcwd(), OUTPUT_DIR, self.load("result_table"))
        results.to_csv(output_file, sep="\t")
        self.log(f"Results written to {output_file}")
        self.log("Done!")
        return 'terminal'
