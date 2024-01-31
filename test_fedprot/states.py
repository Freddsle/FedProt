from __future__ import annotations
from FeatureCloud.app.engine.app import AppState, app_state, Role
import time
import os
import bios

import pandas as pd
import numpy as np
from scipy import linalg

from utils import Client

CONFIG_FILE = "config.yml"
APP_DIR = "/app"
MOUNT_DIR = "/mnt/input"
EXPERIMENT_TYPE = 'DIA'
LOG_TRANSFORMED = False
REMOVE_SINGLE_PEPTIDE_PROT = False
USE_SMPC = True


@app_state(name='initial', role=Role.BOTH)
class InitialState(AppState):
    def register(self):
        self.register_transition('common_proteins', Role.COORDINATOR)
        self.register_transition('validation', Role.PARTICIPANT)

    def run(self):
        # read config
        self.read_config()

        if USE_SMPC:
            self.configure_smpc()
        # defining the client
        cohort_name = self.id
        
        # # initializing the client includes loading and preprocessing of the data
        client = Client(
            cohort_name,
            intensities_file_path = os.path.join(MOUNT_DIR, self.load('intensities_name')),
            count_file_path = os.path.join(MOUNT_DIR, self.load('counts_name')),
            annotation_file_path = os.path.join(MOUNT_DIR, self.load('design_name')),
            experiment_type = EXPERIMENT_TYPE,
            log_transformed = LOG_TRANSFORMED
        ) 
        self.log(f"Data read! {client.intensities.shape[1]} samples, {client.intensities.shape[0]} proteins")

        # store client
        self.store(key='client', value=client)
        # self.configure_smpc()

        # send list of protein names (genes) to coordinator
        self.log("[initial] Sending the protein names list to the coordinator")
        self.send_data_to_coordinator(client.prot_names, send_to_self=True, use_smpc=False)

        # initialize app
        if self.is_coordinator:
            self.log("Transition to creation common proteins list...")
            return 'common_proteins'
        self.log("Transition to validation...")
        return 'validation'
        
    def read_config(self):
        self.log("Reading config...")
        config = bios.read(os.path.join(APP_DIR, CONFIG_FILE))
        config = config["fedprot"]
            
        self.store('intensities_name', config['intensities'])
        self.store('counts_name', config['counts'])
        self.store('design_name', config['design'])
        self.store('sep', config['sep'])
        self.store('target_classes', config['target_classes'])
        self.store('covariates', config['covariates'])
        self.store('max_na_rate', config['max_na_rate'])
        self.store('result_table', config['result_table'])
        


@app_state(name='common_proteins')
class CommonProteinsState(AppState):
    def register(self):
        self.register_transition('validation', role=Role.COORDINATOR,
                                 label='Broadcasting common proteins list')

    def run(self):
        self.log("[common_proteins] Gathering proteins lists from all clients")
        list_of_features_lists = self.gather_data(is_json=False)
        self.log("[common_proteins] Gathering done!")
        
        prot_names = list()
        for features_list in list_of_features_lists:
            if len(prot_names) == 0:
                prot_names =  sorted(set(features_list))
            else:
                prot_names = sorted(list(set(prot_names) & set(features_list)))

        print(f"[common_proteins] Common proteins list was created. Number of proteins: {len(prot_names)}")
        self.store('common_proteins', prot_names)

        # load target classes and broadcast them to all clients
        target_classes = self.load("target_classes")
        covariates = self.load("covariates")
        variables = target_classes + covariates

        self.store('variables', variables)

        self.log("Transition to validation ...")
        self.broadcast_data((prot_names, variables), send_to_self=True, memo="commonProteins")
        return 'validation'


@app_state(name='validation')
class ValidationState(AppState):
    def register(self):
        self.register_transition('prot_na_counting', role=Role.BOTH)

    def run(self):
        self.log("Start validation...")
        # get list of common proteins from coordinator
        stored_features, variables = self.await_data(n=1, is_json=False, memo="commonProteins")
        
        # load client data
        client = self.load('client')
        client.validate_inputs(stored_features, variables)
        # add cohort effects to design
        client.add_cohort_effects_to_design(self.clients[1:])

        self.log(f"Samples in {client.cohort_name} data: {len(client.sample_names)}")        
        self.log(f"Protein groups in {client.cohort_name} data:  {len(client.prot_names)}")
        self.log(f"Design {client.cohort_name} has beed updated with cohort effects")

        self.log("Transition to filtering based on NA values...")
        return 'prot_na_counting'
    

@app_state(name='prot_na_counting')
class NACountState(AppState):
    def register(self):
        self.register_transition('prot_na_filtering', role=Role.COORDINATOR)
        self.register_transition('compute_XtY_XTX', role=Role.PARTICIPANT)
    
    def run(self):
        self.log("Start NA counting...")
        client = self.load('client')
        na_count_in_variable, samples_per_class = client.apply_filters(min_f=self.load('max_na_rate'), 
                                                                       remove_single_peptide_prots=REMOVE_SINGLE_PEPTIDE_PROT)
        self.send_data_to_coordinator([na_count_in_variable.to_dict(orient='index'), samples_per_class], 
                                      send_to_self=True, 
                                      use_smpc=USE_SMPC)        
        if self.is_coordinator:
            self.log(f"Transition to creation of list of proteins passed NA filter ({round(self.load('max_na_rate') * 100, 0)}%)...")
            return 'prot_na_filtering'
        self.log("Transition to computation state...")
        return 'compute_XtY_XTX'
    

@app_state(name='prot_na_filtering')
class NAFilterState(AppState):
    def register(self):
        self.register_transition('compute_XtY_XTX', role=Role.COORDINATOR)

    def run(self):
        self.log("Start NA filtering...")
        self.log(f"Number of proteins before NA filtering: {len(self.load('common_proteins'))}")
        list_of_na_counts_tuples = self.gather_data(is_json=False, use_smpc=USE_SMPC)
        

        # merge na counts
        samples_per_target = dict()
        prot_na_table = pd.DataFrame()

        for na_pair in list_of_na_counts_tuples:
            na_count_client = pd.DataFrame.from_dict(na_pair[0], orient='index')
            samples_per_class_client = na_pair[1]

            for key in samples_per_class_client:
                samples_per_target[key] = samples_per_class_client[key] + samples_per_target.get(key, 0)
            if prot_na_table.empty:
                prot_na_table = na_count_client
            else:
                prot_na_table = na_count_client.add(prot_na_table)
    
        # filter proteins based on NA counts
        prot_na_table = prot_na_table.loc[:, samples_per_target.keys()]
        samples_series = pd.Series(samples_per_target)
        na_perc = prot_na_table.div(samples_series, axis=1)

        keep_proteins = na_perc[na_perc.apply(lambda row: all(row < self.load('max_na_rate')), axis=1)].index.values

        self.log(f"Number of proteins after NA filtering: {len(keep_proteins)}")
        self.store('stored_features', keep_proteins)
        self.broadcast_data(keep_proteins, send_to_self=True, memo="KeepProteins")
        self.log("Transition to computation state...")
        return 'compute_XtY_XTX'  


@app_state(name='compute_XtY_XTX')
class ComputeXtState(AppState):
    def register(self):
        self.register_transition('compute_SSE', role=Role.PARTICIPANT,
                                 label='Compute SSE and cov_coef')
        self.register_transition('compute_beta', role=Role.COORDINATOR,
                                 label='Compute beta and beta stdev')

    def run(self):
        keep_proteins = self.await_data(n=1, is_json=False, memo="KeepProteins")
        client = self.load('client')
        client.update_prot_names(keep_proteins)

        self.log("Start computation of XtY and XTX...")
        client.prepare_for_limma(keep_proteins)
        XtX, XtY = client.compute_XtX_XtY()

        self.log(f"k for XTX computation: {client.design.shape[1]}")

        self.log("XtX and XtY are computed, sending to coordinator...")

        if USE_SMPC:
            XtX_3d_reshaped = {[i, j, k]: XtX[i, j, k] for i in range(XtX.shape[0]) for j in range(XtX.shape[1]) for k in range(XtX.shape[2])}
            XtY_2d_reshaper = {[i, j]: XtY[i, j] for i in range(XtY.shape[0]) for j in range(XtY.shape[1])}
            self.send_data_to_coordinator([XtX_3d_reshaped, XtY_2d_reshaper],
                                          send_to_self=True,
                                          use_smpc=USE_SMPC,
                                          )
        else:
            self.send_data_to_coordinator([XtX, XtY],
                                        send_to_self=True,
                                        use_smpc=USE_SMPC,
                                        )
        
        if self.is_coordinator:
            self.store('variables', client.design.columns.values)
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
        list_of_xt_lists = self.gather_data(is_json=False, use_smpc=USE_SMPC)

        self.log("XtX and XtY are received, computing beta and beta stdev...")

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))
        XtX_glob = np.zeros((n, k, k))
        XtY_glob = np.zeros((n, k))
        stdev_unscaled = np.zeros((n, k))
        beta = np.zeros((n, k))

        self.log(f"Size for computing global beta and beta stdev, k = {k}, n = {n}...")
    
        # non-smpc case, need to aggregate
        XtX_list = list()
        XtY_list = list()
        self.log(f"Size of list of XtX and XtY lists: {len(list_of_xt_lists)}")
        
        for pair in list_of_xt_lists:
            if USE_SMPC:
                shape_3d = XtX_glob.shape
                shape_2d = XtY_glob.shape
                
                XtX_local = np.array([pair[0][(i, j, k)] for i in range(shape_3d[0]) for j in range(shape_3d[1]) for k in range(shape_3d[2])]).reshape(shape_3d)
                XtY_local = np.array([pair[1][(i, j)] for i in range(shape_2d[0]) for j in range(shape_2d[1])]).reshape(shape_2d)
            else:
                XtX_local = pair[0]
                XtY_local = pair[1]
            XtX_list.append(XtX_local)
            XtY_list.append(XtY_local)

        for i in range(0, len(self.clients)):
            self.log(f"Size of XtX and XtY for i={i} : {XtX_list[i].shape}, {XtY_list[i].shape}")
            XtX_glob += XtX_list[i]
            XtY_glob += XtY_list[i]
            self.log(f"Size of XtX_glob and XtY_glob after i={i} : {XtX_glob.shape}, {XtY_glob.shape}")

        self.log("Computing beta and beta stdev...")
        for i in range(0, n):
            if linalg.det(XtX_glob[i, :, :]) == 0:
                self.log(f"XtX is singular for protein {i}, determinant is 0.")
            invXtX = linalg.inv(XtX_glob[i, :, :])
            beta[i, :] = invXtX @ XtY_glob[i, :]
            stdev_unscaled[i, :] = np.sqrt(np.diag(invXtX))

        self.store('stdev_unscaled', stdev_unscaled)
        self.log("Beta and beta stdev are computed, sending to participants...")
        self.broadcast_data(beta, 
                            send_to_self=True, 
                            memo="BetaStdev")
        return 'compute_SSE'
    

@app_state(name='compute_SSE')
class ComputeSSEState(AppState):
    def register(self):
        self.register_transition('aggregate_SSE', role=Role.COORDINATOR)
        self.register_transition('get_counts', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start computation of SSE and cov_coef...")
        beta = self.await_data(n=1, is_json=False, memo="BetaStdev")

        client = self.load('client')
        SSE, cov_coef = client.compute_SSE_and_cov_coef(beta)
        intensities_sum = client.sum_intensities()
        n_measurements = client.get_not_na()

        self.log("SSE and cov_coef are computed, sending to coordinator...")
        self.send_data_to_coordinator([SSE, cov_coef, intensities_sum, n_measurements],
                                        send_to_self=True,
                                        use_smpc=USE_SMPC)
        
        if self.is_coordinator:
            self.log("Transition to aggregation of SSE...")
            return 'aggregate_SSE'
        self.log("Transition to getting counts...")
        return 'get_counts'
    

@app_state(name='aggregate_SSE')
class AggregateSSEState(AppState):
    def register(self):
        self.register_transition('make_contrasts', role=Role.COORDINATOR)

    def run(self):
        self.log("Start aggregation of SSE...")
        list_of_sse_cov_coef = self.gather_data(is_json=False, use_smpc=USE_SMPC)

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))

        SSE = np.zeros(n)
        cov_coef = np.zeros((k, k))
        n_measurements = np.zeros(n)
        Amean = np.zeros(n)

        # if not USE_SMPC:
        #     # non-smpc case, need to aggregate
        #     XtX_list = list()
        #     XtY_list = list()

        #     for pair in list_of_xt_lists:
        #         XtX_list.append(pair[0])
        #         XtY_list.append(pair[1])

        #     for i in range(0, len(self.clients)):
        #         XtX_glob += XtX_list[i]
        #         XtY_glob += XtY_list[i]

        # else:
        #     # smpc case, already aggregated
        #     list_of_sse_cov_coef = list_of_sse_cov_coef[0]
        #     SSE = list_of_sse_cov_coef[0]
        #     cov_coef = list_of_sse_cov_coef[1]
        #     intensities_sum = list_of_sse_cov_coef[2]
        #     n_measurements = list_of_sse_cov_coef[3]

        return 'make_contrasts'



@app_state(name='make_contrasts')
class MakeContrastsState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.COORDINATOR)

    def run(self):
        self.log("Start making contrasts...")
        
        return 'computing'


@app_state(name='get_counts')
class GetCountsState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start getting counts...")
        return 'computing'



@app_state(name='computing')
class ComputingState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.PARTICIPANT,
                                 label='Wait for another round of local training')
        self.register_transition('global_aggregation', role=Role.COORDINATOR,
                                 label='Collect and aggregate local models')
        self.register_transition('write_results', role=Role.BOTH, label='Broadcasting results')

    def run(self):
        self.log("Start computation...")
        time.sleep(3)
        self.log("Transition to aggregation or writing...")

        self.log(f"Aggregated: {self.load('aggregated')}")
        self.log(f"Coordinator: {self.is_coordinator}")
        
        if self.load('aggregated') or not self.is_coordinator:
            return 'write_results'
        else:
            return 'global_aggregation'


@app_state(name='computing')
class ComputingState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.PARTICIPANT,
                                 label='Wait for another round of local training')
        self.register_transition('global_aggregation', role=Role.COORDINATOR,
                                 label='Collect and aggregate local models')
        self.register_transition('write_results', role=Role.BOTH, label='Broadcasting results')

    def run(self):
        self.log("Start computation...")
        time.sleep(3)
        self.log("Transition to aggregation or writing...")

        self.log(f"Aggregated: {self.load('aggregated')}")
        self.log(f"Coordinator: {self.is_coordinator}")
        
        if self.load('aggregated') or not self.is_coordinator:
            return 'write_results'
        else:
            return 'global_aggregation'


@app_state(name='computing')
class ComputingState(AppState):
    def register(self):
        self.register_transition('computing', role=Role.PARTICIPANT,
                                 label='Wait for another round of local training')
        self.register_transition('global_aggregation', role=Role.COORDINATOR,
                                 label='Collect and aggregate local models')
        self.register_transition('write_results', role=Role.BOTH, label='Broadcasting results')

    def run(self):
        self.log("Start computation...")
        time.sleep(3)
        self.log("Transition to aggregation or writing...")

        self.log(f"Aggregated: {self.load('aggregated')}")
        self.log(f"Coordinator: {self.is_coordinator}")
        
        if self.load('aggregated') or not self.is_coordinator:
            return 'write_results'
        else:
            return 'global_aggregation'
            

@app_state(name='global_aggregation')
class GlobalAggregation(AppState):
    def register(self):
        self.register_transition('computing', label='Broadcas final results')

    def run(self):
        self.log("Aggregating intermediate results...")
        self.store('aggregated', True)
        return 'computing'


@app_state(name='write_results')
class WriteResults(AppState):
    def register(self):
        self.register_transition('terminal', label='Done!')

    def run(self):
        self.log("Writing results...")
        self.log("Done!")
        return 'terminal'
