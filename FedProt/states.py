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
        self.register_transition('common_proteins', Role.COORDINATOR)
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
            log_transformed = self.load('log_transformed')
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
        config = bios.read(os.path.join(MOUNT_DIR, CONFIG_FILE))
        config = config["fedprot"]
            
        self.store(key='experiment_type', value=config['experiment_type'])
        self.store(key='remove_single_pep_protein', value=config['remove_single_pep_protein'])
        
        self.store(key='intensities_name', value=config['intensities'])
        self.store(key='use_counts', value=config['use_counts'])
        if config['use_counts']:
            self.store(key='counts_name', value=config['counts'])
        self.store(key='design_name', value=config['design'])
        self.store(key='sep', value=config['sep'])

        self.store(key='target_classes', value=config['target_classes'])
        self.store(key='covariates', value=config['covariates'])

        self.store(key='max_na_rate', value=config['max_na_rate'])
        self.store(key='use_smpc', value=config['use_smpc'])
        self.store(key='log_transformed', value=config['log_transformed'])
        
        self.store(key='result_table', value=config['result_table'])
        

@app_state(name='common_proteins')
class CommonProteinsState(AppState):
    def register(self):
        self.register_transition('validation', role=Role.COORDINATOR,
                                 label='Broadcasting common proteins list')

    def run(self):
        self.log("[common_proteins] Gathering proteins lists from all clients")
        list_of_features_lists = self.gather_data(is_json=False)
        self.log("[common_proteins] Gathering done!")

        # move the coodinator cohort name to the 0 place in the list of self.clients
        cohort_name = self.id
        client_list = [cohort_name] + [client for client in self.clients if client != cohort_name]
        self.store(key='client_list', value=client_list)
        self.log(f'List of clients: {client_list}')

        prot_names = utils.get_common_proteins(list_of_features_lists)

        print(f"[common_proteins] Common proteins list was created. Number of proteins: {len(prot_names)}")
        self.store(key='common_proteins', value=prot_names)

        # load target classes and broadcast them to all clients
        target_classes = self.load("target_classes")
        covariates = self.load("covariates")
        variables = target_classes + covariates

        self.store(key='variables', value=variables)

        self.log("Transition to validation ...")
        self.broadcast_data((client_list, prot_names, variables), send_to_self=True, memo="commonProteins")
        return 'validation'


@app_state(name='validation')
class ValidationState(AppState):
    def register(self):
        self.register_transition('prot_na_counting', role=Role.BOTH)

    def run(self):
        self.log("Start validation...")
        # get list of common proteins from coordinator
        client_list, stored_features, variables = self.await_data(n=1, is_json=False, memo="commonProteins")
        # load client data
        client = self.load('client')
        client.validate_inputs(stored_features, variables)
        # add cohort effects to design
        client.add_cohort_effects_to_design(client_list[1:])

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
        
        self.log("Transition to computation state...")
        return 'compute_XtY_XTX'
    

@app_state(name='prot_na_filtering')
class NAFilterState(AppState):
    def register(self):
        self.register_transition('compute_XtY_XTX', role=Role.COORDINATOR)

    def run(self):
        self.log("Start NA filtering...")
        self.log(f"Number of proteins before NA filtering: {len(self.load('common_proteins'))}")
        
        list_of_na_counts_tuples = self.gather_data(is_json=False, use_smpc=self.load('use_smpc'))
        # sort keep_proteins
        keep_proteins = utils.filter_features_na_rate(list_of_na_counts_tuples, self.load('max_na_rate'))
        
        self.log(f"Number of proteins after NA filtering: {len(keep_proteins)}")
        self.store(key='stored_features', value=keep_proteins)
        
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

        self.send_data_to_coordinator(
            [XtX, XtY],
            send_to_self=True,
            use_smpc=self.load('use_smpc'),
        )
        self.log(f'Design colnames: {client.design.columns.values}')

        if self.is_coordinator:
            self.store(key='variables', value=client.design.columns.values)
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
    
        beta, stdev_unscaled = utils.compute_beta_and_stdev(XtX_glob, XtY_glob, n, k)
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
        self.register_transition('get_counts', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start computation of SSE and cov_coef...")
        beta = self.await_data(n=1, memo="BetaStdev")

        client = self.load('client')
        SSE, cov_coef = client.compute_SSE_and_cov_coef(beta)
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
        else:
            self.log("Transition to getting counts...")
            return 'get_counts'
    

@app_state(name='aggregate_SSE')
class AggregateSSEState(AppState):
    def register(self):
        self.register_transition('make_contrasts', role=Role.COORDINATOR)

    def run(self):
        self.log("Start aggregation of SSE...")
        list_of_sse_cov_coef = self.gather_data(use_smpc=self.load('use_smpc'))

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))

        SSE, cov_coef, Amean, n_measurements = utils.aggregate_SSE_and_cov_coef(
            list_of_sse_cov_coef, n, k, self.load('use_smpc'), len(self.load('client_list'))
        )
        
        self.log("Aggregation of SSE is done, start computing global parameters...")

        sigma, cov_coef, df_residual, Amean, var = utils.compute_SSE_and_cov_coef_global(
            cov_coef, SSE, Amean, n_measurements, n, k
        )

        self.store(key='cov_coef', value=cov_coef)
        self.store(key='var', value=var)
        self.store(key='sigma', value=sigma)
        self.store(key='df_residual', value=df_residual)
        self.store(key='Amean', value=Amean)

        self.log("Global SSE and cov_coef are computed...")
        self.log("Transition to making contrasts...")

        return 'make_contrasts'


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
        self.store(key='beta', value=beta)
        self.store(key='stdev_unscaled', value=stdev_unscaled)

        beta, cov_coef, stdev_unscaled = utils.fit_contrasts(
            self.load('beta'), 
            contrast_matrix,
            self.load('cov_coef'), 
            self.load('var'), 
            ncoef
        )
        self.store(key='beta', value=beta)
        self.store(key='cov_coef', value=cov_coef)
        self.store(key='stdev_unscaled', value=stdev_unscaled)
        
        self.log("Contrasts are fitted...")
        self.log("Transition to eBayes stage...")
        return 'ebayes'


@app_state(name='ebayes')
class eBayesState(AppState):
    def register(self):
        self.register_transition('get_counts', role=Role.COORDINATOR)
        self.register_transition('write_results', role=Role.BOTH)

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
                                  self.load('beta'), self.load('stdev_unscaled'), df_total,
                                  adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
        self.log("P-values are computed...")
        self.store(key='results', value=results)
        
        if self.load('use_counts'):
            self.log("Transition to getting counts...")
            return 'get_counts'
        else:
            self.log("Transition to writing results...")
            return 'write_results'


@app_state(name='get_counts')
class GetCountsState(AppState):
    def register(self):
        self.register_transition('aggregate_counts', role=Role.COORDINATOR)
        self.register_transition('write_results', role=Role.PARTICIPANT)

    def run(self):
        self.log("Start getting counts...")
        client = self.load('client')
        counts = client.get_min_count()
        self.log("Counts are computed, sending to coordinator...")
        self.send_data_to_coordinator(counts.to_dict(),  send_to_self=True, use_smpc=False)
        
        if self.is_coordinator:
            self.log("Transition to aggregation of counts...")
            return 'aggregate_counts'
        
        self.log("Transition to writing results...")
        return 'write_results'


@app_state(name='aggregate_counts')
class AggregateCountsState(AppState):
    def register(self):
        self.register_transition('spectral_count_ebayes', role=Role.COORDINATOR)

    def run(self):
        self.log("Start aggregation of counts...")
        list_of_counts = self.gather_data(is_json=False, use_smpc=False)
        min_counts = list()

        for local_counts in list_of_counts:
            min_counts.append(pd.DataFrame.from_dict(local_counts, orient='index'))
        global_min_counts = pd.concat(min_counts, axis=1)
        # else:
        #     # smpc case, already aggregated
        #     global_min_counts = pd.DataFrame.from_dict(list_of_counts[0], orient='index')

        self.log("Aggregation of counts is done...")
        global_min_counts = global_min_counts.min(axis=1).loc[self.load('stored_features')] + 1
        self.store(key='min_counts', value=global_min_counts)
        self.log(f"Size of global min counts: {global_min_counts.shape}")

        self.log("Transition to calculation of statistics based on counts...")
        return 'spectral_count_ebayes'


@app_state(name='spectral_count_ebayes')
class SpectralCounteBayesState(AppState):
    def register(self):
        self.register_transition('write_results', role=Role.COORDINATOR)

    def run(self):
        self.log("Start spectral count eBayes stage...")
        results = self.load('results')
        results = utils.spectral_count_ebayes(
            results, self.load('min_counts'), self.load('stored_features'),
            self.load('beta'), self.load('stdev_unscaled'),
            self.load('df_residual'), self.load('sigma'),
            return_sorted=False, fit_method="loess")
        
        self.store(key='results', value=results)
        self.log("Spectral count eBayes stage is done...")
        
        self.log("Broadcasting results...")
        self.broadcast_data(results.to_dict(orient = 'index'), send_to_self=True, memo="Results")

        self.log("Transition to writing results...")
        return 'write_results'


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
