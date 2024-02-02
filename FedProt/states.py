from __future__ import annotations
from FeatureCloud.app.engine.app import AppState, app_state, Role
import time
import os
import bios

import pandas as pd
import numpy as np
from scipy import linalg

from client import Client
import utils

CONFIG_FILE = "config.yml"
APP_DIR = "/app"
MOUNT_DIR = "/mnt/input"
OUTPUT_DIR = "/mnt/output"
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
            
        self.store(key='intensities_name', value=config['intensities'])
        self.store(key='counts_name', value=config['counts'])
        self.store(key='design_name', value=config['design'])
        self.store(key='sep', value=config['sep'])
        self.store(key='target_classes', value=config['target_classes'])
        self.store(key='covariates', value=config['covariates'])
        self.store(key='max_na_rate', value=config['max_na_rate'])
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
        
        prot_names = list()
        for features_list in list_of_features_lists:
            if len(prot_names) == 0:
                prot_names =  sorted(set(features_list))
            else:
                prot_names = sorted(list(set(prot_names) & set(features_list)))

        print(f"[common_proteins] Common proteins list was created. Number of proteins: {len(prot_names)}")
        self.store(key='common_proteins', value=prot_names)

        # load target classes and broadcast them to all clients
        target_classes = self.load("target_classes")
        covariates = self.load("covariates")
        variables = target_classes + covariates

        self.store(key='variables', value=variables)

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

        self.send_data_to_coordinator([XtX, XtY],
                                    send_to_self=True,
                                    use_smpc=USE_SMPC,
                                    )
        
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
        list_of_xt_lists = self.gather_data(use_smpc=USE_SMPC)

        self.log("XtX and XtY are received, computing beta and beta stdev...")

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))
        XtX_glob = np.zeros((n, k, k))
        XtY_glob = np.zeros((n, k))
        stdev_unscaled = np.zeros((n, k))
        beta = np.zeros((n, k))

        self.log(f"Size for computing global beta and beta stdev, k = {k}, n = {n}...")
    
        # non-smpc case, need to aggregate
        
        self.log(f"Size of list of XtX and XtY lists: {len(list_of_xt_lists)}")
        XtX_list = list()
        XtY_list = list()

        if not USE_SMPC:
            # non-smpc case, need to aggregate            
            for pair in list_of_xt_lists:
                XtX_list.append(pair[0])
                XtY_list.append(pair[1])
            for i in range(0, len(self.clients)):
                XtX_glob += XtX_list[i]
                XtY_glob += XtY_list[i]   
        else:
            # smpc case, already aggregated
            XtX_XtY_list = list_of_xt_lists[0]
            XtX_glob += XtX_XtY_list[0]
            XtY_glob += XtX_XtY_list[1]

        self.log("Computing beta and beta stdev...")
        for i in range(0, n):
            if linalg.det(XtX_glob[i, :, :]) == 0:
                self.log(f"XtX is singular for protein {i}, determinant is 0.")
            invXtX = linalg.inv(XtX_glob[i, :, :])
            beta[i, :] = invXtX @ XtY_glob[i, :]
            stdev_unscaled[i, :] = np.sqrt(np.diag(invXtX))

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
        self.send_data_to_coordinator([SSE, cov_coef, intensities_sum, n_measurements],
                                        send_to_self=True,
                                        use_smpc=USE_SMPC)
        
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
        list_of_sse_cov_coef = self.gather_data(use_smpc=USE_SMPC)

        k = len(self.load('variables'))
        n = len(self.load('stored_features'))

        SSE = np.zeros(n)
        cov_coef = np.zeros((k, k))
        n_measurements = np.zeros(n)
        Amean = np.zeros(n)

        if not USE_SMPC:
            # non-smpc case, need to aggregate
            SSE_list = []
            cov_coef_list = []
            n_measurements_list = []
            intensities_sum = []

            for pair in list_of_sse_cov_coef:
                SSE_list.append(pair[0])
                cov_coef_list.append(pair[1])
                n_measurements_list.append(pair[2])
                intensities_sum.append(pair[3])

            for c in range(0, len(self.clients)):
                cov_coef += cov_coef_list[c]
                Amean += intensities_sum[c]
                n_measurements += n_measurements_list[c]
                for i in range(0, n):
                    self.SSE[i] += SSE_list[c][i]   

        else:
            # smpc case, already aggregated
            sse_cov_coef = list_of_sse_cov_coef[0]
            for i in range(0, n):
                SSE[i] += sse_cov_coef[0][i]   
            cov_coef += sse_cov_coef[1]
            Amean += sse_cov_coef[2]
            n_measurements += sse_cov_coef[3]

        self.log("Aggregation of SSE is done, start computing global parameters...")
        # estimated covariance matrix of beta
        cov_coef = linalg.inv(cov_coef)
        # estimated residual variance
        var = SSE / (n_measurements - k)
        # estimated residual standard deviations
        sigma = np.sqrt(var)
        # degrees of freedom
        df_residual = np.ones(n) * (n_measurements - k)
        # mean log-intensity
        Amean = Amean / n_measurements

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
        contrasts=[([target_classes[0]], [target_classes[1]])]
        df = {}

        for contr in contrasts:
            group1, group2 = contr
            for name in group1 + group2:
                if not name in self.load('variables'):
                    raise Exception(f"{name} not found in the design matrix.")
            contr_name = "".join(map(str, group1)) + "_vs_" + "".join(map(str, group2))
            c = pd.Series(data=np.zeros(len(self.load('variables'))), index=self.load('variables'))
            c[group1] = 1
            c[group2] = -1
            df[contr_name] = c
        
        self.store(key='contrasts', value=pd.DataFrame.from_dict(df))
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
        # 	Correlation matrix of estimable coefficients
        # 	Test whether design was orthogonal
        if not np.any(self.load('cov_coef')):
            self.log(f".... no coef correlation matrix found in fit - assuming orthogonal...")
            cormatrix = np.identity(ncoef)
            orthog = True
        else:
            cormatrix = utils.cov2cor(self.load('cov_coef'))
            if cormatrix.shape[0] * cormatrix.shape[1] < 2:
                orthog = True
            else:
                if np.sum(np.abs(np.tril(cormatrix, k=-1))) < 1e-12:
                    orthog = True
                else:
                    orthog = False

        self.log(f".... orthogonality of design matrix is {orthog}")
        if np.any(np.isnan(self.load('beta'))):
            self.log(f"NA coefficients found in fit - replacing with large (but finite) standard deviations")
            self.store(key='beta', value=np.nan_to_num(self.load('beta'), nan=0))
            self.store(key='stdev_unscaled', value=np.nan_to_num(self.load('stdev_unscaled'), nan=1e30))

        self.store(key='beta', value=self.load('beta').dot(contrast_matrix))
        # New covariance coefficiets matrix
        self.store(key='cov_coef', value=contrast_matrix.T.dot(self.load('cov_coef')).dot(contrast_matrix))

        if orthog:
            self.store(key='stdev_unscaled', value=np.sqrt((self.load('stdev_unscaled')**2).dot(contrast_matrix**2)))
        else:
            n_genes = self.load('beta').shape[0]
            U = np.ones((n_genes, contrast_matrix.shape[1]))  # genes x contrasts
            o = np.ones(ncoef)
            R = np.linalg.cholesky(cormatrix).T

            for i in range(0, n_genes):
                RUC = R @ (self.load('stdev_unscaled')[i,] * contrast_matrix.T).T
                U[i,] = np.sqrt(o @ RUC**2)
            self.store(key='stdev_unscaled', value=U)
        
        self.log("Contrasts are fitted...")
        self.log("Transition to eBayes stage...")
        return 'ebayes'


@app_state(name='ebayes')
class eBayesState(AppState):
    def register(self):
        self.register_transition('get_counts', role=Role.COORDINATOR)

    def run(self):
        self.log("Start eBayes stage...")

        self.log("Calculating moderated t-statistics...")
        results, df_total = utils.moderatedT(self.load('var'), self.load('df_residual'),
                                             self.load('beta'), self.load('stdev_unscaled'))
        results["AveExpr"] = self.load('Amean')
        self.log("Result table is pre-computed...")

        self.log("Computing B statistic...")
        results = utils.Bstat(df_total, self.load('stdev_unscaled'), results,
                              stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)
        self.log("B statistic is computed...")

        self.log("Computing p-values...")        
        results = utils.topTableT(results,
                                  self.load('stored_features'), 
                                  self.load('beta'), self.load('stdev_unscaled'), df_total,
                                  adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
        self.log("P-values are computed...")
        self.store(key='results', value=results)
        self.log("Transition to getting counts...")
        return 'get_counts'


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
        self.send_data_to_coordinator(counts.to_dict(), send_to_self=True, use_smpc=USE_SMPC)
        
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
        list_of_counts = self.gather_data(use_smpc=USE_SMPC)
        min_counts = list()

        if not USE_SMPC:
            # non-smpc case, need to aggregate            
            for local_counts in list_of_counts:
                min_counts.append(pd.DataFrame.from_dict(local_counts, orient='index'))
            global_min_counts = pd.concat(min_counts, axis=1)
        else:
            # smpc case, already aggregated
            global_min_counts = pd.DataFrame.from_dict(list_of_counts[0], orient='index')

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
