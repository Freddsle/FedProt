import pandas as pd
import numpy as np

# import sys
from copy import copy
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# from scipy.interpolate import interp1d
from scipy import linalg
from scipy.special import digamma, polygamma
from scipy.stats import t

from itertools import combinations

import matplotlib.pyplot as plt

import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

#pd.set_option('display.float_format', '{:.3e}'.format)


def cov2cor(cov_coef):
    cor = np.diag(cov_coef) ** -0.5 * cov_coef
    cor = cor.T * np.diag(cov_coef) ** -0.5
    np.fill_diagonal(cor, 1)
    return cor


class Server:
    def __init__(self, target_classes, covariates):
        self.target_classes = sorted(target_classes)
        self.covariates = sorted(covariates)
        self.variables = self.target_classes + self.covariates
        self.client_names = []
        self.n_samples_per_cli = []
        # self.n_tmt_per_cli = []
        self.stored_features = []

        # attributes for fedDEqMS filtering
        self.prot_na_table = pd.DataFrame()
        self.samples_per_target = None

        # count
        self.pep_counts_table = None

        # attributes for fedLmFit
        self.XtX_glob = None
        self.Xty_glob = None
        self.beta = None
        self.cov_coef = None
        self.stdev_unscaled = None
        self.var = None
        self.sigma = None
        self.df_residual = None
        self.df_total = None
        self.Amean = None
        self.results = None
        self.table = None

        # attributes for fed median
        self.lb = None
        self.ub = None
        self.approx_median = None  # approximate median found in binary search
        self.tol = 0.1  # binary search of approx_medina stops when ub-lb<tol
        self.is_converged = None  # whether approx_median is found for each gene
        self.le_gt = None  # dataframe with numbers of values lesser or equal and greater than approx_median. "le":<= , "gt":>
        self.precise_median = None

    def join_client(self, client_cohort_name, client_prot_names, client_n_samples):
        """
        Collects names of genes and variables, and the number of samples from client.
        """
        if client_cohort_name in self.client_names:
            logging.error(f"Choose client name other than {self.client_names}")
            logging.error(f"Failed to join client {client_cohort_name}")
            return False

        # keep only shared features - intersection of features from all clients
        if len(self.client_names) == 0:
            self.stored_features = sorted(set(client_prot_names))
        else:
            self.stored_features = sorted(list(set(self.stored_features) & set(client_prot_names)))

        self.n_samples_per_cli.append(client_n_samples)
        # self.n_tmt_per_cli.append(client.n_tmt)
        self.client_names.append(client_cohort_name)
        logging.info(f"Server: joined client  {client_cohort_name}")

        # self.prots_na_table = self.prots_na_table.loc[self.stored_features, :]

        return True

    def filter_by_median(self, shuffled_count_matrices, threshold=1, func=np.median):
        """
        Accepts a list of matrices with shuffled expressions and computes medians.
        """
        medians = pd.concat(shuffled_count_matrices, axis=1).apply(lambda x: func(x), axis=1)
        keep_genes = list(medians[medians > threshold].index.values)  # send to clients

        logging.info(
            f"Genes with median counts < {threshold} will be dropped: {medians.shape[0] - len(keep_genes)}"
        )  # with  median  CPM > 1
        return keep_genes

    def compute_F(self, f):
        F = np.exp(np.mean(np.log(f)))
        return F

    def create_na_df(self, na_count_in_variable, samples_per_target):
        """
        Create the df with the number of NA value per target class
        """
        if self.prot_na_table.empty:
            self.prot_na_table = na_count_in_variable
            self.samples_per_target = samples_per_target

        else:
            # if already exisct - uodate df and dict
            for key in self.samples_per_target:
                self.samples_per_target[key] = self.samples_per_target[key] + samples_per_target.get(key)
            self.prot_na_table = self.prot_na_table.add(na_count_in_variable)

    def update_prot_names(self, min_f):
        """
        Updates shared protein list using data from the clients after filtering
        """

        #caclulate the percentage of NA values per target class
        self.prot_na_table = self.prot_na_table.loc[:, self.samples_per_target.keys()]
        samples_series = pd.Series(self.samples_per_target)
        na_perc = self.prot_na_table.div(samples_series, axis=1)

        keep_proteins = na_perc[na_perc.apply(lambda row: all(row < min_f), axis=1)].index.values

        # updated shared protein list
        self.stored_features = sorted(list(set(keep_proteins)))
        return self.stored_features

    ###### fedLmFit #####
    def compute_beta_and_beta_stdev(self, XtX_list, XtY_list):
        """Calcualtes global beta and variance of beta"""
        k = len(self.variables)
        n = len(self.stored_features)
        self.XtX_glob = np.zeros((n, k, k))
        self.XtY_glob = np.zeros((n, k))
        self.stdev_unscaled = np.zeros((n, k))

        logging.info(f"Server: computing global beta and beta stdev, k = {k}, n = {n}")

        for i in range(0, len(self.client_names)):
            self.XtX_glob += XtX_list[i]
            self.XtY_glob += XtY_list[i]
        self.beta = np.zeros((n, k))
        self.rank = np.ones(n) * k

        for i in range(0, n):
            invXtX = linalg.inv(self.XtX_glob[i, :, :])
            self.beta[i, :] = invXtX @ self.XtY_glob[i, :]
            self.stdev_unscaled[i, :] = np.sqrt(np.diag(invXtX))  # standart err for b coefficients

    def aggregate_SSE_and_cov_coef(self, SSE_list, cov_coef_list, 
                                   intensities_sum, n_measurements_list):
        n = len(self.stored_features)  # proteins
        k = len(self.variables)
        self.SSE = np.zeros(n)
        self.cov_coef = np.zeros((k, k))
        self.n_measurements = np.zeros(n)  # number of not nan measurements per protein
        self.Amean = np.zeros(n)

        for c in range(0, len(self.client_names)):
            self.cov_coef += cov_coef_list[c]
            self.Amean += intensities_sum[c]
            self.n_measurements += n_measurements_list[c]
            for i in range(0, n):
                self.SSE[i] += SSE_list[c][i]

        self.cov_coef = linalg.inv(self.cov_coef)
        # estimated residual variance
        self.var = self.SSE / (self.n_measurements - k)
        # estimated residual standard deviations
        self.sigma = np.sqrt(self.var)
        # degrees of freedom
        self.df_residual = np.ones(n) * (self.n_measurements - k)
        # mean log-intensity
        self.Amean = self.Amean / self.n_measurements

    ### apply contrasts

    def make_contrasts(self, contrasts=[]):
        """Creates contrast matrix given deisgn matrix and pairs or columns to compare.\n
        For example:\n
        contrasts = [([A],[B]),([A,B],[C,D])] defines two contrasts:\n
        A-B and (A and B) - (C and D)."""
        df = {}

        for contr in contrasts:
            group1, group2 = contr
            for name in group1 + group2:
                if not name in self.variables:
                    raise Exception(f"{name} not found in the design matrix.")
            contr_name = "".join(map(str, group1)) + "_vs_" + "".join(map(str, group2))
            c = pd.Series(data=np.zeros(len(self.variables)), index=self.variables)
            c[group1] = 1
            c[group2] = -1
            df[contr_name] = c
        return pd.DataFrame.from_dict(df)

    def fit_contasts(self, contrast_matrix):
        ncoef = self.cov_coef.shape[1]
        # 	Correlation matrix of estimable coefficients
        # 	Test whether design was orthogonal
        if not np.any(self.cov_coef):
            logging.info(f"no coef correlation matrix found in fit - assuming orthogonal")
            cormatrix = np.identity(ncoef)
            orthog = True
        else:
            cormatrix = cov2cor(self.cov_coef)
            if cormatrix.shape[0] * cormatrix.shape[1] < 2:
                orthog = True
            else:
                if np.sum(np.abs(np.tril(cormatrix, k=-1))) < 1e-12:
                    orthog = True
                else:
                    orthog = False
        # print("is design orthogonal:",orthog)
        # 	Replace NA coefficients with large (but finite) standard deviations
        # 	to allow zero contrast entries to clobber NA coefficients.
        if np.any(np.isnan(self.beta)):
            logging.info(f"NA coefficients found in fit - replacing with large (but finite) standard deviations")
            np.nan_to_num(self.beta, nan=0)
            np.nan_to_num(self.stdev_unscaled, nan=1e30)

        self.beta = self.beta.dot(contrast_matrix)
        # New covariance coefficiets matrix
        self.cov_coef = contrast_matrix.T.dot(self.cov_coef).dot(contrast_matrix)

        if orthog:
            self.stdev_unscaled = np.sqrt((self.stdev_unscaled**2).dot(contrast_matrix**2))
        else:
            n_genes = self.beta.shape[0]
            U = np.ones((n_genes, contrast_matrix.shape[1]))  # genes x contrasts
            o = np.ones(ncoef)
            R = np.linalg.cholesky(cormatrix).T

            for i in range(0, n_genes):
                RUC = R @ (self.stdev_unscaled[i,] * contrast_matrix.T).T
                U[i,] = np.sqrt(o @ RUC**2)
            self.stdev_unscaled = U

    #### e-Bayes ############

    def trigamma(self, x):
        """
        Calculate trigamma function using polygamma function. It is a second derivative of log gamma function
        """
        return polygamma(1, x)

    def psigamma(self, x, deriv=2):
        """
        Calculate polygamma function of order deriv.
        """
        return polygamma(deriv, x)

    def trigammaInverse(self, x):
        """
        Define inverse of trigamma function using Newton's method.
        """
        if not hasattr(x, "__iter__"):
            x_ = np.array([x])
        for i in range(0, x_.shape[0]):
            if np.isnan(x_[i]):
                x_[i] = np.NaN
            elif x_[i] > 1e7:
                x_[i] = 1.0 / np.sqrt(x[i])
            elif x_[i] < 1e-6:
                x_[i] = 1.0 / x[i]
        # Newton's method
        y = 0.5 + 1.0 / x_
        for i in range(0, 50):
            tri = self.trigamma(y)
            dif = tri * (1.0 - tri / x_) / self.psigamma(y, deriv=2)
            y = y + dif
            if np.max(-dif / y) < 1e-8:  # tolerance
                return y

        logging.warning("Iteration limit exceeded")
        return y

    def fitFDist(self, var, df_residual, covariate=False):
        """Given x (sigma^2) and df1 (df_residual), fits x ~ scale * F(df1,df2) and returns estimated df2 and scale (s0^2)"""
        if covariate:
            # TBD
            logging.info(f"Covariate not implemented: Set covariate=False.")
            return        
        # Avoid negative variances
        var = [max(x, 0) for x in var]
        # calculates the median of the variances
        m = np.median(var)
        if m == 0:
            logging.warning(f"More than half of residual variances are exactly zero: eBayes unreliable")
            m = 1
        else:
            if 0 in var:
                logging.warning(f"Zero sample variances detected, have been offset (+1e-5) away from zero")
        
                # avoid very small variances (close to zero)
        var = [max(x, 1e-5 * m) for x in var]

        z = np.log(var)
        e = z - digamma(df_residual * 1.0 / 2) + np.log(df_residual * 1.0 / 2)
        emean = np.nanmean(e)
        evar = np.nansum((e - emean) ** 2) / (len(var) - 1)

        # Estimate scale and df2
        evar = evar - np.nanmean(self.trigamma(df_residual * 1.0 / 2))
        
        if evar > 0:
            df_prior = 2 * self.trigammaInverse(evar)
            var_prior = np.exp(emean + digamma(df_prior * 1.0 / 2) - np.log(df_prior * 1.0 / 2))
        else:
            df2 = np.Inf
            var_prior = np.exp(emean)

        return var_prior, df_prior

    def posterior_var(self, var_prior=np.ndarray([]), df_prior=np.ndarray([])):
        var = self.var.values
        df = self.df_residual.values
        ndxs = np.argwhere(np.isfinite(var)).reshape(-1)
        # if no infinite vars
        if len(ndxs) == len(var):  # np.isinf(df_prior).any():
            return (df * var + df_prior * var_prior) / (df + df_prior)  # var_post
        # For infinite df.prior, set var_post = var_prior
        var_post = np.repeat(var_prior, len(var))
        for ndx in ndxs:
            var_post[ndx] = (df[ndx] * var[ndx] + df_prior * var_prior) / (df[ndx] + df_prior)
        return var_post

    def squeezeVar(self, covariate=False, robust=False, winsor_tail_p=(0.05, 0.1)):
        """Estimates df and var priors and computes posterior variances."""
        if robust:
            # TBD fitFDistRobustly()
            logging.info("Set robust=False.")
            return
        else:
            var_prior, df_prior = self.fitFDist(self.var, self.df_residual, covariate=covariate)

        if np.isnan(df_prior):
            logging.error("Error: Could not estimate prior df.")
            return

        var_post = self.posterior_var(var_prior=var_prior, df_prior=df_prior)
        self.results = {
            "df_prior": df_prior,
            "var_prior": var_prior,
            "var_post": var_post,
        }

    def moderatedT(self, covariate=False, robust=False, winsor_tail_p=(0.05, 0.1)):
        # var,df_residual,coefficients,stdev_unscaled,
        self.squeezeVar(covariate=covariate, robust=robust, winsor_tail_p=winsor_tail_p)

        self.results["s2_prior"] = self.results["var_prior"]
        self.results["s2_post"] = self.results["var_post"]
        del self.results["var_prior"]
        del self.results["var_post"]
        self.results["t"] = self.beta / self.stdev_unscaled
        self.results["t"] = self.results["t"].T / np.sqrt(self.results["s2_post"])
        self.df_total = self.df_residual.values + self.results["df_prior"]

        df_pooled = sum(self.df_residual)
        self.df_total = np.minimum(self.df_total, df_pooled)  # component-wise min

        self.results["p_value"] = 2 * t.cdf(-np.abs(self.results["t"]), df=self.df_total)
        self.results["p_value"] = self.results["p_value"].T
        self.results["t"] = self.results["t"].T
        return self.results

    def tmixture_matrix(self, var_prior_lim=False, proportion=0.01):
        tstat = self.results["t"].copy()
        stdev_unscaled = self.stdev_unscaled.copy()
        df_total = self.df_total.copy()
        ncoef = self.results["t"].shape[1]
        v0 = np.zeros(ncoef)
        for j in range(0, ncoef):
            v0[j] = self.tmixture_vector(tstat[:, j], stdev_unscaled[:, j], df_total, proportion, var_prior_lim)
        return v0

    def tail_p_calculation(self, tstat, df):
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr

        stats = importr('stats')
        
        tstat_r = robjects.FloatVector(tstat)
        df_r = robjects.FloatVector(df)
        
        TailP = stats.pt(tstat_r, df=df_r, lower_tail=False, log_p=True)
        logging.info("Calculating tail p-values")
        return np.array(TailP)

    def t_ppf_calculation(self, TailP, MaxDF):
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr

        stats = importr('stats')
        
        TailP_r = robjects.FloatVector(TailP)
        MaxDF_r = robjects.FloatVector([MaxDF] * len(TailP))
        
        tstat = stats.qt(TailP_r, df=MaxDF_r, lower_tail=False, log_p=True)
        return np.array(tstat)

    def tmixture_vector(self, tstat, stdev_unscaled, df, proportion, var_prior_lim):
        ngenes = len(tstat)
        # Remove missing values
        notnan_ndx = np.where(~np.isnan(tstat))[0]
        if len(notnan_ndx) < ngenes:
            tstat = tstat[notnan_ndx]
            stdev_unscaled = stdev_unscaled[notnan_ndx]
            df = df[notnan_ndx]

        # ntarget t-statistics will be used for estimation
        ntarget = int(np.ceil(proportion / 2 * ngenes))
        if ntarget < 1:  #
            return

        # If ntarget is v small, ensure p at least matches selected proportion
        # This ensures ptarget < 1
        p = np.maximum(ntarget * 1.0 / ngenes, proportion)
        # Method requires that df be equal
        tstat = abs(tstat)
        MaxDF = np.max(df)
        i = np.where(df < MaxDF)[0]

        if len(i) > 0:
            # epsilon = 1e-30
            # TailP = np.log(np.clip(1 - t.cdf(tstat[i], df=df[i]), epsilon, 1))
            TailP = self.tail_p_calculation(tstat[i], df[i])
            #TailP = t.logcdf(tstat[i], df=df[i])  
            # PT: CDF of t-distribution: pt(tstat[i],df=df[i],lower.tail=FALSE,log.p=TRUE)
            # QT - qunatile funciton - returns a threshold value x
            # below which random draws from the given CDF would fall p percent of the time. [wiki]
            #tstat[i] = t.ppf(np.exp(TailP), df=MaxDF)  # QT: qt(TailP,df=MaxDF,lower.tail=FALSE,log.p=TRUE)

            tstat[i] = self.t_ppf_calculation(TailP, MaxDF)
            df[i] = MaxDF

        # Select top statistics
        order = tstat.argsort()[::-1][:ntarget]  # TBD: ensure the order is decreasing
        tstat = tstat[order]

        v1 = stdev_unscaled[order] ** 2

        # Compare to order statistics
        rank = np.array(range(1, ntarget + 1))
        p0 = 2 * t.sf(tstat, df=MaxDF)  # PT
        ptarget = ((rank - 0.5) / ngenes - (1.0 - p) * p0) / p
        v0 = np.zeros(ntarget)
        pos = np.where(ptarget > p0)[0]
        if len(pos) > 0:
            qtarget = -t.ppf(ptarget[pos] / 2, df=MaxDF)  # qt(ptarget[pos]/2,df=MaxDF,lower.tail=FALSE)
            # print(qtarget[:5])
            v0[pos] = v1[pos] * ((tstat[pos] / qtarget) ** 2 - 1)

        if var_prior_lim[0] and var_prior_lim[1]:
            v0 = np.minimum(np.maximum(v0, var_prior_lim[0]), var_prior_lim[1])
        return np.mean(v0)

    def Bstat(self, stdev_coef_lim=np.array([0.1, 4]), proportion=0.01):
        var_prior_lim = stdev_coef_lim**2 / np.median(self.results["s2_prior"])
        # print("Limits for var.prior:",var_prior_lim)

        self.results["var_prior"] = self.tmixture_matrix(proportion=0.01, var_prior_lim=var_prior_lim)

        nan_ndx = np.argwhere(np.isnan(self.results["var_prior"]))

        if len(nan_ndx) > 0:
            self.results["var.prior"][nan_ndx] < -1.0 / self.results["s2_prior"]
            logging.warning("Estimation of var.prior failed - set to default value")
        
        r = np.outer(np.ones(self.results["t"].shape[0]), self.results["var_prior"])
        r = (self.stdev_unscaled**2 + r) / self.stdev_unscaled**2
        t2 = self.results["t"] ** 2

        valid_df_ndx = np.where(self.results["df_prior"] <= 1e6)[0]
        
        if len(valid_df_ndx) < len(self.results["df_prior"]):
            logging.info("Large (>1e6) priors for DF:" + str(len(valid_df_ndx)))
            kernel = t2 * (1 - 1.0 / r) / 2
            for i in valid_df_ndx:
                kernel[i] = (
                    (1 + self.df_total[i])
                    / 2
                    * np.log((t2[i, :].T + self.df_total[i]) / ((t2[i, :] / r[i, :]).T + self.df_total[i]))
                )
        else:
            kernel = (1 + self.df_total) / 2 * np.log((t2.T + self.df_total) / ((t2 / r).T + self.df_total))

        self.results["B"] = np.log(proportion / (1 - proportion)) - np.log(r) / 2 + kernel.T

    def topTableT(self, adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95):
        feature_names = self.stored_features
        self.results["logFC"] = pd.Series(self.beta[:, 0], index=feature_names)

        # confidence intervals for LogFC
        if confint:
            alpha = (1.0 + confint) / 2
            margin_error = np.sqrt(self.results["s2_post"]) * self.stdev_unscaled[:, 0] * t.ppf(alpha, df=self.df_total)
            self.results["CI.L"] = self.results["logFC"] - margin_error
            self.results["CI.R"] = self.results["logFC"] + margin_error
        # adjusting p-value for multiple testing
        _, adj_pval, _, _ = multipletests(
            self.results["p_value"][:, 0],
            alpha=p_value,
            method=adjust,
            is_sorted=False,
            returnsorted=False,
        )
        self.results["adj.P.Val"] = pd.Series(adj_pval, index=feature_names)
        self.results["P.Value"] = pd.Series(self.results["p_value"][:, 0], index=feature_names)
        # make table
        self.table = self.results.copy()
        # remove 'df_prior', 's2_prior', 's2_post', 'df_total','var_prior'
        for key in ["df_prior", "s2_prior", "s2_post", "var_prior", "p_value"]:
            del self.table[key]
        self.table["t"] = pd.Series(self.table["t"][:, 0], index=feature_names)
        self.table["B"] = pd.Series(self.table["B"][:, 0], index=feature_names)
        self.table = pd.DataFrame.from_dict(self.table)

    def eBayes(self):
        covariate = False  # Amean for limma-trend
        robust = False  #
        winsor_tail_p = (0.05, 0.1)  # needed for fitFDistRobustly()

        self.results = self.moderatedT(covariate=covariate, robust=robust, winsor_tail_p=winsor_tail_p)
        # self.results = moderatedT(self.var,self.df_residual,self.beta,self.stdev_unscaled,
        #                         covariate=covariate,robust=robust, winsor_tail_p=winsor_tail_p)
        self.results["AveExpr"] = self.Amean

        self.Bstat(stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)

        self.topTableT(adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
        # self.table = self.table.sort_values(by="P.Value")

    # def plot_voom_results(self, svg_file=""):
    #     """Plots relationships between average log-count and square root of estimated std.dev for each gene."""
    #     fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    #     axes[0].set_title("Before Voom")
    #     axes[0].scatter(x=self.mean_logC, y=self.sigma**0.5, s=0.5, color="black")
    #     axes[0].set_xlabel(r"average($log_2$(counts + 0.5)) per gene")
    #     axes[0].set_ylabel(r"$\sqrt{s_g}$ - root of estimated std.dev for gene g")
    #     axes[0].scatter(self.lowess_curve[:, 0], self.lowess_curve[:, 1], s=0.5, color="red")

    #     axes[1].set_title("After Voom")
    #     axes[1].scatter(x=self.mean_logC_w, y=self.sigma_w**0.5, s=0.5, color="black")
    #     axes[1].set_xlabel(r"average($log_2$(counts + 0.5)) per gene")
    #     # axes[1].set_ylabel(r'$\sqrt{s_g}$ - root of estimated std.dev for gene g')
    #     axes[1].scatter(self.lowess_curve_w[:, 0], self.lowess_curve_w[:, 1], s=0.5, color="red")
    #     if svg_file:
    #         plt.savefig(svg_file)

    ##### DEqMS #####
    def loess_r_fit(self, logVAR, x):
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr

        stats = importr('stats')

        logVAR_r = robjects.FloatVector(logVAR)
        x_r = robjects.FloatVector(x)

        formula = robjects.Formula("y ~ x")
        env = formula.environment
        env["y"] = logVAR_r
        env["x"] = x_r

        loess_model = stats.loess(formula, span=0.75)
        y_pred = stats.fitted(loess_model)
        return np.array(y_pred)

    def fit_LOWESS_curve(self, x, log_var, use_py=False):
        """
        Fit LOWESS (Locally Weighted Scatterplot Smoothing) curve.
        :param x: input data
        :param log_var: logarithm of the variable
        :return: 2nd column of LOWESS output
        """
        if use_py:
            # Python implementation
            delta = (max(x) - min(x)) * 0.01
            logging.info(f"server: delta for LOWESS fit is set to: {delta}")
            lowess = sm.nonparametric.lowess        
            y_pred = lowess(log_var, x, frac=0.75, delta=delta, return_sorted=False, is_sorted=False)
        else:
            # R implementation
            y_pred = self.loess_r_fit(log_var, x)
        
        return y_pred
   
    def calculate_sca_t(self, sca_t, post_df):
        """
        Calculate spectra count adjusted t values. Using rpy2
        """
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr

        # Import R's stats package
        stats = importr('stats')

        sca_t_r = robjects.FloatVector(np.abs(sca_t))
        post_df_r = robjects.FloatVector(post_df)

        # 2 * pt(abs(sca.t), post.df, lower.tail = FALSE)
        # 2 * (1 - t.cdf(np.abs(sca_t), post_df))
        sca_p_r = stats.pt(sca_t_r, post_df_r, lower_tail=False)

        return 2 * np.array(sca_p_r)


    def spectral_count_ebayes(self, fit_method="loess", return_sorted=False):
        """
        Spectral count model for estimating variance of log-counts.
        :param fit_method: "loess" or "nls"
        :return: None
        """
        log_var = np.log(self.sigma**2)
        df = self.df_residual
        n_genes = log_var.shape[0]
        df[df == 0] = np.nan
        eg = log_var - digamma(df * 1.0 / 2) + np.log(df * 1.0 / 2)
        
        if fit_method == "loess":
            x = np.log2(self.min_counts)
            # y_pred = self.fit_predict_loess(x, log_var)
            logging.info("Fitting LOWESS curve...")
            logging.info(f"min_count: {self.min_counts.shape}, log_var: {log_var.shape}")
            y_pred = self.fit_LOWESS_curve(x, log_var)
        
        elif fit_method == "nls":
            raise NotImplementedError("NLS method not implemented yet.")
        elif fit_method == "spline":
            raise NotImplementedError("Spline method not implemented yet.")
        else:
            raise ValueError(f"Unrecognized fit method: {fit_method}")

        # egpred<-y.pred-digamma(df/2)+log(df/2)
        eg_pred = y_pred - digamma(df * 1.0 / 2) + np.log(df * 1.0 / 2)
        # myfct<- (eg-egpred)^2 - trigamma(df/2)
        myfct = (eg - eg_pred) ** 2 - self.trigamma(df * 1.0 / 2)

        mean_myfct = np.nanmean(myfct)
        priordf = []  # priordf <-vector()
        testd0 = []  # testd0 <-vector()

        for i in range(1, n_genes * 10 + 1):
            testd0.append(i / 10)
            priordf.append(np.abs(mean_myfct - self.trigamma(testd0[i-1] / 2)))
            if i > 2 and priordf[i - 3] < priordf[i - 2]:
                break

        # prior degree found
        # d0<-testd0[match(min(priordf),priordf)]
        d0 = testd0[np.argmin(priordf)]

        # calculate prior variance
        # s02<-exp(egpred + digamma(d0/2) - log(d0/2))
        s02 = np.exp(eg_pred + digamma(d0 / 2) - np.log(d0 / 2))

        post_var = (d0 * s02 + df * self.sigma**2) / (d0 + df)
        post_df = d0 + df

        # sca.t and scc.p stands for spectra count adjusted t and p values.
        sca_t = self.beta[:, 0] / (self.stdev_unscaled[:, 0] * np.sqrt(post_var))

        # t() is PT: CDF of t-distribution
        # 2 * pt(abs(sca.t), post.df, lower.tail = FALSE)
        #sca_p = 2 * (1 - t.cdf(np.abs(sca_t), post_df))
        sca_p = self.calculate_sca_t(sca_t, post_df)

        self.update_table(post_df, sca_t, sca_p)
        self.apply_multiple_test_correction(return_sorted=False)

    def update_table(self, post_df, sca_t, sca_p):
        self.table["post.df"] = post_df
        self.table["counts"] = self.min_counts
        self.table["sca.t"] = sca_t
        self.table["sca.P.Value"] = sca_p

    def apply_multiple_test_correction(self,  return_sorted=False):
        _, adj_pval, _, _ = multipletests(
            self.table["sca.P.Value"].values,
            alpha=1,
            method="fdr_bh",
            is_sorted=False,
            returnsorted=False,
        )
        self.table["sca.adj.pval"] = pd.Series(adj_pval, index=self.stored_features)
        if return_sorted:
            self.table.sort_values(by=["sca.adj.pval", "sca.P.Value"], inplace=True)

    # DEqMS pep counts
    def create_pep_counts(self, client_pep_counts):
        """Creates a table with unique protein-peptide rows from all clients."""
        if self.pep_counts_table is None:
            self.pep_counts_table = client_pep_counts
        else:
            self.pep_counts_table = pd.concat([self.pep_counts_table, client_pep_counts], axis=0)
            self.pep_counts_table = self.pep_counts_table.drop_duplicates()
            self.pep_counts_table = self.pep_counts_table.reset_index(drop=True)

    def summarize_pep_counts(self):
        """Summarizes peptide counts for each protein group to get protein counts."""
        self.min_counts = self.pep_counts_table.groupby('Protein.Group')['Precursor.Id'].nunique().reset_index().set_index('Protein.Group')
        self.min_counts = self.min_counts.loc[self.stored_features,:]
        # tranform to series
        self.min_counts = self.min_counts.iloc[:, 0]


    # def variance_deqms_plot(self, n=20, xlab="count", ylab="log(Variance)", main=''):   
    #     #plotting
    #     x = self.min_counts
    #     y = np.log(self.sigma**2)

    #     df_temp = pd.DataFrame({"pep_count": x, "variance": y})
    #     df_temp_filter = df_temp[df_temp["pep_count"] <= n]

    #     x_values = df_temp_filter["pep_count"]
    #     y_pred = self.fit_LOWESS_curve(x_values, df_temp_filter["variance"])
    #     print(y_pred)
    #     fig, ax = plt.subplots(figsize=(10, 6))

    #     boxplot_data = [group["variance"].tolist() for _, group in df_temp_filter.groupby("pep_count")]
    #     ax.boxplot(boxplot_data, positions=df_temp_filter["pep_count"].unique(), widths=0.6)
    #     ax.set_xticklabels(df_temp_filter["pep_count"].unique())

    #     ax.set_xlabel(xlab)
    #     ax.set_ylabel(ylab)
    #     ax.set_title(main)
    #     ax.plot(x_values, y_pred, color='red', linewidth=3)
        
    #     return fig, ax