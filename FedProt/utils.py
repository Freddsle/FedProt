import numpy as np
import pandas as pd

from scipy.special import digamma, polygamma
from scipy.stats import t

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import conversion, default_converter

import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

# Functions from server.py

def cov2cor(cov_coef):
    cor = np.diag(cov_coef) ** -0.5 * cov_coef
    cor = cor.T * np.diag(cov_coef) ** -0.5
    np.fill_diagonal(cor, 1)
    return cor


# eBays functions
def moderatedT(var, df_residual, beta, stdev_unscaled):
    results = squeezeVar(var, df_residual)

    # reorganize the results
    results["s2_prior"] = results["var_prior"]
    results["s2_post"] = results["var_post"]
    del results["var_prior"]
    del results["var_post"]
    results["t"] = beta / stdev_unscaled
    results["t"] = results["t"].T / np.sqrt(results["s2_post"])
    df_total = df_residual + results["df_prior"]

    df_pooled = sum(df_residual)
    df_total = np.minimum(df_total, df_pooled)  # component-wise min

    results["p_value"] = 2 * t.cdf(-np.abs(results["t"]), df=df_total)
    results["p_value"] = results["p_value"].T
    results["t"] = results["t"].T
    return results, df_total


def squeezeVar(var, df_residual):
    """Estimates df and var priors and computes posterior variances."""
    var_prior, df_prior = fitFDist(var, df_residual)

    if np.isnan(df_prior):
        logging.error("Error: Could not estimate prior df.")
        return

    var_post = posterior_var(var=var, df_residual=df_residual,
                             var_prior=var_prior, df_prior=df_prior)
    results = {
        "df_prior": df_prior,
        "var_prior": var_prior,
        "var_post": var_post,
    }
    return results


def fitFDist(var, df_residual):
    """Given x (sigma^2) and df1 (df_residual), fits x ~ scale * F(df1,df2) and returns estimated df2 and scale (s0^2)"""  
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
    evar = evar - np.nanmean(trigamma(df_residual * 1.0 / 2))
    
    if evar > 0:
        df_prior = 2 * trigammaInverse(evar)
        var_prior = np.exp(emean + digamma(df_prior * 1.0 / 2) - np.log(df_prior * 1.0 / 2))
    else:
        df_prior = np.Inf
        logging.warning("Could not estimate prior df, setting to Inf")
        var_prior = np.exp(emean)

    return var_prior, df_prior


def trigamma(x):
    """
    Calculate trigamma function using polygamma function. It is a second derivative of log gamma function
    """
    return polygamma(1, x)


def trigammaInverse(x):
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
        tri = trigamma(y)
        dif = tri * (1.0 - tri / x_) / psigamma(y, deriv=2)
        y = y + dif
        if np.max(-dif / y) < 1e-8:  # tolerance
            return y

    logging.warning("Iteration limit exceeded")
    return y


def psigamma(x, deriv=2):
    """
    Calculate polygamma function of order deriv.
    """
    return polygamma(deriv, x)


def posterior_var(var, df_residual, var_prior=np.ndarray([]), df_prior=np.ndarray([])):
    var = var
    df = df_residual
    ndxs = np.argwhere(np.isfinite(var)).reshape(-1)
    # if no infinite vars
    if len(ndxs) == len(var):  # np.isinf(df_prior).any():
        return (df * var + df_prior * var_prior) / (df + df_prior)  # var_post
    # For infinite df.prior, set var_post = var_prior
    var_post = np.repeat(var_prior, len(var))
    for ndx in ndxs:
        var_post[ndx] = (df[ndx] * var[ndx] + df_prior * var_prior) / (df[ndx] + df_prior)
    return var_post


def Bstat(df_total, stdev_unscaled, results, 
          stdev_coef_lim=np.array([0.1, 4]), proportion=0.01):
    var_prior_lim = stdev_coef_lim**2 / np.median(results["s2_prior"])

    results["var_prior"] = tmixture_matrix(df_total, stdev_unscaled, results,
                                           proportion=0.01, var_prior_lim=var_prior_lim)

    nan_ndx = np.argwhere(np.isnan(results["var_prior"]))

    if len(nan_ndx) > 0:
        results["var.prior"][nan_ndx] < -1.0 / results["s2_prior"]
        logging.warning("Estimation of var.prior failed - set to default value")
    
    r = np.outer(np.ones(results["t"].shape[0]), results["var_prior"])
    r = (stdev_unscaled**2 + r) / stdev_unscaled**2
    t2 = results["t"] ** 2

    valid_df_ndx = np.where(results["df_prior"] <= 1e6)[0]
    
    if len(valid_df_ndx) < len(results["df_prior"]):
        logging.info("Large (>1e6) priors for DF:" + str(len(valid_df_ndx)))
        kernel = t2 * (1 - 1.0 / r) / 2
        for i in valid_df_ndx:
            kernel[i] = (
                (1 + df_total[i])
                / 2
                * np.log((t2[i, :].T + df_total[i]) / ((t2[i, :] / r[i, :]).T + df_total[i]))
            )
    else:
        kernel = (1 + df_total) / 2 * np.log((t2.T + df_total) / ((t2 / r).T + df_total))

    results["B"] = np.log(proportion / (1 - proportion)) - np.log(r) / 2 + kernel.T
    return results


def tmixture_matrix(df_total, stdev_unscaled, results, 
                    var_prior_lim=False, proportion=0.01):
    tstat = results["t"].copy()
    stdev_unscaled = stdev_unscaled.copy()
    df_total = df_total.copy()
    ncoef = results["t"].shape[1]
    
    v0 = np.zeros(ncoef)
    for j in range(0, ncoef):
        v0[j] = tmixture_vector(tstat[:, j], stdev_unscaled[:, j], df_total, proportion, var_prior_lim)
    
    return v0


def tmixture_vector(tstat, stdev_unscaled, df, proportion, var_prior_lim):
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
        TailP = tail_p_calculation(tstat[i], df[i])
        tstat[i] = t_ppf_calculation(TailP, MaxDF)
        df[i] = MaxDF

    # Select top statistics
    order = tstat.argsort()[::-1][:ntarget]  # TODO: ensure the order is decreasing
    tstat = tstat[order]

    v1 = stdev_unscaled[order] ** 2

    # Compare to order statistics
    rank = np.array(range(1, ntarget + 1))
    p0 = 2 * t.sf(tstat, df=MaxDF)  # PT
    ptarget = ((rank - 0.5) / ngenes - (1.0 - p) * p0) / p
    v0 = np.zeros(ntarget)
    pos = np.where(ptarget > p0)[0]
    if len(pos) > 0:
        qtarget = -t.ppf(ptarget[pos] / 2, df=MaxDF)
        v0[pos] = v1[pos] * ((tstat[pos] / qtarget) ** 2 - 1)

    if var_prior_lim[0] and var_prior_lim[1]:
        v0 = np.minimum(np.maximum(v0, var_prior_lim[0]), var_prior_lim[1])

    return np.mean(v0)


def tail_p_calculation(tstat, df):
    with conversion.localconverter(default_converter):
        stats = importr('stats')

        tstat_r = robjects.FloatVector(tstat)
        df_r = robjects.FloatVector(df)
    
        TailP = np.array(stats.pt(tstat_r, df=df_r, lower_tail=False, log_p=True))
    logging.info("Calculating tail p-values")
    return TailP


def t_ppf_calculation(TailP, MaxDF):
    with conversion.localconverter(default_converter):
        stats = importr('stats')
        
        TailP_r = robjects.FloatVector(TailP)
        MaxDF_r = robjects.FloatVector([MaxDF] * len(TailP))
        
        tstat = np.array(stats.qt(TailP_r, df=MaxDF_r, lower_tail=False, log_p=True))
    return tstat



def topTableT(results, feature_names, beta, stdev_unscaled, df_total,
              adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95):
    results["logFC"] = pd.Series(beta[:, 0], index=feature_names)

    # confidence intervals for LogFC
    if confint:
        alpha = (1.0 + confint) / 2
        margin_error = np.sqrt(results["s2_post"]) * stdev_unscaled[:, 0] * t.ppf(alpha, df=df_total)
        results["CI.L"] = results["logFC"] - margin_error
        results["CI.R"] = results["logFC"] + margin_error
    # adjusting p-value for multiple testing
    _, adj_pval, _, _ = multipletests(
        results["p_value"][:, 0],
        alpha=p_value,
        method=adjust,
        is_sorted=False,
        returnsorted=False,
    )
    results["adj.P.Val"] = pd.Series(adj_pval, index=feature_names)
    results["P.Value"] = pd.Series(results["p_value"][:, 0], index=feature_names)
    # make table
    table = results.copy()
    # remove 'df_prior', 's2_prior', 's2_post', 'df_total','var_prior'
    for key in ["df_prior", "s2_prior", "s2_post", "var_prior", "p_value"]:
        del table[key]
    table["t"] = pd.Series(table["t"][:, 0], index=feature_names)
    table["B"] = pd.Series(table["B"][:, 0], index=feature_names)
    table = pd.DataFrame.from_dict(table)

    return table

# DEqMS functions
def spectral_count_ebayes(results, min_counts, stored_features,
                          beta, stdev_unscaled,
                          df_residual, sigma,
                          return_sorted=False,
                          fit_method="loess"):
    """
    Spectral count model for estimating variance of log-counts.
    :param fit_method: "loess" or "nls"
    :return: None
    """
    log_var = np.log(sigma**2)
    df = df_residual
    n_genes = log_var.shape[0]
    df[df == 0] = np.nan
    eg = log_var - digamma(df * 1.0 / 2) + np.log(df * 1.0 / 2)
    
    if fit_method == "loess":
        x = np.log2(min_counts)
        # y_pred = self.fit_predict_loess(x, log_var)
        logging.info("Fitting LOWESS curve...")
        logging.info(f"min_count: {min_counts.shape}, log_var: {log_var.shape}")
        y_pred = fit_LOWESS_curve(x, log_var)

    elif fit_method == "nls":
        raise NotImplementedError("NLS method not implemented yet.")
    elif fit_method == "spline":
        raise NotImplementedError("Spline method not implemented yet.")
    else:
        raise ValueError(f"Unrecognized fit method: {fit_method}")

    # egpred<-y.pred-digamma(df/2)+log(df/2)
    eg_pred = y_pred - digamma(df * 1.0 / 2) + np.log(df * 1.0 / 2)
    # myfct<- (eg-egpred)^2 - trigamma(df/2)
    myfct = (eg - eg_pred) ** 2 - trigamma(df * 1.0 / 2)

    mean_myfct = np.nanmean(myfct)
    priordf = []  # priordf <-vector()
    testd0 = []  # testd0 <-vector()

    for i in range(1, n_genes * 10 + 1):
        testd0.append(i / 10)
        priordf.append(np.abs(mean_myfct - trigamma(testd0[i-1] / 2)))
        if i > 2 and priordf[i - 3] < priordf[i - 2]:
            break

    # prior degree found
    # d0<-testd0[match(min(priordf),priordf)]
    d0 = testd0[np.argmin(priordf)]

    # calculate prior variance
    # s02<-exp(egpred + digamma(d0/2) - log(d0/2))
    s02 = np.exp(eg_pred + digamma(d0 / 2) - np.log(d0 / 2))

    post_var = (d0 * s02 + df * sigma**2) / (d0 + df)
    post_df = d0 + df

    # sca.t and scc.p stands for spectra count adjusted t and p values.
    sca_t = beta[:, 0] / (stdev_unscaled[:, 0] * np.sqrt(post_var))

    # t() is PT: CDF of t-distribution
    # 2 * pt(abs(sca.t), post.df, lower.tail = FALSE)
    #sca_p = 2 * (1 - t.cdf(np.abs(sca_t), post_df))
    sca_p = calculate_sca_t(sca_t, post_df)

    results = update_table(results, min_counts, post_df, sca_t, sca_p)
    results = apply_multiple_test_correction(results, stored_features, return_sorted=return_sorted)
    
    return results


def fit_LOWESS_curve(x, log_var, use_py=False):
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
        y_pred = loess_r_fit(log_var, x)
    
    return y_pred


def loess_r_fit(logVAR, x):
    with conversion.localconverter(default_converter):
        stats = importr('stats')

        logVAR_r = robjects.FloatVector(logVAR)
        x_r = robjects.FloatVector(x)

        formula = robjects.Formula("y ~ x")
        env = formula.environment
        env["y"] = logVAR_r
        env["x"] = x_r

        loess_model = stats.loess(formula, span=0.75)
        y_pred = np.array(stats.fitted(loess_model))

    return y_pred


def calculate_sca_t(sca_t, post_df):
        """
        Calculate spectra count adjusted t values. Using rpy2
        """
        with conversion.localconverter(default_converter):
            # Import R's stats package
            stats = importr('stats')

            sca_t_r = robjects.FloatVector(np.abs(sca_t))
            post_df_r = robjects.FloatVector(post_df)

            # 2 * pt(abs(sca.t), post.df, lower.tail = FALSE)
            # 2 * (1 - t.cdf(np.abs(sca_t), post_df))
            sca_p_r = np.array(stats.pt(sca_t_r, post_df_r, lower_tail=False))

        return 2 * sca_p_r


def update_table(table, min_counts, post_df, sca_t, sca_p):
    table["post.df"] = post_df
    table["counts"] = min_counts
    table["sca.t"] = sca_t
    table["sca.P.Value"] = sca_p
    return table


def apply_multiple_test_correction(table, stored_features, return_sorted=False):
    _, adj_pval, _, _ = multipletests(
        table["sca.P.Value"].values,
        alpha=1,
        method="fdr_bh",
        is_sorted=False,
        returnsorted=return_sorted,
    )
    table["sca.adj.pval"] = pd.Series(adj_pval, index=stored_features)

    if return_sorted:
        table.sort_values(by=["sca.adj.pval", "sca.P.Value"], inplace=True)

    return table