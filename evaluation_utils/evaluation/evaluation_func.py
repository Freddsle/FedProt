import pandas as pd
import numpy as np 
from statsmodels.stats.multitest import multipletests

from scipy.stats import mannwhitneyu

import matplotlib.pyplot as plt
from matplotlib.table import table
from matplotlib.colors import LinearSegmentedColormap

import logging


def read_results(workdir, 
                fedprot_name="/DPE.csv", 
                deqms_name="/Central_res_irs_on_median.tsv",
                fisher_name="/MA_CM.tsv",
                rem_name="/MA_REM.tsv",
                stouffer_name="/MA_Stouffer.tsv",   
                rankprod_name="/MA_RankProd.tsv",
                corrected_deqms_name=None,
                only_two = False,
                drop_na=True, 
                simulated=False):
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    df = {}
    
    # Central analysis
    DEqMS = pd.read_csv(workdir+deqms_name, sep="\t", index_col=0)
    df["pv_DEqMS"] = DEqMS["sca.adj.pval"]
    # df["pv_DEqMS"] = DEqMS["sca.P.Value"]
    df["lfc_DEqMS"] = DEqMS["logFC"]
    df["AvgExpr_DEqMS"] = DEqMS["AveExpr"]
    logging.info(f"Results loaded for DEqMS with {DEqMS.shape[0]} proteins.")

    # FedProt
    fedprot = pd.read_csv(workdir+fedprot_name, sep="\t", index_col=0)
    df["pv_FedProt"] = fedprot["sca.adj.pval"]
    # df["pv_FedProt"] = fedprot["sca.P.Value"]
    df["lfc_FedProt"] = fedprot["logFC"]
    logging.info(f"Results loaded for FedProt with {fedprot.shape[0]} proteins.")

    if only_two:
        if corrected_deqms_name:
            DEqMS_corrected = pd.read_csv(workdir+corrected_deqms_name, sep="\t", index_col=0)
            df["pv_DEqMS_corrected"] = DEqMS_corrected["sca.adj.pval"]
            df["lfc_DEqMS_corrected"] = DEqMS_corrected["logFC"]
            logging.info(f"Results loaded for DEqMS_corrected with {DEqMS_corrected.shape[0]} proteins.")
        
        df = pd.DataFrame.from_dict(df)
        if drop_na:
            df = df.dropna(axis=0)
        logging.info(f"Results loaded from {workdir} with {df.shape[0]} genes. Adj.p-values were not log-transformed.")
        return df

    # Fisher
    ma_cm = pd.read_csv(workdir+fisher_name, sep="\t")
    if simulated:
        ma_cm.index = ma_cm["Symbol"].values
    else:
        ma_cm.index = ma_cm["ID"].values
    df["lfc_Fisher"] = ma_cm["metafc"]
    _, adj_pval,_,_ = multipletests(ma_cm["metap"].values, alpha=0.05, method='fdr_bh',
                                    is_sorted=False, returnsorted=False)
    df["pv_Fisher"] = pd.Series(adj_pval,index=ma_cm["metap"].index)
    logging.info(f"Results loaded for Fisher with {ma_cm.shape[0]} proteins.")

    # REM
    ma_rem = pd.read_csv(workdir+rem_name, sep="\t")
    if simulated:
        ma_rem.index = ma_rem["Symbol"].values
    else:       
        ma_rem.index = ma_rem["ID"].values
    df["lfc_REM"] = ma_rem["randomSummary"]
    _, adj_pval, _, _ = multipletests(ma_rem["randomP"].values, alpha=0.05, method='fdr_bh',
                                      is_sorted=False, returnsorted=False)
    df["pv_REM"] = pd.Series(adj_pval,index=ma_rem["randomP"].index)
    logging.info(f"Results loaded for REM with {ma_rem.shape[0]} proteins.")

    ### Stoufer 
    if simulated:
        stoufer  = pd.read_csv(workdir+stouffer_name, sep="\t", index_col=0)
    else:
        stoufer  = pd.read_csv(workdir+stouffer_name, sep="\t")
        stoufer.index = stoufer["ID"].values
    df["pv_Stouffer"] = stoufer["FDR"]
    df["lfc_Stouffer"] = df["lfc_Fisher"]  # take logFC from MetaVolcanoR
    logging.info(f"Results loaded for Stouffer with {stoufer.shape[0]} proteins.")

    ### RankProd
    if simulated:
        rankprod  = pd.read_csv(workdir+rankprod_name, sep="\t", index_col=0)
    else:
        rankprod  = pd.read_csv(workdir+rankprod_name, sep="\t")
        rankprod.index = rankprod["ID"].values
    rankprod["FDR"] = rankprod.loc[:,["down_reg.FDR","up_reg.FDR"]].min(axis=1)
    df["pv_RankProd"] = rankprod["FDR"]
    df["lfc_RankProd"] = rankprod["avgL2FC"] 
    logging.info(f"Results loaded for RankProd with {rankprod.shape[0]} proteins.")

    df = pd.DataFrame.from_dict(df)
    if drop_na:
        df = df.dropna(axis=0)

    # warning if adj.p-values are smaller than 0
    logging.info(f"Results loaded from {workdir} with {df.shape[0]} genes. Adj.p-values were not log-transformed.")
    return df


def calculate_correlations(df, methods, top_genes, column_name="pv_"):
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    if column_name == "pv_":
        max_p_val = np.max(np.abs(df['pv_DEqMS']))
        logging.info(f"Calculating corrs. Using p-vals - {'log-transformed' if max_p_val > 1.1 else 'not log-transformed'}.")

    correlations = {}
    df_sorted = df.sort_values(by="pv_DEqMS", ascending=False).head(top_genes)
    pearson_corr = df_sorted[[column_name + "DEqMS"] + [column_name + m for m in methods]].corr(
        ).loc[[column_name + "DEqMS"],]
    spearman_corr = df_sorted[[column_name + "DEqMS"] + [column_name + m for m in methods]].corr(method="spearman"
        ).loc[[column_name + "DEqMS"],]

    pearson_corr.rename(lambda x: x.replace(column_name, ""), axis="columns", inplace=True)
    spearman_corr.rename(lambda x: x.replace(column_name, ""), axis="columns", inplace=True)

    correlations['r'] = pearson_corr.T[column_name+"DEqMS"]
    correlations['ρ'] = spearman_corr.T[column_name+"DEqMS"]
    
    logging.info(f"Correlations computed for {'all' if top_genes == -1 else top_genes} genes from {column_name} columns.")
    return correlations


def calculate_rmse(df, methods, top_genes, column_name="pv_"):
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    if column_name == "pv_":
        max_p_val = np.max(np.abs(df['pv_DEqMS']))
        logging.info(f"Calculating RMSE. Using p-vals - {'log-transformed' if max_p_val > 1.1 else 'not log-transformed'}.")

    rmse_results = {}
    df_sorted = df.sort_values(by="pv_DEqMS", ascending=False).head(top_genes)
    for m in methods:        
        x = df_sorted[column_name + "DEqMS"].values
        y = df_sorted[column_name + m].values

        rmse = np.sqrt(np.sum((x - y) ** 2) / len(x))

        # calculate min max variant of RMSE
        # nrmse = rmse / np.var(x)
        nrmse = rmse / (np.max(x) - np.min(x))
        rmse_results[m] = {"RMSE": rmse, "NRMSE": nrmse}

    logging.info(f"RMSE and NRMSE computed for {'all' if top_genes == -1 else top_genes} genes from {column_name} columns.")
    return rmse_results


def calculate_performance_metrics(
    df, de, 
    methods, 
    lfc_thr, adj_pval_thr, 
    top_genes, all_genes
    ):

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    max_p_val = np.max(np.abs(de['pv_DEqMS']))
    if max_p_val < 1.1:
        logging.error("Adj.p-values are not log-transformed. Please, log-transform them before calculating performance metrics.")

    results = {}
    T = set(de.index.values)
    F = all_genes.difference(T)

    for m in methods: 
        de2 = df.loc[:, ["pv_" + m, "lfc_" + m]].sort_values(by="pv_" + m, ascending=False)
        de2 = de2.loc[(de2["pv_" + m] > adj_pval_thr) & (np.abs(de2["lfc_" + m]) >= lfc_thr), :]
        
        de2 = de2.head(top_genes if top_genes > 0 else de2.shape[0])

        P = set(de2.index.values)
        N = all_genes.difference(P)

        TP = len(T.intersection(P))
        FP = len(F.intersection(P))
        TN = len(F.intersection(N))
        FN = len(T.intersection(N))

        Prec = TP / (TP + FP) if (TP + FP) > 0 else 0
        Rec = TP / (TP + FN) if (TP + FN) > 0 else 0
        F1 = 2 * (Prec * Rec) / (Prec + Rec) if Prec and Rec else 0
        MCC = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        Jaccard_i = len(T.intersection(P)) / len(T.union(P)) if len(T.union(P)) > 0 else 0

        results[m] = {"Number": len(T), "TP": TP, "FP": FP, "TN": TN, "FN": FN, 
                      "Precision": Prec, "Recall": Rec, "F1": F1, "MCC": MCC, "Jaccard": Jaccard_i}

    logging.info(f"Performance metrics calculated for {'all' if top_genes == -1 else top_genes} genes.")
    return results


def calc_stats(df_input, lfc_thr=1, adj_pval_thr=0.05,
               stats=['Number', "TP", "TN", "FP", "FN", 
                      "Precision", "Recall", "F1", "MCC", 
                      "Jaccard"
                      "r", "ρ", "RMSE", 
                      "MeanDiff", "MaxDiff", "MinDiff"],
               methods=["FedProt", "Fisher", "Stouffer", "REM", "RankProd"],
               column_name="pv_",
               top_genes=-1):

    df_prepared = df_input.copy()
    adjusted_pval_thr = -np.log10(adj_pval_thr)
    all_genes = set(df_prepared.index.values)

    de = df_prepared.loc[(df_prepared["pv_DEqMS"] > adjusted_pval_thr) & (np.abs(df_prepared["lfc_DEqMS"]) >= lfc_thr), :]
    de = de.head(top_genes if top_genes > 0 else df_prepared.shape[0])

    # create a dictionary to store the results
    results = {}

    if any([s in stats for s in ["Number", "TP", "TN", "FP", "FN", "Precision", "Recall", "F1", "MCC", "Jaccard"]]):
        performance_resuts = calculate_performance_metrics(
            df_prepared, de, methods, lfc_thr, adjusted_pval_thr, top_genes, all_genes
        )
        for m in methods:
            if m in results:
                results[m].update(performance_resuts[m])
            else:
                results[m] = performance_resuts[m]

    if "MeanDiff" in stats:
         for m in methods:
            res_value = np.mean(np.abs(df_prepared[column_name + "DEqMS"] - df_prepared[column_name + m]))

            if m in results:
                results[m].update({"MeanDiff": res_value})
            else:
                results[m] = {"MeanDiff": res_value}
    
    if "MaxDiff" in stats:
         for m in methods:
            res_value = np.max(np.abs(df_prepared[column_name + "DEqMS"] - df_prepared[column_name + m]))

            if m in results:
                results[m].update({"MaxDiff": res_value})
            else:
                results[m] = {"MaxDiff": res_value}
    
    if "MinDiff" in stats:
        for m in methods:
            res_value = np.min(np.abs(df_prepared[column_name + "DEqMS"] - df_prepared[column_name + m]))

            if m in results:
                results[m].update({"MinDiff": res_value})
            else:
                results[m] = {"MinDiff": res_value}

    if any([s in stats for s in ["RMSE", "NRMSE"]]):
        rmse_results = calculate_rmse(df_prepared, methods, top_genes, column_name)
        
        for m in methods:
            if m in results:
                results[m].update(rmse_results[m])
            else:
                results[m] = rmse_results[m]

    if "r" in stats or "ρ" in stats:
        correlation_results = calculate_correlations(df_prepared, methods, top_genes, column_name)
        for m in methods:
            if m in results:
                results[m].update({k: correlation_results[k][m] for k in correlation_results})
            else:
                results[m] = {k: correlation_results[k][m] for k in correlation_results}

    results_df = pd.DataFrame.from_dict(results).T
    return results_df.loc[:, stats]


def plt_results(dfs, methods=["FedProt","Fisher","Stouffer","REM","RankProd"],
                color_dict={"Methods": {"FedProt": "blue", "Fisher": "green", "Stouffer": "orange", "REM": "purple", "RankProd": "brown"}},
                what="pv_", 
                text="", dotsize=1,
                datasets=["Balanced", "Imbalanced"],
                add_table=True, sharey=True, sharex=True,
                comparsions=["pyr/glu", "pyr/glu"],
                use_RMSE=False,
                figsize=(11,4.5), after_comma=3,
                set_lims=None,
                titles=None):
    """
    Function to plot results based on different datasets and methods.
    """
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    fig, axes = plt.subplots(1, len(datasets), figsize=figsize, sharey=sharey)
    # plt.subplots_adjust(bottom=0.2)

    if what == "pv_":
        max_p_val = np.max(np.abs(dfs[datasets[0]]['pv_DEqMS']))
        suptitle = "$-log_{10}(adj.p-val.)$" if max_p_val > 1.1 else "adj.p-val."
        logging.info(f"Plotting corrs using p-vals - {'log-transformed' if max_p_val > 1.1  else 'not log-transformed'}.")
    elif what == "lfc_":
        suptitle = "Log2FC"
        logging.info(f"Plotting corrs using logFC values.")

    for i, dataset in enumerate(datasets):
        df = dfs[dataset].filter([f'{what}DEqMS']+[what+method for method in methods])
        axes[i].set_title(titles[i] if titles else dataset, fontsize=16)
        axes[i].set_xlabel(f'{suptitle} {comparsions[i]}, \nDEqMS',fontsize=12)
        axes[i].set_ylabel(f'{suptitle} {comparsions[i]}, \nother methods',fontsize=12)

        mins = []
        maxs = []

        for method in methods:
            x = df[what+"DEqMS"].values
            y = df[what+method].values
            if method == "FedProt":
                # use triangle marker for FedProt
                axes[i].scatter(x, y, s=1.5, color=color_dict["Methods"][method], alpha=0.6, edgecolors=color_dict["Methods"][method], 
                                marker='^', label=method if i == 0 else "")
            else:
                axes[i].scatter(x, y, s=dotsize, color=color_dict["Methods"][method], alpha=0.6, edgecolors=color_dict["Methods"][method],
                                label=method if i == 0 else "")
            
            mins.append(np.min(y))
            maxs.append(np.max(y))

        y_min, y_max = min(mins) - max(maxs) * 0.01, max(maxs) + max(maxs) * 0.02
        x_min, x_max = np.min(df[what+"DEqMS"]) - np.max(df[what+"DEqMS"]) * 0.01, np.max(df[what+"DEqMS"]) + np.max(df[what+"DEqMS"]) * 0.01

        axes[i].set_xlim(x_min, x_max)
        axes[i].set_ylim(y_min, y_max)

        # Adjust limits if necessary
        if set_lims:
            if len(set_lims[i]) == 2:
                axes[i].set_xlim(set_lims[i][0], set_lims[i][1])
                y_min, y_max = set_lims[i][0], set_lims[i][1]
                x_min, x_max = set_lims[i][0], set_lims[i][1]
                axes[i].set_ylim(set_lims[i][0], set_lims[i][1])

        # Display table conditionally
        if add_table:
            display_table(axes[i], df, methods, color_dict, what, use_RMSE, after_comma)

        # Plot identity line
        axes[i].plot([x_min, x_max], [x_min, x_max], color="gray", ls="--", lw=0.2)

    
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.1), 
               title="\nMethods", fontsize="large", markerscale=5, frameon=False, title_fontsize="large", ncol=len(methods))

    if text:
        plt.figtext(0.5, 0.01, text, ha="center", fontsize=12)

    plt.tight_layout()
    plt.show()


def display_table(ax, df, methods, color_dict, what, use_RMSE, after_comma):
    """
    Function to display a table of statistics on the given axis.
    """
    data = {}
    colLabels = ["r", "ρ"] + (['NRMSE*'] if use_RMSE else [])
    for method in methods:
        corr = df[[what+"DEqMS", what+method]].corr().iloc[0, 1]
        rho = df[[what+"DEqMS", what+method]].corr(method='spearman').iloc[0, 1]
        row = [round(corr, after_comma), round(rho, after_comma)]
        if use_RMSE:
            rmse = np.sqrt(np.mean((df[what+"DEqMS"].values - df[what+method].values)**2))
            nrmse = rmse / (df[what+"DEqMS"].max() - df[what+"DEqMS"].min())
            row.append(round(nrmse, after_comma))
        data[method] = row

    cell_colors = [[color_dict["Methods"][method] for _ in colLabels] for method in methods]
    the_table = ax.table(cellText=list(data.values()), colLabels=colLabels, rowLabels=methods, cellLoc='center', loc='top', cellColours=cell_colors)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(8)
    the_table.scale(1, 1.5)


def calc_stats_TOP(
    df,
    stats=["TP","TN","FP","FN","Precision","Recall","F1", "Jaccard"],
    methods=["FedProt","Fisher","Stouffer","REM","RankProd"],
    top_genes=-1):

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    max_p_val = np.max(np.abs(df['pv_DEqMS']))
    if max_p_val < 1.1:
        logging.error("Adj.p-values are not log-transformed. Please, log-transform them before calculating performance metrics.")


    results={}
    all_genes = set(df.index.values)

    if top_genes<=0:
        top_genes = df.shape[0]
    de = df.sort_values(by="pv_DEqMS", ascending = False)
    de = de.head(top_genes)

    T = set(de.index.values)
    F = all_genes.difference(T)
    if len(set(stats).intersection(set(["TP", "TN", "FP", "FN", "Precision", "Recall", "F1", "Jaccard"]))) > 0:
        for m in methods:
            de2 = df.sort_values(by="pv_" + m, ascending=False)
            de2 = de2.head(top_genes)
            P = set(de2.index.values)
            N = all_genes.difference(P)

            TP = len(T.intersection(P))
            FP = len(F.intersection(P))
            TN = len(F.intersection(N))
            FN = len(T.intersection(N))

            # Calculate Precision, Recall, and F1
            Prec = TP / (TP + FP) if (TP + FP) > 0 else 0
            Rec = TP / (TP + FN) if (TP + FN) > 0 else 0
            F1 = 2 * (Prec * Rec) / (Prec + Rec) if Prec and Rec else 0

            # Calculate Jaccard Similarity
            Jaccard = len(T.intersection(P)) / len(T.union(P)) if len(T.union(P)) > 0 else 0
            results[m] = {"TP": TP, "FP": FP, "TN": TN, "FN": FN, "Precision": Prec, "Recall": Rec, "F1": F1, "Jaccard": Jaccard}
    
    if len(results.keys())>0:
        results = pd.DataFrame.from_dict(results).T
    if type(results)==dict:
        results = pd.DataFrame.from_dict(results)
    return results.loc[:,stats]



def plot_stats_for_topN(dfs,
                        datasets=["Balance", "Mild Imbalance", "Strong Imbalance"],
                        metrics=["F1", "Jaccard"],
                        methods=["FedProt", "Fisher", "Stouffer", "REM", "RankProd"],
                        color_dict={"Methods": {"FedProt": "blue", "Fisher": "green", "Stouffer": "orange", "REM": "purple", "RankProd": "brown"}},
                        min_n_genes=10, max_n_genes=1000, step=10,
                        text="",
                        figfile="", suptitle="", sharey=False,
                        titles=None):

    """Calculated and plots statisctics for top N genes ordered by p-value.
    Top genes are chosen based on a sliding threshold, starting from 'min_n_genes' and moving to 'max_n_genes' with 'step'."""
    
    cmap = LinearSegmentedColormap.from_list("from_dict", list(color_dict["Methods"].values()))
    fig, all_axes = plt.subplots(len(metrics), len(datasets), figsize=(13,4*len(metrics)), sharey=sharey)
    all_stats = {}
    min_ylim = {}
    max_ylim = {}

    for k in range(len(metrics)):
        metric = metrics[k]
        all_stats[metric] = {}
        if len(metrics) == 1:
            axes = all_axes
        else:
            axes = all_axes[k]
        for i in range(len(datasets)):
            ds = datasets[i]
            df = dfs[ds]
            df = df.sort_values(by="pv_DEqMS", ascending=False)
            stats = {}
            if isinstance(max_n_genes, list):
                i_max = max_n_genes[i]
                top_n_genes = np.arange(min_n_genes, i_max, step)
            else:
                i_max = max_n_genes
                top_n_genes = np.arange(min_n_genes, i_max, step)

            for j in range(len(top_n_genes)):  #
                confusion_matrix = calc_stats_TOP(df, 
                                                  stats=[metric],
                                                  methods=methods, 
                                                  top_genes=top_n_genes[j])
                stats[top_n_genes[j]] = confusion_matrix[metric]
            stats = pd.DataFrame.from_dict(stats)
            # print(stats.T)
            stats.T.plot(ax=axes[i], cmap=cmap)
            min_ylim[metric] = min(min_ylim.get(metric, 0), stats.values.min())
            max_ylim[metric] = max(max_ylim.get(metric, 0), stats.values.max())

            # axes[i].set_yscale('log')
            if k == len(metrics) - 1:
                tmp = axes[i].set_xlabel("number of top-ranked proteins", fontsize=14)
            if i == 0:
                tmp = axes[i].set_ylabel(f"{metric}", fontsize=14)
                if text:
                    tmp = axes[0].text(-0.15 * i_max, np.max(stats.values) * 1.0, text, fontsize=24)
            if i > 0 or k != len(metrics) - 1:
                axes[i].get_legend().remove()
            if k == 0:
                if titles:
                    tmp = axes[i].set_title(titles[i], fontsize=20)
                else:
                    tmp = axes[i].set_title(ds, fontsize=20)
            all_stats[metric][ds] = stats

    for k in range(len(metrics)):
        metric = metrics[k]
        for i in range(len(datasets)):
            if len(metrics) > 1:
                all_axes[i].set_ylim(min_ylim[metric] - 0.1 * (max_ylim[metric] - min_ylim[metric]),
                                     max_ylim[metric] + 0.1 * (max_ylim[metric] - min_ylim[metric]))              

    for ax_row in np.atleast_1d(all_axes):
        for ax in np.atleast_1d(ax_row):  # Make sure it works for both 1D and 2D arrays of axes
            # ax.plot([min_n_genes, max_n_genes], [1, 1], linestyle='--', color='grey', linewidth=1)
            ax.set_yticks(np.arange(0, 1.1, 0.1))

    if suptitle:
        fig.suptitle(suptitle, fontsize=24)
    if figfile:
        fig.savefig(figfile)
    return all_stats


def plot_with_confidence(jaccard_dfs, methods, color_dict, sharey=True,
                        num_top_genes=range(5, 700, 5),
                        figfile="", figsize=(13, 4),
                        titles=None):
    fig, axes = plt.subplots(1, len(jaccard_dfs), figsize=figsize, sharey=sharey)
    datasets = list(jaccard_dfs.keys())

    # Convert column names to integers to plot as numeric x-axis
    for k, df in jaccard_dfs.items():
        new_columns = [int(col) for col in df.columns]
        df.columns = new_columns
        df.sort_index(axis=1, inplace=True)  # Ensures columns are in numeric order for plotting

    for i, ds in enumerate(datasets):
        df = jaccard_dfs[ds]
        for method in methods:
            method_data = df.loc[method].mean()  # mean across rows for each method
            std_dev = df.loc[method].std()       # standard deviation across rows for each method

            # Getting the numerical range for x-axis from column names
            num_top_genes = df.columns
            mean_scores = method_data
            std_deviation = std_dev
            
            axes[i].plot(num_top_genes, mean_scores, label=method, color=color_dict["Methods"][method])
            axes[i].fill_between(num_top_genes, mean_scores - std_deviation, mean_scores + std_deviation, color=color_dict["Methods"][method], alpha=0.1)
        
        axes[i].set_title(f"Simulated {titles[i].lower()}", fontsize=16)
        axes[i].set_xlabel("Number of top-ranked proteins", fontsize=14)
        axes[i].set_yticks(np.arange(0, 1.1, 0.1))
        if i == 0:
            axes[i].set_ylabel("Jaccard similarity", fontsize=14)
            axes[i].legend(title="Method")

        # add y ticks for the second plot (because they are shared and was removed by sharey=True)

    if figfile:
        fig.savefig(figfile)

    plt.show()



def calc_diffs_for_plotting(log_dfs, methods, what):
    plot_data = []
    for dataset_name, dataset in log_dfs.items():
        for method in methods:
            for q in [0.25, 0.75]:
                low_expr = dataset[dataset['AvgExpr_DEqMS'] < dataset['AvgExpr_DEqMS'].quantile(q)]
                high_expr = dataset[dataset['AvgExpr_DEqMS'] >= dataset['AvgExpr_DEqMS'].quantile(q)]

                diff_high = abs(high_expr[f'{what}{method}'] - high_expr[f'{what}DEqMS'])
                for value in diff_high:
                    plot_data.append((dataset_name, method, 'High', q, value))

                diff_low = abs(low_expr[f'{what}{method}'] - low_expr[f'{what}DEqMS'])
                for value in diff_low:
                    plot_data.append((dataset_name, method, 'Low', q, value))

    return plot_data

def get_significance_level(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''  # not significant


def plot_exp_diffs(log_dfs, what="pv_", methods=["FedProt", "Fisher", "Stouffer", "REM", "RankProd"], figsize=(17, 6), figfile=None):
    plot_data = calc_diffs_for_plotting(log_dfs, methods, what)
    plot_df = pd.DataFrame(plot_data, columns=['Dataset', 'Method', 'Expr_Level', 'Quantile', 'Diffs'])
    datasets = plot_df['Dataset'].unique()

    fig, axes = plt.subplots(len(datasets), len(methods), figsize=figsize, sharey=False)
    p_values = []
    mean_diffs = []

    for row, dataset in enumerate(datasets):
        for col, method in enumerate(methods):
            ax = axes[row, col]
            method_data = plot_df[(plot_df['Method'] == method) & (plot_df['Dataset'] == dataset)]
            high_diff = method_data[method_data['Expr_Level'] == 'High']['Diffs'].values
            low_diff = method_data[method_data['Expr_Level'] == 'Low']['Diffs'].values
            
            stat, p_value = mannwhitneyu(high_diff, low_diff,    alternative='two-sided')

            # Calculate effect size (r)
            N = len(high_diff) + len(low_diff)
            Z = (stat - (len(high_diff)*len(low_diff)/2)) / np.sqrt(len(high_diff)*len(low_diff)*(N + 1)/12)
            r = Z / np.sqrt(N)

            p_values.append((dataset, method, p_value))

            mean_high_diff = np.mean(high_diff)
            mean_low_diff = np.mean(low_diff)
            median_high_diff = np.median(high_diff)
            median_low_diff = np.median(low_diff)
            significance_level = get_significance_level(p_value)
            mean_diffs.append((significance_level, dataset, method, mean_high_diff, mean_low_diff, median_high_diff, median_low_diff, r))

            ax.boxplot([high_diff, low_diff], labels=['High Expr\n(25% quantile)', 'Low Expr\n(75% quantile)'])
            if row == len(datasets) - 1:
                ax.set_xlabel('Expression Level (DEqMS-based)')
            if col == 0:
                if what == 'lfc_':
                    ax.set_ylabel('Abs Log2FC Diff\n{}'.format(dataset))
                elif what == 'pv_':
                    ax.set_ylabel('Abs -log10(adj.p-value) Diff\n{}'.format(dataset))
            ax.set_title(f"{method}{significance_level}")

    plt.tight_layout()
    if figfile:
        plt.savefig(figfile)
    plt.show()

    return p_values, mean_diffs


def plot_ma_plots(log_dfs, what="lfc_", methods=["FedProt", "Fisher", "REM"], lfc_thr=0.5, adj_pval_thr=0.05, 
                  figsize=(15, 5), figfile=None):
    num_datasets = len(log_dfs)
    num_methods = len(methods)
    fig, axes = plt.subplots(num_datasets, num_methods, figsize=figsize)
    
    if num_datasets == 1 and num_methods == 1:
        axes = np.array([[axes]])
    elif num_datasets == 1:
        axes = np.array([axes])
    elif num_methods == 1:
        axes = np.array([[ax] for ax in axes])
    
    for i, (dataset_name, dataset) in enumerate(log_dfs.items()):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            x = dataset['AvgExpr_DEqMS']
            y = dataset[f'{what}{method}']
            pval = dataset[f'pv_{method}']
            color = np.where((pval > -np.log10(adj_pval_thr[i])) & (abs(y) > lfc_thr[i]), 'red', 'blue')
            
            ax.scatter(x, y, c=color, s=1, alpha=0.5)
            ax.set_xlabel('log2 Average Protein Intensity')
            if what == 'lfc_':
                ax.set_ylabel(f'Log2FC, {method}')
            elif what == 'pv_':
                ax.set_ylabel(f'-log10(p-value), {method}')
            ax.set_title(f'MA Plot, {dataset_name}')
            ax.axhline(0, color='grey', lw=0.5)
    
    plt.tight_layout()
    if figfile:
        fig.savefig(figfile)
    plt.show()