import pandas as pd
import numpy as np 
from statsmodels.stats.multitest import multipletests

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
                rankprod_name="/MA_RankProd.tsv"):
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    df = {}
    
    # Central analysis
    DEqMS = pd.read_csv(workdir+deqms_name, sep="\t", index_col=0)
    df["pv_DEqMS"] = DEqMS["sca.adj.pval"]
    # df["pv_DEqMS"] = DEqMS["sca.P.Value"]
    df["lfc_DEqMS"] = DEqMS["logFC"]

    # FedProt
    fedprot = pd.read_csv(workdir+fedprot_name, sep="\t", index_col=0)
    df["pv_FedProt"] = fedprot["sca.adj.pval"]
    # df["pv_FedProt"] = fedprot["sca.P.Value"]
    df["lfc_FedProt"] = fedprot["logFC"]

    # Fisher
    ma_cm = pd.read_csv(workdir+fisher_name, sep="\t")
    ma_cm.index = ma_cm["Symbol"].values
    df["lfc_Fisher"] = ma_cm["metafc"]
    _, adj_pval,_,_ = multipletests(ma_cm["metap"].values, alpha=0.05, method='fdr_bh',
                                    is_sorted=False, returnsorted=False)
    df["pv_Fisher"] = pd.Series(adj_pval,index=ma_cm["metap"].index)

    # REM
    ma_rem = pd.read_csv(workdir+rem_name, sep="\t")
    ma_rem.index = ma_rem["Symbol"].values
    df["lfc_REM"] = ma_rem["randomSummary"]
    _, adj_pval, _, _ = multipletests(ma_rem["randomP"].values, alpha=0.05, method='fdr_bh',
                                      is_sorted=False, returnsorted=False)
    df["pv_REM"] = pd.Series(adj_pval,index=ma_rem["randomP"].index)

    ### Stoufer 
    stoufer  = pd.read_csv(workdir+stouffer_name, sep="\t", index_col=0)
    df["pv_Stouffer"] = stoufer["FDR"]
    df["lfc_Stouffer"] = df["lfc_Fisher"]  # take logFC from MetaVolcanoR

    ### RankProd
    rankprod  = pd.read_csv(workdir+rankprod_name, sep="\t", index_col=0)
    rankprod["FDR"] = rankprod.loc[:,["down_reg.FDR","up_reg.FDR"]].min(axis=1)
    df["pv_RankProd"] = rankprod["FDR"]
    df["lfc_RankProd"] = rankprod["avgL2FC"] 
    
    df = pd.DataFrame.from_dict(df)
    df = df.dropna(axis=0)

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
        nrmse = rmse / np.var(x)
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

        results[m] = {"Number": len(T), "TP": TP, "FP": FP, "TN": TN, "FN": FN, 
                      "Precision": Prec, "Recall": Rec, "F1": F1, "MCC": MCC}

    logging.info(f"Performance metrics calculated for {'all' if top_genes == -1 else top_genes} genes.")
    return results


def calc_stats(df_input, lfc_thr=1, adj_pval_thr=0.05,
               stats=['Number', "TP", "TN", "FP", "FN", "Precision", "Recall", "F1", "MCC", 
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

    if any([s in stats for s in ["Number", "TP", "TN", "FP", "FN", "Precision", "Recall", "F1", "MCC"]]):
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
                add_table=True
    ):
    
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    fig, axes = plt.subplots(1, len(datasets), figsize=(10,4.5), sharey=False)
    i=0
    se = 0 
    results = {}
    max_xlim = 0
    min_xlim = 0

    if what == "pv_":
        max_p_val = np.max(np.abs(dfs[datasets[0]]['pv_DEqMS']))
        suptitle = "$-log_{10}(sca.adj.p-val.)$" if max_p_val > 1.1 else "sca.adj.p-val."
        logging.info(f"Plotting corrs using p-vals - {'log-transformed' if max_p_val > 1.1  else 'not log-transformed'}.")
    else:
        suptitle = "logFC"
        logging.info(f"Plotting corrs using logFC values.")


    for k in datasets:
        df = dfs[k].filter([f'{what}DEqMS']+[what+i for i in methods])
        axes[i].set_title(k, fontsize=16)
        rmse = {}
        
        for j in range(len(methods)):
            method = methods[j]
            col = color_dict["Methods"][methods[j]]
            x = df[what+"DEqMS"].values
            y = df[what+method].values
            rmse[method] = np.sqrt(np.sum((x-y)**2)/len(x))
            axes[i].scatter(x = np.abs(x), y= np.abs(y),s=dotsize, color=col, alpha=0.5)
        
        axes[i].set_xlabel(f'{suptitle} (pyr/glu), \nDEqMS',fontsize=12)
        axes[i].set_ylabel(f'{suptitle} (pyr/glu), \nother methods',fontsize=12)
        axes[i].plot([np.min(np.abs(df.values)), np.max(np.abs(df.values))+5], 
                     [np.min(np.abs(df.values)), np.max(np.abs(df.values))+5],
                   color = "red",ls="--",lw=0.5)
     
        # Calculate correlations
        corrs = df[[what+"DEqMS"]+[what+m for m in methods]].corr().loc[[what+"DEqMS"],]
        corrs.rename(lambda x: x.replace(what,""), axis="columns",inplace = True)
        corrs = corrs.T.to_dict()[what+'DEqMS']
        rank_corrs = df[[what+"DEqMS"]+[what+m for m in methods]].corr(method="spearman").loc[[what+"DEqMS"],]
        rank_corrs.rename(lambda x: x.replace(what,""), axis="columns",inplace = True)
        rank_corrs = rank_corrs.T.to_dict()[what+'DEqMS']

        if add_table:
            # Prepare data for table
            data = {}
            colLabels = ["r", "ρ"]

            for j, method in enumerate(methods):
                data[method] = [f"{round(corrs[method],3)}", f"{round(rank_corrs[method],3)}"]
                
            # Create table for each axes
            the_table = table(
                axes[i], 
                cellText=list(data.values()),
                colLabels=colLabels,
                rowLabels=list(data.keys()),
                cellLoc = 'center', rowLoc = 'right',
                bbox=[0.6, 0.66, 0.4, 0.34],
                cellColours=[['white'] * len(colLabels) for col in color_dict["Methods"]]
            )
            
            # Set font size for the entire table
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(10)  # Change this value as required

            for j, label in enumerate(data.keys()):
                cell = the_table.get_celld()[(j + 1, -1)]  # get the row label cell
                cell.set_facecolor(color_dict["Methods"][label])  # set the row label background color
                if label == "FedProt" or label == "RankProd":
                    cell.get_text().set_color('white')  # set the row label text color 
                else:
                    cell.get_text().set_color('black')

        else:
            axes[i].legend(methods, loc='upper right', title='Methods', fontsize=10, markerscale=8, markerfirst=False)
            

        i += 1
        results[(k,"r")] = corrs
        results[(k,"ρ")] = rank_corrs
        results[(k,"RMSE")] = pd.Series(rmse)

        max_xlim_method = np.max(np.abs(df[what+"DEqMS"].values))
        if what == "pv_":
            max_xlim_method = max_xlim_method + max_xlim_method*0.2
        else:
            max_xlim_method = max_xlim_method + max_xlim_method*0.05
        if max_xlim < max_xlim_method:
            max_xlim = max_xlim_method

        min_xlin_method = np.min(np.abs(df[what+"DEqMS"].values))
        min_xlin_method = min_xlin_method - max_xlim_method*0.01

        if min_xlim > min_xlin_method:
            min_xlim = min_xlin_method

    
    axes[0].set_xlim([min_xlim, max_xlim])
    axes[1].set_xlim([min_xlim, max_xlim])
    
    results = pd.DataFrame.from_dict(results)
    
    if text:
        tmp = axes[0].text(-0.2*np.max(df.values), np.max(df.values), text, fontsize=24)
        
    plt.tight_layout()
    return results.loc[methods,]


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
                        figfile="", suptitle="", sharey=False):

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
            top_n_genes = np.arange(min_n_genes, max_n_genes, step)
            for j in range(len(top_n_genes)):  #
                confusion_matrix = calc_stats_TOP(df, 
                                                  stats=[metric],
                                                  methods=methods, 
                                                  top_genes=top_n_genes[j])
                stats[top_n_genes[j]] = confusion_matrix[metric]
            stats = pd.DataFrame.from_dict(stats)
            stats.T.plot(ax=axes[i], cmap=cmap)
            min_ylim[metric] = min(min_ylim.get(metric, 0), stats.values.min())
            max_ylim[metric] = max(max_ylim.get(metric, 0), stats.values.max())

            # axes[i].set_yscale('log')
            if k == len(metrics) - 1:
                tmp = axes[i].set_xlabel("number of top-ranked proteins", fontsize=14)
            if i == 0:
                tmp = axes[i].set_ylabel(f"{metric}", fontsize=14)
                if text:
                    tmp = axes[0].text(-0.15 * max_n_genes, np.max(stats.values) * 1.0, text, fontsize=24)
            if i > 0 or k != len(metrics) - 1:
                axes[i].get_legend().remove()
            if k == 0:
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
