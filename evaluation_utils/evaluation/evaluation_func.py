import pandas as pd
import numpy as np 
from statsmodels.stats.multitest import multipletests

import logging


def read_results(workdir, 
                 fedprot_name="/DPE.csv", 
                 deqms_name="/Central_res_irs_on_median.tsv"):
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
    ma_cm = pd.read_csv(workdir+"/MA_CM.tsv", sep="\t")
    ma_cm.index = ma_cm["Symbol"].values
    df["lfc_Fisher"] = ma_cm["metafc"]
    _, adj_pval,_,_ = multipletests(ma_cm["metap"].values, alpha=0.05, method='fdr_bh',
                                    is_sorted=False, returnsorted=False)
    df["pv_Fisher"] = pd.Series(adj_pval,index=ma_cm["metap"].index)

    # REM
    ma_rem = pd.read_csv(workdir+"/MA_REM.tsv", sep="\t")
    ma_rem.index = ma_rem["Symbol"].values
    df["lfc_REM"] = ma_rem["randomSummary"]
    _, adj_pval, _, _ = multipletests(ma_rem["randomP"].values, alpha=0.05, method='fdr_bh',
                                      is_sorted=False, returnsorted=False)
    df["pv_REM"] = pd.Series(adj_pval,index=ma_rem["randomP"].index)

    ### Stoufer 
    stoufer  = pd.read_csv(workdir+"/MA_Stouffer.tsv", sep="\t", index_col=0)
    df["pv_Stouffer"] = stoufer["FDR"]
    df["lfc_Stouffer"] = df["lfc_Fisher"]  # take logFC from MetaVolcanoR

    ### RankProd
    rankprod  = pd.read_csv(workdir+"/MA_RankProd.tsv", sep="\t", index_col=0)
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
        logging.info(f"Calculating corrs. Using p-vals - {'log-transformed' if max_p_val > 50 else 'not log-transformed'}.")

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
        logging.info(f"Calculating RMSE. Using p-vals - {'log-transformed' if max_p_val > 50 else 'not log-transformed'}.")

    rmse_results = {}
    df_sorted = df.sort_values(by="pv_DEqMS", ascending=False).head(top_genes)
    for m in methods:        
        x = df_sorted[column_name + "DEqMS"].values
        y = df_sorted[column_name + m].values
        rmse = np.sqrt(np.sum((x - y) ** 2) / len(x))
        rmse_results[m] = {"RMSE": rmse}

    logging.info(f"RMSE computed for {'all' if top_genes == -1 else top_genes} genes from {column_name} columns.")
    return rmse_results


def calculate_performance_metrics(
    df, de, 
    methods, 
    lfc_thr, adj_pval_thr, 
    top_genes, all_genes
    ):

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    max_p_val = np.max(np.abs(de['pv_DEqMS']))
    logging.info(f"Calculating performance. Using p-vals - {'log-transformed' if max_p_val > 50 else 'not log-transformed'}.")

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
               stats=['Number', "TP", "TN", "FP", "FN", 
                      "Precision", "Recall", "F1", "MCC", 
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

    if "RMSE" in stats:
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











# FROM HERE ON, THE FUNCTIONS WERE NOT MODIFIED!

from matplotlib.table import table

def plt_results(dfs, methods=["FedProt","Fisher","Stouffer","REM","RankProd"],
                colors=["D44500","2E5EAA","FFFB0A","47A025","010B13"], 
                what="pv_", suptitle="$-log_{10}(adj.p-values)$", text="",dotsize=1,
                datasets=["Balanced", "Imbalanced"]):
    
    fig, axes = plt.subplots(1, 2, figsize=(10,4.5), sharey=False)
    i=0
    se = 0
    results = {}
    max_xlim = 0
    min_xlim = 0

    for k in datasets:
        df = dfs[k].filter([f'{what}DEqMS']+[what+i for i in methods])
        axes[i].set_title(k,fontsize=16)
        rmse = {}
        
        for j in range(len(methods)):
            method = methods[j]
            col = color_dict["Methods"][methods[j]]
            x = df[what+"DEqMS"].values
            y = df[what+method].values
            rmse[method] = np.sqrt(np.sum((x-y)**2)/len(x))
            axes[i].scatter(x = np.abs(x), y= np.abs(y),s=dotsize, color=col, alpha=0.5)
        
        axes[i].set_xlabel('-log10 adj.p-values (pyr/glu), DEqMS',fontsize=12)
        axes[i].set_ylabel('-log10 adj.p-values (pyr/glu), other methods',fontsize=12)
        axes[i].plot([np.min(np.abs(df.values)), np.max(np.abs(df.values))+5], 
                     [np.min(np.abs(df.values)), np.max(np.abs(df.values))+5],
                   color = "red",ls="--",lw=0.5)
     
        corrs = df[[what+"DEqMS"]+[what+m for m in methods]].corr().loc[[what+"DEqMS"],]
        corrs.rename(lambda x: x.replace(what,""), axis="columns",inplace = True)
        corrs = corrs.T.to_dict()[what+'DEqMS']
        rank_corrs = df[[what+"DEqMS"]+[what+m for m in methods]].corr(method="spearman").loc[[what+"DEqMS"],]
        rank_corrs.rename(lambda x: x.replace(what,""), axis="columns",inplace = True)
        rank_corrs = rank_corrs.T.to_dict()[what+'DEqMS']

        # Prepare data for table
        data = {}
        for j, method in enumerate(methods):
            if method == "FedProt":
                if rmse[method] < 1:
                    data[method] = [f"{round(corrs[method],3)}", f"{round(rank_corrs[method],3)}", f"{rmse[method]:.0e}"]
                else:
                    data[method] = [f"{round(corrs[method],3)}", f"{round(rank_corrs[method],3)}", f"{round(rmse[method], 2)}"]
            else:
                data[method] = [f"{round(corrs[method],3)}", f"{round(rank_corrs[method],3)}", f"{round(rmse[method],2)}"]

        # Create table for each axes
        colLabels = ["r", "ρ", "RMSE"]
        the_table = table(axes[i], cellText=list(data.values()),
                        colLabels=colLabels,
                        rowLabels=list(data.keys()),
                        cellLoc = 'center', rowLoc = 'right',
                        bbox=[0.5, 0.66, 0.5, 0.34],  
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






def calculate_stats_for_topN(dfs, datasets=["Balanced"],
                             metrics=["F1"], methods=["FedProt", "Fisher", "Stouffer", "REM", "RankProd"],
                             min_n_genes=10, max_n_genes=1000, step=10,
                             lfc_thr=1.0, adj_pval_thr=0.05):
    """
    Calculates statistics for top N genes ordered by p-value.
    Top genes are chosen based on a sliding threshold, starting from 'min_n_genes'
    and moving to 'max_n_genes' with 'step'.
    """
    all_stats = {}
    for metric in metrics:
        all_stats[metric] = {}
        for ds in datasets:
            df = dfs[ds]
            df = df.sort_values(by="pv_DEqMS", ascending=False)
            top_n_genes = np.arange(min_n_genes, max_n_genes, step)
            stats = {}
            for n_genes in top_n_genes:
                confusion_matrix = calc_stats(df, lfc_thr, adj_pval_thr,
                                                  stats=[metric], methods=methods, top_genes=n_genes)
                stats[n_genes] = confusion_matrix[metric]
            all_stats[metric][ds] = pd.DataFrame.from_dict(stats)

    return all_stats
