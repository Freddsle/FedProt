import pandas as pd
import numpy as np 
from statsmodels.stats.multitest import multipletests


def read_results(workdir, i):
    df = {}
    
    rlimma = pd.read_csv(f'{workdir}/limma_central_{i}.tsv', sep="\t", index_col=0)
    df["pv_Rlimma"] = -np.log10(rlimma["adj.P.Val"])
    df["lfc_Rlimma"] = rlimma["logFC"]

    # Fisher
    ma_cm = pd.read_csv(workdir+"/MA_CM.tsv", sep="\t")
    ma_cm.index = ma_cm["Symbol"].values
    df["lfc_Fisher"] = ma_cm["metafc"]
    _, adj_pval,_,_ = multipletests(ma_cm["metap"].values, alpha=0.05, method='fdr_bh',
                                           is_sorted=False, returnsorted=False)
    df["pv_Fisher"] = -np.log10(pd.Series(adj_pval,index=ma_cm["metap"].index))

    # REM
    ma_rem = pd.read_csv(workdir+"/MA_REM.tsv", sep="\t")
    ma_rem.index = ma_rem["Symbol"].values
    df["lfc_REM"] = ma_rem["randomSummary"]
    _, adj_pval,_,_ = multipletests(ma_rem["randomP"].values, alpha=0.05, method='fdr_bh',
                                           is_sorted=False, returnsorted=False)
    df["pv_REM"] = -np.log10(pd.Series(adj_pval,index=ma_rem["randomP"].index))

    # fedprot
    fedprot = pd.read_csv(f"{workdir}/results.FedProt_{i}.tsv", sep="\t", index_col=0)
    df["pv_FedProt"] = -np.log10(fedprot["adj.P.Val"])
    df["lfc_FedProt"] = fedprot["logFC"]

    ### Stoufer 
    stoufer  = pd.read_csv(workdir+"/MA_Stouffer.tsv", sep="\t", index_col=0)
    df["pv_Stouffer"] = -np.log10(stoufer["FDR"])
    df["lfc_Stouffer"] = df["lfc_Fisher"]  # take logFC from MetaVolcanoR

    ### RankProd
    rankprod  = pd.read_csv(workdir+"/MA_RankProd.tsv", sep="\t", index_col=0)
    rankprod["FDR"] = rankprod.loc[:,["down_reg.FDR","up_reg.FDR"]].min(axis=1)
    df["pv_RankProd"] = -np.log10(rankprod["FDR"])
    df["lfc_RankProd"] = rankprod["avgL2FC"] 
    
    df = pd.DataFrame.from_dict(df)
    df = df.dropna(axis=0)
    return df


def calculate_performance_metrics(df, de, methods, lfc_thr, adj_pval_thr, top_genes, all_genes):
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

    return results


def calculate_rmse(df, methods, top_genes):
    rmse_results = {}
    df_sorted = df.sort_values(by="pv_Rlimma", ascending=False).head(top_genes)
    for m in methods:        
        x = df_sorted["pv_Rlimma"].values
        y = df_sorted["pv_" + m].values
        rmse = np.sqrt(np.sum((x - y) ** 2) / len(x))
        rmse_results[m] = {"RMSE": rmse}
    return rmse_results


def calculate_correlations(df, methods, top_genes):
    correlations = {}
    df_sorted = df.sort_values(by="pv_Rlimma", ascending=False).head(top_genes)
    pearson_corr = df_sorted[["pv_" + "Rlimma"] + ["pv_" + m for m in methods]].corr().loc[["pv_" + "Rlimma"],]
    spearman_corr = df_sorted[["pv_" + "Rlimma"] + ["pv_" + m for m in methods]].corr(method="spearman").loc[["pv_" + "Rlimma"],]

    pearson_corr.rename(lambda x: x.replace("pv_", ""), axis="columns", inplace=True)
    spearman_corr.rename(lambda x: x.replace("pv_", ""), axis="columns", inplace=True)

    correlations['r'] = pearson_corr.T['pv_Rlimma']
    correlations['ρ'] = spearman_corr.T['pv_Rlimma']
    return correlations


def calc_stats(df_input, lfc_thr=1, adj_pval_thr=0.05,
               stats=['Number', "TP", "TN", "FP", "FN", "Precision", "Recall", "F1", "MCC", "r", "ρ", "RMSE", "MeanLFCdif", "MeanPValdif"],
               methods=["FedProt", "Fisher", "Stouffer", "REM", "RankProd"],
               top_genes=-1):

    df_prepared = df_input.copy()
    adjusted_pval_thr = -np.log10(adj_pval_thr)
    all_genes = set(df_prepared.index.values)

    de = df_prepared.loc[(df_prepared["pv_Rlimma"] > adjusted_pval_thr) & (np.abs(df_prepared["lfc_Rlimma"]) >= lfc_thr), :]
    de = de.head(top_genes if top_genes > 0 else df_prepared.shape[0])

    results = calculate_performance_metrics(df_prepared, de, methods, lfc_thr, adjusted_pval_thr, top_genes, all_genes)

    if "MeanLFCdif" in stats:
        for m in methods:
            if m in results:
                results[m]["MeanLFCdif"] = np.mean(np.abs(df_prepared["lfc_Rlimma"] - df_prepared["lfc_" + m]))
    
    if "MeanPValdif" in stats:
        for m in methods:
            if m in results:
                results[m]["MeanPValdif"] = np.mean(np.abs(df_prepared["pv_Rlimma"] - df_prepared["pv_" + m]))


    if "RMSE" in stats:
        rmse_results = calculate_rmse(df_prepared, methods, top_genes)
        for m in methods:
            if m in results:
                results[m].update(rmse_results[m])
            else:
                results[m] = rmse_results[m]

    if "r" in stats or "ρ" in stats:
        correlation_results = calculate_correlations(df_prepared, methods, top_genes)
        for m in methods:
            if m in results:
                results[m].update({k: correlation_results[k][m] for k in correlation_results})
            else:
                results[m] = {k: correlation_results[k][m] for k in correlation_results}

    results_df = pd.DataFrame.from_dict(results).T
    return results_df.loc[:, stats]


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
            df = df.sort_values(by="pv_Rlimma", ascending=False)
            top_n_genes = np.arange(min_n_genes, max_n_genes, step)
            stats = {}
            for n_genes in top_n_genes:
                confusion_matrix = calc_stats(df, lfc_thr, adj_pval_thr,
                                                  stats=[metric], methods=methods, top_genes=n_genes)
                stats[n_genes] = confusion_matrix[metric]
            all_stats[metric][ds] = pd.DataFrame.from_dict(stats)

    return all_stats
