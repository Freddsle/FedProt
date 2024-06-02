# Rscript run_MetaDE_and_RankProd.R working_dir/ cli1 cli2 ...
# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_RankProd.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/balanced/ lab_A lab_E lab_B lab_D lab_C &&
# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_RankProd.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/imbalanced/ lab_A lab_E lab_D lab_B lab_C  


suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))

args <- commandArgs(trailingOnly = TRUE)
w_dir <- args[1] #
cohorts <- args[2:length(args)] # all but 1st arguments are cohort names

lfcs <- NULL
for (cohort in cohorts) {
  fname <- paste0(w_dir, cohort, "_res.tsv")
  res <- read.table(fname, row.names = 1, sep="\t", header = TRUE)
  res$ID <- rownames(res)
  lfc <- res[, c("ID", "logFC")]
  if (is.null(lfcs)) {
    lfcs <- lfc
  } else {
    lfcs <- full_join(lfcs, lfc, by = "ID")
  }
}
rownames(lfcs) <- lfcs$ID
lfcs$ID <- NULL
colnames(lfcs) <- cohorts
lfcs <- as.matrix(lfcs)

classes <- rep(1, length(cohorts))
origins <- rep(1, length(cohorts))#c(1:length(cohorts))

RP_obj_lfc <-  RP.advance(lfcs,classes,origins, logged=TRUE, gene.names=rownames(lfcs), rand=1)

pvals <- RP_obj_lfc$pval
colnames(pvals) <- c("down_reg.pval","up_reg.pval")
fdr <- RP_obj_lfc$pval
colnames(fdr) <- c("down_reg.FDR","up_reg.FDR")
avgL2FC <- RP_obj_lfc$AveFC
colnames(avgL2FC) <- c("avgL2FC")
result <- cbind(pvals,fdr,avgL2FC)
dim(result)
write.table(result,paste0(w_dir,"/MA_RankProd.tsv"),row.names=TRUE,sep="\t", quote = FALSE,dec = ".")
	