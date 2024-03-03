# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_MetaVolcanoR.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/balanced/ lab_A lab_E lab_B lab_D lab_C &&
# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_MetaVolcanoR.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/imbalanced/  lab_A lab_E lab_D lab_B lab_C

library("MetaVolcanoR")
#https://bioconductor.org/packages/release/bioc/vignettes/MetaVolcanoR/inst/doc/MetaVolcano.html#random-effect-model-metavolcano

args <- commandArgs(trailingOnly = TRUE)
w_dir <- args[1] #
cohorts <- args[2:length(args)] # all but 1st arguments are cohort names

diffexplist <- list()

for (cohort in cohorts){ #,"Other")){
    fname <- paste0(w_dir,cohort,"_res.tsv")
    res <- read.table(fname, row.names = 1, sep="\t")
    res["Symbol"] <- rownames(res)
    diffexplist[[cohort]] <- res
    
}

################################################################################################################################
meta_degs_comb <- combining_mv(diffexp=diffexplist,
                   pcriteria="sca.P.Value",
                   foldchangecol='logFC', 
                   genenamecol='Symbol',
                   geneidcol=NULL,
                   metafc='Mean',
                   metathr=0.01, 
                   collaps=TRUE,
                   jobname="MetaVolcano",
                   outputfolder=".",
                   draw='HTML')

result <- meta_degs_comb@metaresult
result <- result[order(result$metap),]
write.table(result,paste0(w_dir,"/MA_CM.tsv"),row.names=TRUE,sep="\t", quote = FALSE)
#head(result,3)


meta_degs_rem <- rem_mv(diffexp=diffexplist,
            pcriteria="sca.P.Value",
            foldchangecol='logFC', 
            genenamecol='Symbol',
            geneidcol=NULL,
            collaps=FALSE,
            llcol='CI.L',
            rlcol='CI.R',
            vcol=NULL, 
            cvar=TRUE,
            metathr=0.01,
            jobname="MetaVolcano",
            outputfolder=".", 
            draw='HTML',
            ncores=1)

result <- meta_degs_rem@metaresult
result  <- result[order(result$randomP),]
write.table(result,paste0(w_dir,"/MA_REM.tsv"),row.names=TRUE,sep="\t", quote = FALSE)
#head(result,3)

