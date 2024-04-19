# conda activate deqms_fed
# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_MetaDE.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/balanced/ lab_A lab_E lab_B lab_D lab_C &&
# Rscript /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/run_MetaDE.R /home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/imbalanced/ lab_A lab_E lab_D lab_B lab_C  


args <- commandArgs(trailingOnly = TRUE)
w_dir <- args[1] #
cohorts <- args[2:length(args)] # all but 1st arguments are cohort names

#suppressPackageStartupMessages(library("MetaDE"))

log2FC <- data.frame()
pvals <- data.frame()
pvals_sca <- data.frame()

for (cohort in cohorts){ #,"Other")){
    fname <- paste0(w_dir, cohort, "_res.tsv")
    res <- read.table(fname, row.names = 1, sep="\t", header=TRUE)

    lfc <- res["logFC"]
    # pv <- res["P.Value"]
    pv <- res["sca.P.Value"]  
    #res$rank <- c(1:dim(res)[[1]])
    #rank <- res["rank"]
    if (dim(log2FC)[[2]] == 0) {
        log2FC <- lfc
        pvals <- pv
    } else {
        lfc <- lfc[rownames(log2FC),]
        log2FC <-  cbind(log2FC, lfc)
        pv <- pv[rownames(pvals),] 
        pvals <-  cbind(pvals, pv)
    }
}

colnames(log2FC) <- cohorts
colnames(pvals) <- cohorts
log2FC <- as.matrix(log2FC)
pvals <- as.matrix(pvals)

ind_res <- list("log2FC" = log2FC, "p" = pvals)


#-----------------------------------------------------------------------------#
#  MetaDE functions                                                           #
#-----------------------------------------------------------------------------#

# MetaDE functions:

MetaDE.pvalue <-function(x, meta.method, rth=NULL, parametric=TRUE) {
  
  check.parametric(meta.method, parametric)
  K <- ncol(x$p)
  if (parametric) x$bp <- NULL     
  nm <- length(meta.method)
  meta.res <- list(stat=NA, pval=NA, FDR=NA, AW.weight=NA)
  meta.res$stat <- meta.res$pval <-meta.res$FDR <-matrix(NA, nrow(x$p), nm)
  for(i in 1:nm){
    temp <- switch(meta.method[i],
                 maxP={get.maxP(x$p,x$bp)}, minP={get.minP(x$p,x$bp)},
                 Fisher={get.fisher(x$p,x$bp)}, roP={get.roP(x$p,x$bp,rth=rth)},
                 AW={get.AW(x$p)},
                 Fisher.OC={get.fisher.OC(x$p,x$bp)},
                 maxP.OC={get.maxP.OC(x$p,x$bp)},
                 minP.OC={get.minP.OC(x$p,x$bp)},
                 roP.OC={get.roP.OC(x$p,x$bp,rth=rth)},
                 Stouffer={get.Stouff(x$p,x$bp)},
                 Stouffer.OC={get.Stouff.OC(x$p,x$bp)},
                 SR={get.SR(x$p,x$bp)},PR={get.PR(x$p,x$bp)})
    meta.res$stat[,i]<-temp$stat
    meta.res$pval[,i]<-temp$pval
    meta.res$FDR[,i]<-temp$FDR
    if(meta.method[i] == "AW"){
      meta.res$AW.weight <- temp$AW.weight
    }
  }
  colnames(meta.res$stat) <- colnames(meta.res$pval) <- colnames(meta.res$FDR) <- meta.method
  rownames(meta.res$stat) <- rownames(meta.res$pval) <- rownames(meta.res$FDR) <- rownames(x$p)   
  attr(meta.res,"nstudy") <- K
  attr(meta.res,"meta.method") <- meta.method 
  res<-list(meta.analysis=meta.res, ind.p=x$p)	 
  #class(res)<-"MetaDE.pvalue"
  return(res)
}

#-----------------------------------------------------------------------------#
#   check if parametric is ok                                                 #
#-----------------------------------------------------------------------------#
check.parametric <- function(meta.method,parametric) {
  if (parametric == TRUE & sum(meta.method %in% c("SR","PR","rankProd", "Fisher.OC","minMCC")) > 0) {
    stop(paste("There is no parametric result for", meta.method))
  }
}

get.Stouff<-function(p, bp=NULL, miss.tol=0.3){
	k <- ncol(p)
	pval <- stat <- rep(NA, ncol(p))
	if(!is.null(bp)){
		rnum <- which(apply(p, 1, function(x) !any(is.na(x))))
		Ug <- matrix(NA, nrow(p), 1)
		Ug[rnum, 1] <- apply(p[rnum,], 1, function(x)sum(qnorm(x))/sqrt(k))
		Ubg <- matrix(apply(bp, 1, function(x)sum(qnorm(x))/sqrt(k)), nrow(p), nrow(bp)/nrow(p))
		pval[rnum] <- perm.p(Ug[rnum,1], Ubg[rnum,], tail="abs")
		qval <- p.adjust(pval, method="BH")
		res <- list(stat=c(Ug), pval=pval, FDR=qval)
	} else {
		rnum <- which(apply(p, 1, function(x) sum(is.na(x))/k) <=miss.tol)
		pval[rnum] <- apply(p[rnum,], 1, function(x) {
		    ## calculate Stouffer p-value 
        2*pnorm(abs(sum(qnorm(x), na.rm=T) / sqrt(sum(!is.na(x)))), lower.tail=FALSE) 
		})
		stat[rnum] <- apply(p[rnum,], 1, function(x) {
        sum(qnorm(x), na.rm=T) / sqrt(sum(!is.na(x)))
		})
		qval <- p.adjust(pval, method="BH")
		res <- list(stat=stat, pval=pval, FDR=qval)
	}
	  names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}

get.fisher<-function(p, bp=NULL, miss.tol=0.3) {
  k <- ncol(p)
	pval <- stat <- rep(NA, nrow(p))
	if(!is.null(bp)){
		rnum <- which(apply(p, 1, function(x) !any(is.na(x))))
		Ug <- matrix(NA, nrow(p), 1)
		Ug[rnum, 1] <- as.matrix(rowSums(-2*log(p[rnum,])))
		Ubg <- matrix(rowSums(-2*log(bp)), nrow(p), nrow(bp)/nrow(p))
		pval[rnum] <- perm.p(Ug[rnum,1], Ubg[rnum,], tail="high")
		qval <- p.adjust(pval, method="BH")
		res <- list(stat=Ug, pval=pval, FDR=qval)
	} else {
		rnum <- which(apply(p, 1, function(x) sum(is.na(x))/k) <miss.tol)
		pval[rnum] <- apply(p[rnum,], 1, function(x) {
        x <- as.numeric(x)  # Convert to numeric
        ## calculate Fisher p-value 
        pchisq(-2*sum(log(x), na.rm=T), lower.tail=FALSE, 2*sum(!is.na(x)))
		})
		qval <- p.adjust(pval, method="BH")
		stat[rnum] <- apply(p[rnum,], 1, function(x)-2*sum(log(x), na.rm=T))
		res <- list(stat=stat, pval=pval, FDR=qval)
	}
	names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
  return(res)
}

get.roP <- function(p,bp=NULL,rth) {
	k <- ncol(p)
	rnum <- which(apply(p, 1, function(x) !any(is.na(x))))
	pval <- stat <- rep(NA, nrow(p))
	if(!is.null(bp)){
		p <- t(apply(p, 1, sort, na.last = T))
		bp <- t(apply(bp, 1, sort, na.last = T))
		Ug <- matrix(NA, nrow(p), 1)
		Ug[rnum,1] <- p[rnum,rth]	
		Ubg <- matrix(bp[,rth], nrow(p), nrow(bp)/nrow(p))
		pval[rnum] <- perm.p(Ug[rnum,1], Ubg[rnum,], tail="low")
		qval <- p.adjust(pval, method="BH")
		res <- list(stat=Ug, pval=pval, FDR=qval)
	} else {
		pval[rnum] <- apply(p[rnum,], 1, function(x) {
		  ## calculate rOP p-value 
      pbeta(x[order(x)][rth], rth, k - rth + 1)
		}) 
		stat[rnum] <- apply(p[rnum,], 1, function(x)x[order(x)][rth])			
		qval <- p.adjust(pval, method="BH")
		res <- list(stat=stat, pval=pval, FDR=qval)	
	}
	  names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}


#-----------------------------------------------------------------------------#
#   Continue our analysis                                                     #
#-----------------------------------------------------------------------------#

methods <- c("Stouffer", "Fisher", "roP")

for (i in 1:3){ # Stouffer,roP(rth=2), Fisher
    meta.method <- methods[i]
    print(meta.method)
    rth <- NULL
    if (meta.method == "roP"){
        rth <- max(2, length(cohorts) %/% 2)
    }
    meta.res <- MetaDE.pvalue(ind_res, meta.method, rth=rth, parametric=T)
    result <- data.frame(ind.p = meta.res$ind.p,
                         stat = meta.res$meta.analysis$stat,
                         pval = meta.res$meta.analysis$pval,
                         FDR = meta.res$meta.analysis$FDR,
                         weight = meta.res$meta.analysis$AW.weight
                        )
    colnames(result)[seq(length(cohorts) + 1, length(cohorts) + 3)] <- c("stat", "pval", "FDR")
    write.table(result, paste0(w_dir, "/MA_", meta.method, ".tsv"), row.names=TRUE, sep="\t", quote = FALSE, dec = ".")

    rm(meta.res, result)
}

