library(DEqMS)
library(tidyverse)


# Customize function from DEqMS for meta analysis
outputResult <-function(fit, coef_col=1){
    results.table = limma::topTable(fit, coef=coef_col, n= Inf, confint=TRUE)
    
    results.table$gene = rownames(results.table)
    results.table$count = fit$count[results.table$gene]
    
    results.table$sca.t = fit$sca.t[results.table$gene,coef_col]
    results.table$P.Value = fit$p.value[results.table$gene,coef_col]
    results.table$sca.P.Value = fit$sca.p[results.table$gene,coef_col]
    results.table$adj.P.Val = p.adjust(results.table$P.Value, method = "BH")
    results.table$sca.adj.pval = p.adjust(results.table$sca.P.Value, method = "BH")
    results.table = results.table[order(results.table$sca.P.Value), ]
}

#############################################################
#############################################################
# create design
make_design <- function(batch_info, condition_column, labs = NULL){
    if(is.null(labs)){
    design <- model.matrix(~0 + batch_info[[condition_column]])
    colnames(design) <- levels(as.factor(batch_info[[condition_column]]))
  } else {
    design <- model.matrix(~0 + batch_info[[condition_column]] + batch_info[[labs]])
    colnames(design) <- c(levels(batch_info[[condition_column]]), levels(batch_info[[labs]])[-1])
  }
  return(design)
}

# Function to run DE analysis
run_DE <- function(intensities, counts, design, contrasts){
  fit <- lmFit(intensities, design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit3 <- eBayes(fit2)

  if(is.null(counts)){
    results.table = limma::topTable(fit3, coef=1, n=Inf, confint=TRUE)
    cat("Count information is not available\n")
    cat("Using P.Value and adj.P.Val as sca.P.Value and sca.adj.pval\n")
    # rename
    results.table$sca.adj.pval = results.table$adj.P.Val
    results.table$sca.P.Value = results.table$P.Value
    return(results.table)
  }

  fit3$count <- counts[rownames(fit3$coefficients), "count"] + 1
  if(is.na(min(fit3$count)) | min(fit3$count) == 0){
    print("No DEqMS results")
  }
      
  fit4 <- spectraCounteBayes(fit3)
  result <- outputResult(fit4, coef_col = 1)
  
    return(result)
}

