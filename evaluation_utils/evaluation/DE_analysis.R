library(DEqMS)
library(tidyverse)


run_DE <- function(result_two, batch_info, labs = NULL){

  if(is.null(labs)){
    design <- model.matrix(~0 + batch_info$condition)
    colnames(design) <- levels(as.factor(batch_info$condition))
    contrasts <- makeContrasts(A - B, levels = colnames(design))
  } else {
    design <- model.matrix(~0 + batch_info$condition + batch_info$batch)
    colnames(design) <- c(levels(batch_info$condition), levels(batch_info$batch)[-1])
    contrasts <- makeContrasts(A - B, levels = colnames(design))
  }
  

    fit <- lmFit(result_two, design)
    fit2 <- contrasts.fit(fit, contrasts)
    fit3 <- eBayes(fit2)
    result <- topTable(fit3, coef=1, adjust="BH", sort.by="none", number=Inf, confint = TRUE)

    n_de <- result %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% nrow()
    return(list(result, n_de))
}

