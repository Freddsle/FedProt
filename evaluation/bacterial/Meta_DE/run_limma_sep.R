# For each lab perform the similar check and save results

# Load libraries
library(tidyverse)
library(stringr)
library(DEqMS)

library(VennDiagram)

############################################################################################################
# functions for plots #
############################################################################################################

plot_limma_results <- function(limma_result, plot_path, title){
  limma_result$log.adj.P.Val = -log10(limma_result$adj.P.Val)
  limma_result$PG <-  rownames(limma_result)
  
  plot <- ggplot(limma_result, aes(x = logFC, y = log.adj.P.Val )) + 
      geom_point(size = 0.5 )+
      theme_bw(base_size = 16) + # change theme
      xlab(expression("log2(A/B)")) + # x-axis label
      ylab(expression("-log10(P-value)")) + # y-axis label
      geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
      geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
      geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
      scale_colour_gradient(low = "black", high = "black", guide = FALSE) +
      ggrepel::geom_text_repel(data = subset(limma_result, abs(logFC) > 1 & log.adj.P.Val > 3),
                      aes(logFC, log.adj.P.Val ,label=PG)) + # add PG label
      ggtitle(title) # add title

    ggsave(plot_path, plot)
}

plot_deqms_results <- function(DEqMS_result, plot_path, title) {
    DEqMS_result$log.sca.pval <- -log10(DEqMS_result$sca.P.Value)
    DEqMS_result$gene <- rownames(DEqMS_result)

    plot <- ggplot(DEqMS_result, aes(x=logFC, y=log.sca.pval )) + 
        geom_point(size=0.5 )+
        theme_bw(base_size = 16) + # change theme
        xlab(expression("log2(A/B)")) + # x-axis label
        ylab(expression(" -log10(P-value)")) + # y-axis label
        geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
        geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
        geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
        scale_colour_gradient(low = "black", high = "black", guide = FALSE)+
        ggrepel::geom_text_repel(data = subset(DEqMS_result, abs(logFC) > 1 & log.sca.pval > 3),
                        aes(logFC, log.sca.pval, label=gene)) + # add PG label
      ggtitle(title) # add title

    ggsave(plot_path, plot)
}

##################
# Customize function from DEqMS for meta analysis
outputResult <-function(fit,coef_col=1){
    results.table = limma::topTable(fit, coef=coef_col, n= Inf, confint=TRUE)
    
    results.table$gene = rownames(results.table)
    results.table$count = fit$count[results.table$gene]
    
    results.table$sca.t = fit$sca.t[results.table$gene,coef_col]
    results.table$P.Value = fit$sca.p[results.table$gene,coef_col]
    results.table$sca.P.Value = fit$sca.p[results.table$gene,coef_col]
    results.table$adj.P.Val = p.adjust(results.table$sca.P.Value, method = "BH")
    results.table$sca.adj.pval = p.adjust(results.table$sca.P.Value, method = "BH")
    results.table = results.table[order(results.table$sca.P.Value), ]
}


############################################################################################################
# Analysis #
############################################################################################################

datasets_list = c('balanced', 'imbalanced')

for(dataset in datasets_list){
    before_vars <- ls()
    # open the file "my_log.txt" for writing
    sink(paste0("/home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/", dataset, "/", dataset, "_log.txt"))
    
    counts_table <- data.frame(Steps = c("Samples", "Raw", "D_Filtered", "DEqMS", "DEqMSLFC"))
    limma_all_PLF <- list()

    labs_list = c('lab_A', 'lab_B', 'lab_C', 'lab_D' , 'lab_E')  

    tmp_rownames <- read.csv(paste0("//home/yuliya/repos/cosybio/FedProt/evaluation/", dataset, "/results/DPE_fedprot_app.csv"), sep="\t")$X

    batch_info_all <- read.csv(paste0("/home/yuliya/repos/cosybio/FedProt/bacterial_data/", dataset, "/bath_info_all.tsv"), check.names = FALSE, sep="\t") %>%
      column_to_rownames('rowname') %>%
      mutate(lab = factor(lab), condition = factor(condition))

    for (name in labs_list) {
        print(name)
        file_name_prefix <- paste0('/home/yuliya/repos/cosybio/FedProt/bacterial_data/', dataset, '/', name, '/')
        batch_info <- batch_info_all[batch_info_all$lab == name, ]

        pg_matrix <- read.csv(paste0(file_name_prefix, 'protein_groups_matrix.tsv'), check.names = FALSE, sep="\t") 
        rownames(pg_matrix) <- pg_matrix$rowname
        pg_matrix$rowname <- NULL
        count_raw <- dim(pg_matrix)[1]
        
        print(dim(pg_matrix))
        pg_matrix <- pg_matrix[tmp_rownames, batch_info$file]
        
        pg_matrix <- pg_matrix[!apply(pg_matrix, 1, function(x) all(is.na(x))), ]
        print(dim(pg_matrix))
        pg_matrix <- log2(pg_matrix + 1)
            
        prec_counts_table <- read.csv(paste0(file_name_prefix, 'protein_counts.tsv'), check.names = FALSE, sep="\t") 
        rownames(prec_counts_table) <- prec_counts_table$rowname

        save_results <- paste0('/home/yuliya/repos/cosybio/FedProt/evaluation/Meta_DE/', dataset, '/', name)

        count_filtered <- dim(pg_matrix)[1]
        samples_filtered <- dim(pg_matrix)[2]

        # create design matrix
        design <- model.matrix(~0 + batch_info$condition)
        colnames(design) <- c(levels(batch_info$condition))
        contrasts <- makeContrasts(Glu - Pyr, levels = colnames(design))

        # run limma
        fit <- lmFit(pg_matrix, design)
        fit2 <- contrasts.fit(fit, contrasts)
        fit3 <- eBayes(fit2)
        fit3$count <- prec_counts_table[rownames(fit3$coefficients), "count"]

        if(is.na(min(fit3$count))){
            print("No DEqMS results - NA")
        } else if (min(fit3$count) == 0) {
            print("No DEqMS results - 0")
        }
        else {
            fit4 <- spectraCounteBayes(fit3)
            DEqMS_results <- outputResult(fit4, coef_col = 1)
            # save results to tsv file
            write.table(DEqMS_results, paste0(save_results, "_res.tsv"), sep = "\t")

            # plot Venn diagrams
            deqms_P_list <- filter(DEqMS_results, sca.adj.pval < 0.05)
            count_deqmsP <- length(rownames(deqms_P_list))

            deqms_PLF_list <- filter(DEqMS_results, sca.adj.pval < 0.05 & abs(logFC) > 1)
            count_deqmsPLF <- length(rownames(deqms_PLF_list))
            # add DEP to the vector for the final VENN diagramm
            limma_all_PLF[[name]] <- rownames(deqms_PLF_list)

            #add counts to table
            counts_table[name] <- c(samples_filtered, count_raw, count_filtered, count_deqmsP, count_deqmsPLF)
    }}

    rownames(counts_table) <- counts_table$Steps
    counts_table$Steps <- NULL
    counts_table$MEAN <- rowMeans(counts_table)
    print(t(counts_table))
    
    # turn off the sink 
    sink()

    # Identify the newly created variables in the loop
    after_vars <- ls()
    new_vars <- setdiff(after_vars, before_vars)
    
    # Remove the new variables
    rm(list = new_vars)
}
