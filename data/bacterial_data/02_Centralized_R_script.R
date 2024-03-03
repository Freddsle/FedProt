# Analysis only

library(tidyverse)
library(gtools)
library(limma)

library(VennDiagram)
library(ggrepel)





################################################################################
################################################################################
################################################################################
# Balanced #
################################################################################
################################################################################
################################################################################

# Prepare data
# batch info
batch_info <- read.csv("/home/yuliya/repos/cosybio/FedProt/bacterial_data/balanced/bath_info_all.tsv", check.names = FALSE, sep="\t") %>%
  column_to_rownames('rowname') %>%
  mutate(lab = factor(lab), condition = factor(condition))

# PG matrix and counts
labs_list = c('lab_A', 'lab_E', 'lab_B', 'lab_D', 'lab_C')  
pg_matrix <- NULL
prec_counts_table <- NULL

for (name in labs_list) {
  file_name_prefix <- paste0('/home/yuliya/repos/cosybio/FedProt/bacterial_data/balanced/', name, "/")

  if(is.null(pg_matrix)){
    pg_matrix <- read.csv(paste0(file_name_prefix, 'protein_groups_matrix.tsv'), check.names = FALSE, sep="\t") 
    prec_counts_table <- read.csv(paste0(file_name_prefix, 'protein_counts.tsv'), check.names = FALSE, sep="\t") 
  } else {
    pg_matrix <- full_join(pg_matrix, 
                       read.csv(paste0(file_name_prefix, 'protein_groups_matrix.tsv'), check.names = FALSE, sep="\t"),
                       by = "rowname")
    prec_counts_table <- full_join(prec_counts_table, 
                               read.csv(paste0(file_name_prefix, 'protein_counts.tsv'), check.names = FALSE, sep="\t"),
                       by = "rowname")
  }
}


rownames(pg_matrix) <- pg_matrix$rowname
pg_matrix$rowname <- NULL
tmp_rownames <- read.csv("/home/yuliya/repos/cosybio/FedProt/evaluation/balanced/results/DPE_fedprot_app.csv", sep="\t")$X
pg_matrix_filter <- pg_matrix[tmp_rownames, rownames(batch_info)]
pg_matrix_filter <- log2(pg_matrix_filter + 1)
dim(pg_matrix_filter)

## counts
rownames(prec_counts_table) <- prec_counts_table$rowname
prec_counts_table$rowname <- NULL
prec_counts_table <- prec_counts_table %>% transmute(count = pmin(!!!., na.rm = TRUE))
prec_counts_table$Protein.Group <- rownames(prec_counts_table)
dim(prec_counts_table)

# write central data
write_tsv(pg_matrix %>% rownames_to_column(), 
          '/home/yuliya/repos/cosybio/FedProt/bacterial_data/balanced/centralized/protein_groups_matrix.tsv')
write_tsv(prec_counts_table %>% rownames_to_column(),
          '/home/yuliya/repos/cosybio/FedProt/bacterial_data/balanced/centralized/protein_counts.tsv')

## Analysis
design <- model.matrix(~0 + batch_info$condition + batch_info$lab)
colnames(design) <- c(levels(batch_info$condition), levels(batch_info$lab)[-1])
contrasts <- makeContrasts(Glu - Pyr, levels = colnames(design))

fit <- lmFit(pg_matrix_filter, design)
fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2)

fit3$count <- prec_counts_table[rownames(fit3$coefficients), "count"] + 1
if(is.na(min(fit3$count)) | min(fit3$count) == 0){
  print("No DEqMS results")
}

library(DEqMS)
fit4 <- spectraCounteBayes(fit3)
DEqMS_results <- outputResult(fit4, coef_col = 1)


DEqMS_results %>% rownames_to_column("PG") %>% write_tsv("/home/yuliya/repos/cosybio/FedProt/evaluation/balanced/results/DPE_deqms_central.tsv")


################################################################################
################################################################################
################################################################################
# Imbalanced #
################################################################################
################################################################################
################################################################################

# Prepare data
## batch info
batch_info <- read.csv("/home/yuliya/repos/cosybio/FedProt/bacterial_data/imbalanced/bath_info_all.tsv", check.names = FALSE, sep="\t") %>%
  column_to_rownames('rowname') %>%
  mutate(lab = factor(lab), condition = factor(condition))
rownames(batch_info) <- batch_info$file

# PG matrix and counts
labs_list = c('lab_A', 'lab_E', 'lab_D', 'lab_B', 'lab_C') 
pg_matrix <- NULL
prec_counts_table <- NULL

for (name in labs_list) {
  file_name_prefix <- paste0('/home/yuliya/repos/cosybio/FedProt/bacterial_data/imbalanced/', name, '/')
  if(is.null(pg_matrix)){
    pg_matrix <- read.csv(paste0(file_name_prefix, 'protein_groups_matrix.tsv'), check.names = FALSE, sep="\t")
    prec_counts_table <- read.csv(paste0(file_name_prefix, 'protein_counts.tsv'), check.names = FALSE, sep="\t")
  } else {
    pg_matrix <- full_join(pg_matrix, 
                       read.csv(paste0(file_name_prefix, 'protein_groups_matrix.tsv'), check.names = FALSE, sep="\t"),
                       by = "rowname") 
    prec_counts_table <- full_join(prec_counts_table, 
                               read.csv(paste0(file_name_prefix, 'protein_counts.tsv'), check.names = FALSE, sep="\t"),
                       by = "rowname")
  }
}

rownames(pg_matrix) <- pg_matrix$rowname
pg_matrix$rowname <- NULL
tmp_rownames <- read.csv("/home/yuliya/repos/cosybio/FedProt/evaluation/imbalanced/results/DPE_fedprot_app.csv", sep="\t")$X
pg_matrix_filter <- pg_matrix[tmp_rownames, rownames(batch_info)] 
pg_matrix_filter <- log2(pg_matrix_filter + 1)
dim(pg_matrix_filter)

## counts
rownames(prec_counts_table) <- prec_counts_table$rowname
prec_counts_table$rowname <- NULL
prec_counts_table <- prec_counts_table %>% transmute(count = pmin(!!!., na.rm = TRUE))
prec_counts_table$Protein.Group <- rownames(prec_counts_table)
dim(prec_counts_table)

# write central data
write_tsv(pg_matrix %>% rownames_to_column(), 
          '/home/yuliya/repos/cosybio/FedProt/bacterial_data/imbalanced/centralized/protein_groups_matrix.tsv')
write_tsv(prec_counts_table %>% rownames_to_column(),
          '/home/yuliya/repos/cosybio/FedProt/bacterial_data/imbalanced/centralized/protein_counts.tsv')

## Analysis
design <- model.matrix(~0 + batch_info$condition + batch_info$lab)
colnames(design) <- c(levels(batch_info$condition), levels(batch_info$lab)[-1])
contrasts <- makeContrasts(Glu - Pyr, levels = colnames(design))

fit <- lmFit(pg_matrix_filter, design)
fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2)

fit3$count <- prec_counts_table[rownames(fit3$coefficients), "count"] + 1
if(is.na(min(fit3$count)) | min(fit3$count) == 0){
  print("No DEqMS results")
}

library(DEqMS)
fit4 <- spectraCounteBayes(fit3)
DEqMS_results <- outputResult(fit4, coef_col = 1)


DEqMS_results %>% rownames_to_column("PG") %>% write_tsv("/home/yuliya/repos/cosybio/FedProt/evaluation/imbalanced/results/DPE_deqms_central.tsv")


