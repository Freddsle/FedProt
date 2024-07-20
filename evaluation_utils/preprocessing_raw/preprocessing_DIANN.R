library(tidyverse)

# read and preprocess bacterial data from DIA-NN reports (MBR on)
read_and_preprocc <- function(path, name) {
    # read data
    data <- diann::diann_load(path)
    data$File.Name <- data$Run

    
    # outlier samples based on name
    filter_sample_out <- switch(name,
        'lab_C' = c("BBM_673_P283_01_VEB_008_R2"),
        'lab_A' = c("Ref8537_QC1_20230414_2", 'Ref8537_QC2_20230414_2', 'Ref8537_QC3_20230414_2', 'Ref8537_QC4_20230414_2'),
        # 'lab_B' = c('Clinspect_E_coli_A_S29_Slot1-19_1_8668'),
        NULL)
        
    # filter data if filter_sample_out is not NULL
    if (!is.null(filter_sample_out)) {
        data <- data[!data$File.Name %in% filter_sample_out, ]
    }

    # further filter data
    data <- data %>% filter(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01)

    if(name == 'lab_B'){
        data <- data %>% mutate(File.Name = ifelse(
                File.Name == 'Clinspect_E_coli_B_66_Slot1-13_1_8647', 'Clinspect_E_coli_B_S66_Slot1-13_1_8647', File.Name)) %>%
            mutate(Run = ifelse(Run == 'Clinspect_E_coli_B_66_Slot1-13_1_8647', 'Clinspect_E_coli_B_S66_Slot1-13_1_8647',Run))
    }
    
    return(data)
}

# function to create expression matrix
extract_pg_matrix <- function(data) {
  # extract protein group matrix
  pg_matrix <- data %>% 
    select(Protein.Group, Run, PG.MaxLFQ) %>% 
    unique() %>% 
    pivot_wider(names_from = Run, values_from = PG.MaxLFQ) %>% 
    column_to_rownames("Protein.Group") %>% 
    # transform 0 to NA
    replace(., . == 0, NA)
  return(pg_matrix)
}


get_pep_counts_table <- function(data, sample_names=NULL) {
    if (is.null(sample_names)) {
        pre_prec.count.table <- data %>%
        filter(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01) %>%
        select(c(Protein.Group, Precursor.Id, Run)) %>%
        unique() 
    } else {
       pre_prec.count.table <- data %>%
        filter(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01) %>%
        filter(Run %in% sample_names) %>%
        select(c(Protein.Group, Precursor.Id, Run)) %>%
        unique()
    }
    
    prec.count.table <- pre_prec.count.table %>%
        select(c(Protein.Group, Precursor.Id)) %>%
        unique()
        
    summ.prec.count.table <- pre_prec.count.table %>%
        group_by(Run, Protein.Group) %>% 
        summarise(count = n_distinct(Precursor.Id)) %>% 
        pivot_wider(names_from = Run, values_from = count) %>%
        # add new column with the minimum number per row, with rm.na = TRUE
        mutate(count = pmin(!!!.[-1], na.rm = TRUE)) %>%
        select(Protein.Group, count) %>%
        as.data.frame(.) 

    rownames(summ.prec.count.table) <- summ.prec.count.table$Protein.Group
    summ.prec.count.table$Protein.Group <- NULL

    return(list(summ.prec.count.table, prec.count.table))
}