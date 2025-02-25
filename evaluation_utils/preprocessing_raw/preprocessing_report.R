library(tidyverse)

#' Preprocess Protein or Peptide Data from MaxQuant Output
#'
#' This function reads a tab-separated file from a specified path and performs preprocessing on the data.
#' It is capable of handling both protein and peptide data from MaxQuant output. The function removes
#' decoy matches, matches to contaminants, and, if applicable, entries only identified by modification sites.
#'
#' @param path The file path of the data to be preprocessed.
#' @param metadata Metadata associated with the data which should include necessary column names or indices.
#' @param data_type Type of the data to process: 'protein' or 'peptide'. This parameter determines
#'                  the subset of columns and processing steps to be applied.
#' @param do_filter Logical flag indicating whether to filter out decoy, contaminant, and modification site-only entries.
#'
#' @return A data frame containing the preprocessed data with selected columns as specified in the metadata.
#'         The function also prints the count of processed entries.
#'
#' @examples
#' # Preprocess protein data
#' preprocess_data_mxout("/path/to/protein_data.txt", protein_metadata, "protein")
#'
#' @export
preprocess_data_mxout <- function(path, metadata, data_type, do_filter=TRUE){
    data_report <- read.table(
        path,
        header=TRUE,
        sep="\t",
        stringsAsFactors=FALSE)
    
    # General preprocessing
    if(do_filter){
        cat("Filtering out decoy, contaminant, and modification site-only entries...\n")
        cat(paste0('Initial data count: ', nrow(data_report)), '\n')
        data_report <- data_report[!data_report$Reverse=="+",]
        data_report <- data_report[!data_report$Potential.contaminant=="+",]
        if(!all(is.na(data_report$Only.identified.by.site))) {
            data_report <- data_report[!data_report$Only.identified.by.site=="+",]
        }
        cat(paste0('Filtered data count: ', nrow(data_report)), '\n')
    }
    
    # Data type specific handling
    if(data_type == "protein") {
        selected_columns <- c("Majority.protein.IDs", "Gene.names", rownames(metadata))
        col_names <- c("Majority.protein.IDs", "Gene.names", metadata$Quantitative.column.name)
        counts_columns <- c("Majority.protein.IDs", "Gene.names", "Peptide.IDs", "Razor...unique.peptides", "Peptide.is.razor")
    } else if(data_type == "peptide") {
        # if Proteins is empty - take the value from Leading.razor.protein
        data_report[data_report$Proteins == '', 'Proteins'] <- data_report$Leading.razor.protein[data_report$Proteins == '']
        selected_columns <- c("Sequence", "Proteins", "Gene.names", rownames(metadata), 'Reverse', 'Potential.contaminant')
        col_names <- selected_columns
    }

    # Select and rename columns
    processed_data <- data_report[,selected_columns] %>% setNames(col_names)
    
    cat(paste0('Processed data count: ', nrow(processed_data)), '\n')
    if(data_type == "protein") {
        pg_counts <- data_report[,counts_columns]
        cat(paste0('Counts data count: ', nrow(pg_counts)), '\n')
        return(list(processed_data, pg_counts))
    } else {
        return(processed_data)
    }
}


#' Preprocess Spectronaut Data
#' 
#' This function reads a tab-separated file from a specified path and performs preprocessing on the data.
#' 
#' @param path The file path of the data to be preprocessed.
#' @param name The name of the dataset to be processed.
#' 
#' 
#' @return A data frame containing the preprocessed data with selected columns as specified in the metadata.
#' 
preprocess_spectronaut <- function(path, name, use_filter = TRUE){
    report_data <- read.delim(
        path,
        header=TRUE,
        sep="\t",
        stringsAsFactors=FALSE)
    
    df_long <- report_data %>%
        # pivot the columns that start with "X." which contain the quantity and count values
        pivot_longer(
            cols = starts_with("X."),
            # names_pattern: 
            #   X\\.[0-9]+\\.\\.(.+)\\.PG\\.(.+)
            names_to = c("file", "measure"),
            names_pattern = "X\\.[0-9]+\\.\\.(.+)\\.PG\\.(.+)",
            values_to = "value"
        ) %>%
        # Now pivot wider to have separate columns for the two measures
        pivot_wider(
            names_from = measure,
            values_from = value
        ) %>%
        # Create File.Name and Run (both get the same value extracted above)
        mutate(File.Name = file,
               Run = file)

    # Rename columns for clarity: change the original protein group and gene columns,
    # and rename the measures to the names you want.
    if(name == "lab_E"){
        df_long <- df_long %>%
            rename(Protein.Group = PG.ProteinAccessions,
                   Genes = PG.Genes,
                   PG.Quantity = Quantity,
                   PG.Count = NrOfStrippedSequencesUsedForQuantification)
    } else {
        df_long <- df_long %>%
            rename(Protein.Group = PG.ProteinGroups,
                Genes = PG.Genes,
                PG.Quantity = Quantity,
                PG.Count = NrOfStrippedSequencesUsedForQuantification
                # PG.Count = RunEvidenceCount
                )
    }
    
    # if column PG.Qvalue is present, filter the data based on it
    if("PG.Qvalue" %in% colnames(df_long) && use_filter){
        df_long <- df_long %>%
            filter(PG.Qvalue <= 0.01)
    }
    df_long <- df_long %>%
        # remove .raw from the file name and Run
        mutate(File.Name = gsub("\\.raw", "", File.Name),
               Run = gsub("\\.raw", "", Run)) %>%
        # Finally, select the columns in the desired order.
        select(File.Name, Run, Protein.Group, Genes, PG.Quantity, PG.Count)

    # outlier samples based on name
    filter_sample_out <- switch(name,
        'lab_A' = c("Ref8537_QC1_20230414_2", 'Ref8537_QC2_20230414_2', 'Ref8537_QC3_20230414_2', 'Ref8537_QC4_20230414_2'),
        NULL)
        
    # filter data if filter_sample_out is not NULL
    if (!is.null(filter_sample_out)) {
        df_long <- df_long[!df_long$File.Name %in% filter_sample_out, ]
    }
    return(df_long)
}