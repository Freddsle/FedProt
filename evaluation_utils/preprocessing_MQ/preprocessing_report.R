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
        data_report <- data_report[!data_report$Reverse=="+",]
        data_report <- data_report[!data_report$Potential.contaminant=="+",]
        if(!all(is.na(data_report$Only.identified.by.site))) {
            data_report <- data_report[!data_report$Only.identified.by.site=="+",]
        }
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
