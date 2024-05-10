library(tidyverse)
library(data.table)

#################################################################
# Filtering functions
#################################################################
filter_na_proteins <- function(dt, meta_data, quantitative_column_name) {
    cat('Filtering out features that have NAs in all columns\n')
    cat('\tBefore filtering:', dim(dt), '\n')
    # Filter out proteins that have NAs in all columns - 2 (only two non NA)
    dt <- dt[!rowSums(is.na(dt[, (meta_data[[quantitative_column_name]])])) == length(meta_data[[quantitative_column_name]]),]
    cat('\tAfter filtering:', dim(dt), '\n')
    return(dt)
}

filter_per_center <- function(intensities, metadata, quantitative_column_name, centers, center_column_name) {
  cat('Filtering by center - two not-NA per center\n')
  cat('\tBefore filtering:', dim(intensities), "\n")
  
  # Initialize a list to store the sample indices for each center
  center_samples <- list()
  
  # Loop through each center and extract relevant sample names
  for (center in centers) {
    center_samples[[center]] <- metadata[metadata[[center_column_name]] == center, ][[quantitative_column_name]]
  }
  # Determine rows with at least 2 non-NA values across each center's samples
  conditions <- sapply(center_samples, function(samples) {
    rowSums(!is.na(intensities[, samples, drop = FALSE])) >= 2
  })
  # Filter intensities where all conditions across centers are met
  filtered_intensities <- intensities[rowSums(conditions) == length(centers), ]

  cat('\tAfter filtering:', dim(filtered_intensities), "\n")
  return(filtered_intensities)
}

filter_by_condition <- function(intensities, metadata, quantitative_column_name, groups, groups_column_name, min_f=0.8) {
  cat('Filtering by condition - two not-NA per condition\n')
  cat('\tBefore filtering:', dim(intensities), "\n")

  # Initialize a list to store sample indices for each group
  condition_samples <- list()

  # Extract sample names for each condition using the specified groups_column_name column
  for (group in groups) {
    condition_samples[[group]] <- metadata[metadata[[groups_column_name]] == group, ][[quantitative_column_name]]
  }

  # Determine rows with at least 2 non-NA values for each group's samples
  conditions <- sapply(condition_samples, function(samples) {
    rowSums(is.na(intensities[, samples, drop = FALSE])) / length(samples) <= min_f
  })
  # Filter intensities where all conditions are met
  filtered_intensities <- intensities[rowSums(conditions) == length(groups), ]

  cat('\tAfter filtering:', dim(filtered_intensities), "\n")
  return(filtered_intensities)
}



#################################################################
# Normalization functions
#################################################################
# from PRONE.R
medianNorm <- function(dt){
  # find median of each sample
  sample_med <- apply(dt, 2, stats::median, na.rm=TRUE) # columns
  # find mean of medians
  mean_med <- mean(sample_med, na.rm=TRUE)
  # divide data by median
  norm_dt <- t(t(dt)/sample_med)
  # multiply data by mean of medians
  norm_dt <- norm_dt * mean_med
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  return(norm_dt)
}


# from PRONE.R
irsNorm <- function(dt, md, batch, refs){
  # get md of reference samples
  refs_md <- md[md$Quantitative.column.name %in% refs,]
# separate data by batch
  dt_list <- lapply(unique(md[["Pool"]]), function(b){
      md_chunk <- md[md[["Pool"]] == b,]
      dt_chunk <- subset(dt, select = md_chunk$Quantitative.column.name)
      return(dt_chunk)
  })
  names(dt_list) <- unique(md[["Pool"]])
  
  # take reference sample intensities
  irs <- subset(dt, select = refs_md$Quantitative.column.name)
  colnames(irs) <- as.character(refs_md[refs_md$Quantitative.column.name %in% refs,][["Pool"]])
  
  # get the geometric average intensity for each protein
  irs <- tibble::as_tibble(irs)
  irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
  # normalize data
  dt_irs_list <- lapply(names(dt_list), function(b){
    # compute scaling factor vectors
    fac <- irs$average / irs[,b]
    # normalize
    dt_irs_chunk <- dt_list[[b]] * fac[,1]
    return(dt_irs_chunk)
  })
  # reconstruct data after irs normalization
  dt_irs <- do.call(cbind, dt_irs_list)
  dt_irs <- data.table::as.data.table(dt_irs)
  dt_irs <- subset(dt_irs, select = colnames(dt))
  
  return(dt_irs)
}