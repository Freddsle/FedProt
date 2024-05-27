library(tidyverse)
library(data.table)
library(foreach)

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

filter_per_center <- function(intensities, metadata, quantitative_column_name, centers, center_column_name, min_number=2) {
  cat('Filtering by', center_column_name, ' - ', min_number, ' not-NA per ', center_column_name, '\n')
  cat('\tBefore filtering:', dim(intensities), "\n")
  
  # Initialize a list to store the sample indices for each center
  center_samples <- list()
  
  # Loop through each center and extract relevant sample names
  for (center in centers) {
    center_samples[[center]] <- metadata[metadata[[center_column_name]] == center, ][[quantitative_column_name]]
  }
  # Determine rows with at least 2 non-NA values across each center's samples
  conditions <- sapply(center_samples, function(samples) {
    rowSums(!is.na(intensities[, samples, drop = FALSE])) >= min_number
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
    rowSums(is.na(intensities[, samples, drop = FALSE])) / length(samples) <= min_f & rowSums(!is.na(intensities[, samples, drop = FALSE])) >= 2
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
  norm_dt <- as.data.frame(norm_dt)
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


#######################################################
get_samples <- function(metadata, pool, column_name = "file") {
  filtered_data <- metadata[metadata$Pool == pool, ] %>%
    dplyr::pull(all_of(column_name)) %>%
    as.character()
  return(filtered_data)
}

get_intensities <- function(intensities, metadata, pool) {
  samples <- get_samples(metadata, pool)
  intensities <- intensities[, samples]
  return(intensities)
}

aggregate_data <- function(data, method) {
  if (method == "average") {
    return(rowMeans(data, na.rm = TRUE))
  } else {
    return(rowSums(data, na.rm = TRUE))
  }
}

calculate_geometric_mean <- function(irs) {
  geometric_mean <- apply(irs, 1, function(x) {
    non_zero_values <- x[x > 0]  # Exclude zeros from the computation
    non_zero_values <- non_zero_values[!is.na(non_zero_values)]
    exp(mean(log(non_zero_values), na.rm = TRUE))  # Compute geometric mean of non-zero values
  })
  
  return(geometric_mean)
}

# Function to process data for a single center
irsNorm_in_silico_single_center <- function(data, metadata, pool_col = "Pool",
                                            column_name = "file",
                                            center = NULL,
                                            aggregation_method = "average", add_refs=TRUE) {
  # Initialize IRS table
  irs <- tibble()

  # Process each pool to create IRS references
  pools <- unique(metadata[[pool_col]])
  for (pool in pools) {
    samples <- get_samples(metadata, pool, column_name)
    intensities <- data[, samples]
    cat("Shape of intensities for pool ", pool, ":", dim(intensities), "\n")
    processed_data <- aggregate_data(intensities, aggregation_method)
    
    if (ncol(irs) == 0) {
        irs <- tibble(!!pool := processed_data)  
    } else {
        irs[[pool]] <- processed_data
    }
  }

  # Compute geometric mean across all IRS references
  irs_average <- calculate_geometric_mean(irs)

  # Calculate scaling factors and normalize the data
  if (add_refs) {
    corrected_data <- cbind(data, irs)
  } else {
    corrected_data <- cbind(data)
  }
  for (pool in pools) {
    scaling_factor <- irs_average / irs[[pool]]
    scaling_factor[is.na(scaling_factor) | is.infinite(scaling_factor)] <- 1 
    if (add_refs){
      pool_samples <- c(get_samples(metadata, pool, column_name), pool)
    } else {
      pool_samples <- get_samples(metadata, pool, column_name)
    }
    corrected_data[, pool_samples] <- corrected_data[, pool_samples] * scaling_factor
  }

  if(add_refs) {
    # update metadata with IRS column names
    template <- setNames(as.list(rep("NA", ncol(metadata))), names(metadata))

    new_metadata_rows <- do.call(rbind, lapply(pools, function(pool) {
      condition_value <- "in_silico"
      new_row <- template
      new_row$Pool <- pool
      new_row$file <- pool
      new_row$condition <- condition_value
      new_row$lab <- center
      return(new_row)
    }))
    metadata <- metadata %>% mutate_all(as.character)
    metadata <- rbind(metadata, new_metadata_rows)
    # transform all columns to character
    metadata <- metadata %>% mutate_all(as.character) %>% as.data.frame()
    rownames(metadata) <- metadata$file

    return(list(corrected_data = corrected_data, metadata = metadata))
    
  } else {
    return(list(corrected_data = corrected_data)
    )
  }

}