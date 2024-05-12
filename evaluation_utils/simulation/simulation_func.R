library(invgamma)
library(RobNorm)
library(tidyverse)

################################################################################################
# use RobNorm func, but with different IG parameters
sim.dat.fn = function(row.frac, col.frac, mu.up, mu.down, n, m, nu.fix=TRUE) {
	mu.00 = rnorm(n, 0, 2) 	# --- changes here
    # mu.00 = rnorm(n, 0, 1)
	var.00 = rinvgamma(n, 2, 3) # --- changes here
    # var.00 = rinvgamma(n, 5, 2)

    X.0 = matrix(rnorm(n*m, 
					outer(mu.00, rep(1, m)),
					sqrt(outer(var.00, rep(1, m)))), 
				n, m) # the null matrix
    
	S = matrix(0, n, m) 
	if (row.frac*col.frac > 0) {
	   	   bk.nm = round(n*row.frac * m*col.frac)
		   a = rbinom(1, bk.nm, 0.8)
		   S[ 1:round(n*row.frac), 1:(m*col.frac)] = sample(c(rep(mu.up, a), rep(0, bk.nm-a)), bk.nm) # the shifted mean of the signal mx
		   a = rbinom(1, bk.nm, 0.8)
		   S[ (n-round(n*row.frac)+1):n, (m-m*col.frac+1):m] = sample(c(rep(mu.down, a), rep(0, bk.nm-a)), bk.nm) # the signal mx of shifted mean 		  
	}

	X = X.0 + S 
	
	rownames(X) = paste("prt", 1:nrow(X), sep=".")
	colnames(X) = paste("s", 1:ncol(X), sep=".")
    return(list(dat=X, raw=X.0))
}


select_proportions <- function(mode_version){
    if(mode_version == "balanced"){
        percent_batch1 <- 0.6
        percent_batch2 <- 0.6
        percent_batch3 <- 0.6
    } else if (mode_version == 'mild_imbalanced') {
        percent_batch1 <- 0.4
        percent_batch2 <- 0.5
        percent_batch3 <- 0.66
    } else if(mode_version == "imbalanced"){
        percent_batch1 <- 0.2
        percent_batch2 <- 0.5
        percent_batch3 <- 0.7
    }

    return(list(percent_batch1, percent_batch2, percent_batch3))

}


generate_data <- function(
    col_frac_A, col_frac_B,
    frac_1, frac_7,
    mu_1 = 1.25, mu_4 = 1.25,
    batch_info = NULL, nu.fix = TRUE,
    mode_version = "balanced",
    m = 600
    ){

    # add A condition
    data_mu1 <- sim.dat.fn(row.frac=frac_1, col.frac=col_frac_A, mu.up=mu_1, mu.down=0, n=2500, m=m, nu.fix=TRUE)
    # add B condition
    data_mu2 <- sim.dat.fn(row.frac=frac_1, col.frac=col_frac_B, mu.up=0, mu.down=mu_1, n=2500, m=m, nu.fix=TRUE)

    # % of confounder in batches
    b_proportions <- select_proportions(mode_version)
    # rearrange the data7
    batch1_size <- round(b_proportions[[1]] * length(batch_info[batch_info$batch == "batch1" & batch_info$condition == "B",]$file))
    batch2_size <- round(b_proportions[[2]] * length(batch_info[batch_info$batch == "batch2" & batch_info$condition == "B",]$file))
    batch3_size <- round(b_proportions[[3]] * length(batch_info[batch_info$batch == "batch3" & batch_info$condition == "B",]$file))
    need_to_generate <- (batch1_size + batch2_size + batch3_size) / m
    # add confounder
    data_mu7 <- sim.dat.fn(row.frac=frac_7, col.frac=need_to_generate, mu.up=0, mu.down=-1*mu_4, n=1000, m=m,  nu.fix=TRUE)
    cat(ncol(data_mu7$dat) * need_to_generate, "\n")

    new_names <- c(
        batch_info[batch_info$batch == "batch1" & batch_info$condition == "B",]$file[1:batch1_size],
        batch_info[batch_info$batch == "batch2" & batch_info$condition == "B",]$file[1:batch2_size],
        batch_info[batch_info$batch == "batch3" & batch_info$condition == "B",]$file[1:batch3_size])

    all_other_names <- c(colnames(data_mu7$dat)[1:length(batch_info[batch_info$condition == "A",]$file)], 
        setdiff(batch_info[batch_info$condition == "B",]$file, new_names))
    
    colnames(data_mu7$dat) <- c(all_other_names, new_names)
    data_mu7$dat <- data_mu7$dat[, batch_info$file]

    # Combine all the data
    data_allmu <- rbind(
        data_mu1$dat[c(1:(frac_1*2500)),],
        data_mu2$dat[c((2500+1-frac_1*2500):2500),],
        data_mu7$dat[c((1000+1-frac_7*1000):1000),],

        data_mu1$dat[c((frac_1*2500+1):2500),],
        data_mu2$dat[c(1:(2500-frac_1*2500)),],
        data_mu7$dat[c(1:(1000-frac_7*1000)),]
        )

    return(data_allmu)
}


add_batch_effect <- function(result_two, batch_info){

    # Assuming 'result_two' is your data matrix and 'batch_info' is a vector indicating the batch for each column
    n_batches <- length(unique(batch_info$batch))
    n_cols <- ncol(result_two)
    n_proteins <- nrow(result_two)

    # For Additive Effects
    additive_params <- data.frame(mean = c(0, 0.2, -0.2), sd = c(1, 0.5, 1.5), row.names = levels(batch_info$batch)) 
    # For Multiplicative Effects
    multiplicative_params <- data.frame(shape = c(3, 2.5, 4), scale = c(2, 1, 0.5), row.names = levels(batch_info$batch)) 

    # Step 3: Sample Based on Selected Parameters
    # Generate additive effects for each batch, one per protein
    additive_effects <- mapply(function(mean, sd) rnorm(n_proteins, mean, sd), additive_params$mean, additive_params$sd)
    colnames(additive_effects) <- levels(batch_info$batch)
    # Generate multiplicative effects for each batch, one per protein
    multiplicative_effects <- mapply(function(shape, scale) rinvgamma(n_proteins, shape, scale), multiplicative_params$shape, multiplicative_params$scale)
    colnames(multiplicative_effects) <- levels(batch_info$batch)

    # Create matrices for applying effects to samples
    additive_effects_matrix <- matrix(nrow = n_proteins, ncol = n_cols)
    multiplicative_effects_matrix <- matrix(nrow = n_proteins, ncol = n_cols)

    # Apply the batch-specific effects to each sample
    for (sample in 1:n_cols) {
        batch <- batch_info$batch[sample]
        additive_effects_matrix[, sample] <- additive_effects[, batch]
        multiplicative_effects_matrix[, sample] <- multiplicative_effects[, batch]
    }

    # Noise
    noise_effect <- matrix(rnorm(n_cols * n_proteins, mean = 0, sd = 1), nrow = n_proteins, ncol = n_cols)

    # Apply batch effects
    data_with_batch_effects <- result_two + additive_effects_matrix + multiplicative_effects_matrix * noise_effect

    return(data_with_batch_effects)
}


simulateMissingValues <- function(df, alpha, beta) {
  # Ensure matrixStats is available
  if (!require(matrixStats)) {
    install.packages("matrixStats")
    library(matrixStats)
  }
  df <- as.matrix(df)

  # Number of total values
  N <- nrow(df) * ncol(df)

  # Step 1: Creating a threshold matrix T
  alpha_quantile <- quantile(as.vector(df), probs = alpha, na.rm = TRUE)
  T <- matrix(rnorm(n = N, mean = alpha_quantile, sd = 0.3), nrow = nrow(df), ncol = ncol(df))

  # Step 2: Creating a probability matrix P
  P <- matrix(rbinom(n = N, size = 1, prob = beta), nrow = nrow(df), ncol = ncol(df))

  # Step 3: Vectorized operation to simulate MNAR
  mnar_mask <- (df < T) & (P == 1)
  df[mnar_mask] <- NA

  # Simulate MAR: Randomly replace additional values without touching MNAR values
  num_MAR <- round(N * alpha * (1 - beta))
  available_indices <- which(is.na(df) == FALSE, arr.ind = TRUE)
  selected_indices <- available_indices[sample(nrow(available_indices), num_MAR), ]
  df[selected_indices] <- NA

  return(df)
}


create_plots <- function(pg_matrix, metadata, name, plot_file_prefix){
        # plots
    plot_pca <- pca_plot(pg_matrix, metadata, 
        title=paste0(name, " PCA"), 
        quantitative_col_name='file', 
        col_col='condition', 
        shape_col='batch')
    plot_boxplot <- boxplot_pg(pg_matrix, metadata, 
        title=paste0(name, " Boxplot"), 
        color_col='condition', 
        quantitativeColumnName='file')
    plot_density <- plotIntensityDensityByPool(pg_matrix, metadata, 
        title=paste0(name, " Density"), 
        poolColumnName='condition', 
        quantitativeColumnName='file')
    plot_heatmap <- heatmap_nocor_plot(
        pg_matrix, metadata, name,
        condition="condition", lab="batch")

    layout <- (plot_density | plot_pca) /
            (plot_boxplot | plot_heatmap)
            
    if(plot_file_prefix == ""){
        return(layout)
    }
    # save plot
    ggsave(
        file = paste0(plot_file_prefix, "_plots.svg"), 
        plot = layout, width = 11, height = 10)
}