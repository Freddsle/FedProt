
library(tidyverse)
library(gridExtra)
library(patchwork)
library(grid)

pca_plot <- function(
    df, 
    batch_info, 
    title, 
    path = "", 
    quantitative_col_name = "Quantitative.column.name", 
    col_col = "Group", 
    shape_col = "",
    show_legend = TRUE
    ){
  pca <- prcomp(t(na.omit(df)))
  pca_df <- pca$x %>%
    as.data.frame() %>%
    rownames_to_column(quantitative_col_name) %>% 
    left_join(batch_info, by = quantitative_col_name)
  var_expl <- pca$sdev^2 / sum(pca$sdev^2)
  names(var_expl) <- paste0("PC", 1:length(var_expl))

  if(shape_col != ""){
    pca_plot <- pca_df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = col_col, shape = shape_col))
  } else {
    pca_plot <- pca_df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = col_col))
  }

  pca_plot <- pca_plot + 
    geom_point(size=3) +
    theme_classic() +
    labs(title = title,
         x = glue::glue("PC1 [{round(var_expl['PC1']*100, 2)}%]"),
         y = glue::glue("PC2 [{round(var_expl['PC2']*100, 2)}%]"))

  if(!show_legend){
    pca_plot <- pca_plot + 
      theme(legend.position = "none")
  }

  if (path == "") {
    return(pca_plot)
  } else {
    ggsave(path, pca_plot, width = 5, height = 5)
    return(pca_plot)
  }
}

# boxplot
boxplot_pg <- function(protein_matrix, metadata_df, quantitativeColumnName, color_col, title, path="") {
  # Reshape data into long format
  long_data <- tidyr::gather(protein_matrix, 
                             key = "file", value = "Intensity")
  merged_data <- merge(long_data, metadata_df, by.x = "file", by.y = quantitativeColumnName)
  
  # Log tranformed scale
  boxplot <- ggplot(merged_data, aes(x = file, y = Intensity, fill = .data[[color_col]])) + 
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    # adjust fonsize for the x-axis
    theme(axis.text.x = element_text(size = 8)) +
    labs(title = title) 

  if(path == "") {
    return(boxplot)
  } else {
      ggsave(path, boxplot)
      return(boxplot)
  }
}

plotIntensityDensityByPool <- function(
    intensities_df, metadata_df, quantitativeColumnName, poolColumnName, title
) {
  # Reshape the intensities_df from wide to long format
  long_intensities <- reshape2::melt(intensities_df, 
    variable.name = "Sample", value.name = "Intensity")
  
  # Adjust the merge function based on your metadata column names
  merged_data <- merge(long_intensities, metadata_df, by.x = "Sample", by.y = quantitativeColumnName)
  
  # Plot the data
  ggplot(merged_data, aes(x = Intensity, color = .data[[poolColumnName]])) +  
    geom_density() +
    theme_minimal() +
    labs(title = paste(title, " by", poolColumnName),
         x = "Intensity",
         y = "Density")
}


heatmap_plot <- function(pg_matrix, batch_info, name, condition="condition", lab="lab"){
    cor_matrix <- cor(na.omit(pg_matrix), use = "pairwise.complete.obs")
    resulting_plot <- ggpubr::as_ggplot(grid::grid.grabExpr(
        pheatmap::pheatmap(cor_matrix, 
                        annotation_col = select(batch_info, c(condition, lab)),
                        treeheight_row = 0, treeheight_col = 0, 
                        main = paste0(name, ' heatmap')
        )
      )
    )
    return(resulting_plot)
}


heatmap_nocor_plot <- function(pg_matrix, batch_info, name, condition="condition", lab="lab", use_breaks=TRUE){
    
    if(use_breaks){
      breaks = seq(-3, 3, length.out = 101)
    } else {
      breaks = NA
    }
    resulting_plot <- pheatmap::pheatmap(pg_matrix, 
                        annotation_col = select(batch_info, c(condition, lab)),
                        treeheight_row = 0, treeheight_col = 0,
                        # do not cluster the rows and cols
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        main = paste0(name, ' heatmap'),
                        # set limits to -3,3
                        breaks = breaks,
        )
    return(resulting_plot)
}
