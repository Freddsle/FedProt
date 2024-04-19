library(ggrepel)
library(ggplot2)

volcano_plot <- function(result, title, show_names=TRUE, show_legend=TRUE) {
  result$P.Value.log <- -log10(result$adj.P.Val)

  volcano_plot <- ggplot(result, aes(x = logFC, y = P.Value.log, color = abs(logFC) > 0.25 & P.Value.log > -log10(0.05))) +
    # add name to the color legeng
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                        name = "DE proteins") +
    # add line for p-value threshold (0.05) and logFC threshold (1)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_point() +
    theme_minimal() +
    xlab("Log Fold Change") +
    ylab("-log10 P-value") +
    ggtitle(paste0("Volcano plot for limma results, ", title))

  # if show lgene is False, remove the legend
  if(show_legend == FALSE){
    volcano_plot <- volcano_plot + theme(legend.position = "none")
  }

  if(show_names == TRUE){
    # add names for PG that passed the thresholds -- using rownames as labels
    volcano_plot <- volcano_plot + 
      geom_text_repel(data = result[result$adj.P.Val < 0.05 & abs(result$logFC) > 0.25,], 
                    aes(label = rownames(result[result$adj.P.Val < 0.05 & abs(result$logFC) > 0.25,])), 
                    size = 3) 
  }
  return(volcano_plot)
}