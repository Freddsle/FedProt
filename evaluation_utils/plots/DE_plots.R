library(ggrepel)
library(ggplot2)


volcano_plot <- function(
    result, title, 
    pval_column = 'sca.adj.pval', 
    pval_threshold = 0.05, logfc_threshold = 1, 
    show_names = TRUE, show_legend = TRUE,
    save_report = FALSE, report_path = NULL
  ){
  result$P.Value.log <- -log10(result[[pval_column]])

  volcano_plot <- ggplot(result, aes(
        x = logFC, y = P.Value.log, 
        color = abs(logFC) > logfc_threshold & P.Value.log > -log10(pval_threshold))
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "DE proteins") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed") +
    geom_point(size = 1, alpha = 0.6) +
    theme_minimal() +
    xlab("Log2FC") +
    ylab(paste0("-log10 ", pval_column)) +
    ggtitle(paste0("Volcano plot, ", title))

  if(show_legend == FALSE) {
    volcano_plot <- volcano_plot + theme(legend.position = "none")
  }

  if(show_names == TRUE) {
    volcano_plot <- volcano_plot + 
      geom_text_repel(data = result[result[[pval_column]] < pval_threshold & abs(result$logFC) > logfc_threshold,], 
                      aes(label = rownames(result[result[[pval_column]] < pval_threshold & abs(result$logFC) > logfc_threshold,])), 
                      size = 3)
  }

  if(save_report == TRUE) {
    ggsave(report_path, volcano_plot)
  } else {
    return(volcano_plot)
  }
}
