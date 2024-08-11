# chmod +x central_run.R

# Description: This script runs DE analysis on the central data for the bacterial, human serum and simulated datasets.
# Load necessary libraries
source("../../evaluation_utils/evaluation/DE_analysis.R")
source("../../evaluation_utils/plots/DE_plots.R")
source("../../evaluation_utils/filtering/filtering_normalization.R")


##############################################################################################################3
### Bacterial dataset ### 
##############################################################################################################3

args <- commandArgs(trailingOnly = TRUE)
# Assume args are ordered as dataset1_cohorts, dataset2_cohorts, dataset3_cohorts
cohorts_bacterial <- as.vector(strsplit(args[1], ",")[[1]])
cohorts_TMT <- as.vector(strsplit(args[2], ",")[[1]])
cohorts_simulated <- as.vector(strsplit(args[3], ",")[[1]])


path_to_reports <- "../../data/bacterial_data/balanced/"
  # cohorts <- c("lab_A", "lab_B", "lab_C", "lab_D", "lab_E")
cohorts <- cohorts_bacterial
output_file <- "../../evaluation/fc_results/DEP_bacterial_CENTRAL.tsv"


central_intensities <- NULL
central_counts <- NULL
central_batch_info <- NULL

for (dataset in cohorts) {
  batch_info <- read_tsv(paste0(path_to_reports, dataset, "/metadata_short.tsv"), show_col_types = FALSE)
  intensities <- read_tsv(paste0(path_to_reports, dataset, "/protein_groups_matrix.tsv"), show_col_types = FALSE)
  counts <- read_tsv(paste0(path_to_reports, dataset, "/protein_counts.tsv"), show_col_types = FALSE)

  if (is.null(central_intensities)) {
    central_intensities <- intensities
    central_counts <- counts
    central_batch_info <- batch_info
  } else {
    central_intensities <- full_join(central_intensities, intensities, by = 'rowname')
    central_counts <- full_join(central_counts, counts, by = 'rowname')
    central_batch_info <- rbind(central_batch_info, batch_info)
  }
  cat('Dataset: ', dataset, "\n")
  cat('\tNumber of proteins: ', nrow(central_intensities), '\n')
  cat('\tNumber of samples: ', ncol(central_intensities)-1, '\n\n')
}

central_batch_info <- central_batch_info %>%
  mutate(lab = factor(lab, levels = cohorts),
  condition = as.factor(condition))
cat("Levels of condition: ", levels(central_batch_info$condition), "\n")
cat("Levels of lab: ", levels(central_batch_info$lab), "\n")

central_intensities <- central_intensities %>% column_to_rownames('rowname')
central_counts <- central_counts %>% column_to_rownames('rowname')
central_intensities <- central_intensities[, central_batch_info$file]

central_intensities <- filter_by_condition(central_intensities, central_batch_info, 
    'file', c('Glu', 'Pyr'), 'condition')
central_intensities <- filter_na_proteins(central_intensities, central_batch_info, "file")

central_counts$count <- apply(central_counts, 1, min, na.rm = TRUE)
central_counts <- central_counts %>% select(count) %>% as.data.frame()

cat("Rows after all filters:", nrow(central_intensities), "\n")
central_intensities <- log2(central_intensities + 1)

design <- make_design(central_batch_info, 'condition', 'lab')
contrasts <- makeContrasts(Glu-Pyr, levels = colnames(design))
de_results <- run_DE(central_intensities, central_counts, design, contrasts)
de_results <- de_results %>% rownames_to_column('Protein')
write.table(de_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)



##############################################################################################################3
### Human serum dataset ### 
##############################################################################################################3

path_to_reports <- "../../data/TMT_data/01_smaller_lib_balanced_PG_MajorPG/"
# cohorts <- c("Center1", "Center2", "Center3")
cohorts <- cohorts_TMT
output_file <- "../../evaluation/fc_results/DEP_TMT_CENTRAL.tsv"

central_batch_info = read_tsv(paste0(path_to_reports, 'metadata.tsv'), show_col_types = FALSE)
central_batch_info <- central_batch_info[central_batch_info$Group != "Common Reference" & central_batch_info$Group != "in_silico", ]
central_batch_info <- central_batch_info  %>%
  filter(Group != "Common Reference") %>%
  mutate(
      file = Quantitative.column.name,
      lab = as.factor(Center),
      condition = factor(Group, levels = c('heathy', 'FSGS')),
      Pool = as.factor(Pool)
  )
cat("Number of samples: ", nrow(central_batch_info), "\n")

central_intensities <- NULL
central_counts <- NULL

for (dataset in cohorts) {
  intensities <- read_tsv(paste0(path_to_reports, dataset, "/pg_intensities.tsv"), show_col_types = FALSE)
  counts <- read_tsv(paste0(path_to_reports, dataset, "/pg_counts.tsv"), show_col_types = FALSE)

  if (is.null(central_intensities)) {
    central_intensities <- intensities
    central_counts <- counts
  } else {
    central_intensities <- full_join(central_intensities, intensities, by = 'Majority.protein.IDs')
    central_counts <- full_join(central_counts, counts, by = 'Majority.protein.IDs')
  }

}
cat('Dataset: ', dataset, "\n")
cat('\tNumber of proteins: ', nrow(central_intensities), '\n')
cat('\tNumber of samples: ', ncol(central_intensities)-1, '\n\n')

# transform and order
central_intensities <- central_intensities %>% column_to_rownames('Majority.protein.IDs')
central_counts <- central_counts %>% column_to_rownames('Majority.protein.IDs')
central_intensities <- central_intensities[, central_batch_info$file]

# select minimal count across column for each protein (with na.rm = TRUE)
central_counts$count <- apply(central_counts, 1, min, na.rm = TRUE)
central_counts <- central_counts %>% select(count)

# Filter - remove proteins with 1 count 
central_intensities <- central_intensities[rownames(central_intensities) %in% rownames(filter(central_counts, count > 1)), ]
cat("Rows after filtering 1 count:", nrow(central_intensities), "\n")
# filter by condition
central_intensities <- filter_na_proteins(central_intensities, central_batch_info, "file")
central_intensities <- filter_by_condition(central_intensities, central_batch_info, 
    'file', c('heathy', 'FSGS'), 'condition')

central_intensities <- central_intensities[, central_batch_info$file]
cat("Rows after all filters:", nrow(central_intensities), "\n")

# median normalization
central_intensities <- medianNorm(central_intensities) %>% as.data.frame()

# IRS normalization
corrected_intensities <- NULL

for(center in cohorts){
    cat("IRS normalization for center: ", center, "\n")
    center_intensities <- central_intensities[, central_batch_info$lab == center]
    center_batch_info <- central_batch_info[central_batch_info$lab == center, ]
    cat("Shape of the data: ", dim(center_intensities), "\n")

    irs_res <- irsNorm_in_silico_single_center(
        center_intensities, center_batch_info,
        pool_col = "Pool",
        column_name = "file",
        center = name,
        aggregation_method = "average",
        add_refs = FALSE
    )
    # add to the central data
    if(is.null(corrected_intensities)){
        corrected_intensities = irs_res$corrected_data
    } else {
        corrected_intensities = cbind(corrected_intensities, irs_res$corrected_data)
    }
}

central_intensities <- corrected_intensities[, central_batch_info$file]
central_intensities <- log2(central_intensities + 1)
cat("Number of samples for analysis: ", ncol(central_intensities), "\n")
cat("Number of proteins for analysis: ", nrow(central_intensities), "\n")

central_batch_info <- central_batch_info %>% mutate(
    condition = factor(condition, levels = c('heathy', 'FSGS')),
    Pool = as.factor(Pool)
)

# run DE analysis
design <- make_design(central_batch_info, 'condition', 'Pool')
contrasts <- makeContrasts(heathy-FSGS, levels = colnames(design))
de_results <- run_DE(central_intensities, central_counts, design, contrasts)
de_results <- de_results %>% rownames_to_column('Protein')
write.table(de_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)


##############################################################################################################3
### Simulated balanced, mild_imbalanced and imbalanced ### 
##############################################################################################################3

path_to_reports <- "../../data/simulated_data/"
# cohorts <- c("lab1", "lab2", "lab3")
cohorts <- cohorts_simulated
output_file <- "../../evaluation/fc_results/"

for (mode in c("balanced", "mild_imbalanced", "imbalanced")){
  generated_data_directory <- paste0(path_to_reports, mode)

  batch_info <- read.table(paste0(generated_data_directory, "/data/batch_info.tsv"), sep = "\t", header = TRUE)
  batch_info <- batch_info %>%
      mutate(lab = as.factor(batch), condition = as.factor(condition))
  # read data
  lab1_data <- read.table(paste0(generated_data_directory, "/lab1/intensities.tsv"), sep = "\t", header = TRUE, row.names = 1)
  lab2_data <- read.table(paste0(generated_data_directory, "/lab2/intensities.tsv"), sep = "\t", header = TRUE, row.names = 1)
  lab3_data <- read.table(paste0(generated_data_directory, "/lab3/intensities.tsv"), sep = "\t", header = TRUE, row.names = 1)

  # run DEqMS for joined data    
  design <- model.matrix(~0 + batch_info$condition + batch_info$lab)
  colnames(design) <- c(levels(batch_info$condition), levels(batch_info$lab)[-1])
  rownames(design) <- batch_info$file

# full join to keep all proteins
  joined_data <- merge(
          merge(lab1_data, lab2_data, by = "row.names", all = TRUE) %>% column_to_rownames("Row.names"),
      lab3_data,, by = "row.names", all = TRUE) %>% column_to_rownames("Row.names")
  joined_data <- joined_data[,rownames(design)]

  joined_data <- filter_na_proteins(joined_data, batch_info, 'file')
  joined_data <- filter_by_condition(joined_data, batch_info, 
      'file', c('A', 'B'), 'condition')

  contrasts <- makeContrasts(A - B, levels = colnames(design))
  de_results <- run_DE(joined_data, NULL, design, contrasts)
  de_results <- de_results %>% rownames_to_column('Protein')
  if (mode == "balanced"){
    output_file_name = paste0(output_file, "/DEP_simulBal_CENTRAL.tsv")
    write.table(de_results, file = output_file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (mode == "mild_imbalanced") {
    output_file_name = paste0(output_file, "/DEP_simulMildImbal_CENTRAL.tsv")
    write.table(de_results, file = output_file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (mode == "imbalanced") {
    output_file_name = paste0(output_file, "/DEP_simuImb_CENTRAL.tsv")
    write.table(de_results, file = output_file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
