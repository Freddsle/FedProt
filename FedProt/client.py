import pandas as pd
import numpy as np
import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

EXPERIMENT_TYPE = "TMT"
SAMPLE_TYPE = "sample_type"
REF_SAMPLE = "ref"


class Client:
    def __init__(
        self,
        cohort_name,
        intensities_file_path,
        count_file_path,    
        annotation_file_path,
        experiment_type=EXPERIMENT_TYPE,
        log_transformed=False,
        ref_type=None,
        plex_column=None,
        target_classes=None,
    ):
        self.experiment_type = experiment_type
        self.cohort_name = cohort_name
        self.intensities = None
        self.design = None
        self.prot_names = None
        self.sample_names = None

        self.target_classes = target_classes
        
        self.variables = None

        if self.experiment_type == EXPERIMENT_TYPE:
            self.ref_type = ref_type
            self.tmt_names = None
            if self.ref_type != "in_silico_reference":
                self.plex_column = "TMT-plex"
            else:
                self.plex_column = plex_column       
                
        self.log_transformed = log_transformed     

        # df for filters
        self.counts = None
        self.pep_counts = None

        if not self.open_dataset(intensities_file_path, 
                                 count_file_path, 
                                 annotation_file_path):
            raise Exception("Failed to open dataset")
        
        self.check_collinearity = False
        self.XtX = None
        self.XtY = None
        self.SSE = None
        self.cov_coef = None
        self.fitted_logcounts = None
        self.mu = None

    def open_dataset(self, intensities_file_path, count_file_path, design_file_path, count_pep_file_path=None, ):
        """
        For LFQ-data:
        Reads data and design matrices and ensures that sample names are the same.
        Log2(x + 1) transforms intensities.

        For TMT-data:
        Reads data and design matrices and ensures that sample names are the same,
        Log2(x+1) transforms intensities.
        """
        self.read_files(intensities_file_path, count_file_path, design_file_path, count_pep_file_path)

        if not self.process_files():
            logging.error(f"Client {self.cohort_name}: Failed to process files.")
            return False
        
        self.check_and_reorder_samples()

        # remove rows with less than 2 non-NA values
        self.intensities = self.intensities.dropna(axis=0, thresh=2)
        self.prot_names = list(self.intensities.index.values)

        # Check if the data has only one sample or one protein
        # TODO: add parameter to ignore it
        if len(self.sample_names) < 2:
            logging.error(f"Client {self.cohort_name}: Data has only one sample.")
            return False
        if len(self.prot_names) < 2:
            logging.error(f"Client {self.cohort_name}: Data has only one protein.")
            return False

        logging.info(
            f"Client {self.cohort_name}: Loaded {len(self.sample_names)} samples and {len(self.prot_names)} proteins."
        )
        return True

    def read_files(self, intensities_file_path, count_file_path, design_file_path, count_pep_file_path=None):
        """Read files using pandas and store the information in the class attributes."""
        self.intensities = pd.read_csv(intensities_file_path, sep="\t", index_col=0)
        self.sample_names = list(self.intensities.columns.values)

        if not count_file_path:
            logging.info(f"Client {self.cohort_name}: Not using counts.")
            self.use_counts = False
        else:
            logging.info(f"Client {self.cohort_name}: Using counts.")
            self.use_counts = True
            self.counts = pd.read_csv(count_file_path, sep="\t", index_col=0)
            if count_pep_file_path:
                self.pep_counts = pd.read_csv(count_pep_file_path, sep="\t")

        self.design = pd.read_csv(design_file_path, sep="\t", index_col=0)

    def process_files(self):
        """Process the loaded data based on experiment type."""
        if self.experiment_type == EXPERIMENT_TYPE:
            if not SAMPLE_TYPE in self.design.columns and self.ref_type != "in_silico_reference":
                logging.error(f"Client {self.cohort_name}: Design matrix does not contain '{SAMPLE_TYPE}' column.")
                return False
            if not self.process_tmt_files():
                return False
        
        # If any row contains only one non-NA value in a row, replace it with NA
        logging.info(f"Client {self.cohort_name}: Replacing rows with only one non-NA value with NA.")
        logging.info(f"Client {self.cohort_name}: Rows will be affected : {self.intensities.apply(lambda x: x.count() == 1, axis=1).sum()}")
        self.intensities = self.intensities.apply(lambda x: x if x.count() > 1 else np.nan, axis=1)

        # if intensities not log2 transformed
        if not self.log_transformed and self.experiment_type != EXPERIMENT_TYPE:
            self.intensities = np.log2(self.intensities + 1)
            self.log_transformed = True
            logging.info(f"Client {self.cohort_name}: Log2(x+1) transformed intensities.")
        elif self.experiment_type == EXPERIMENT_TYPE and not self.log_transformed:
            logging.info(f"Client {self.cohort_name}: TMT data will be log2 transformed after normalization.") 
        elif self.log_transformed:
            logging.info(f"Client {self.cohort_name}: Intensities are already log2 transformed.")
        else:
            logging.error(f"Client {self.cohort_name}: Failed to transform intensities.")
        
        return True
    
    def check_and_reorder_samples(self):
        """Check and reorder samples based on design and intensities."""
        design_rows = set(self.design.index.values)
        exprs_cols = set(self.intensities.columns.values)
        presented_samples = list(design_rows.intersection(exprs_cols))
        
        # log message for samples that are not in both design and intensities
        if len(presented_samples) < len(self.design.index.values):
            logging.warning(f"Client {self.cohort_name}: Design matrix contains samples that are not in intensities.")
        if len(presented_samples) < len(self.intensities.columns.values):
            logging.warning(f"Client {self.cohort_name}: Intensities contain samples that are not in design matrix.")

        # keep only sample names that are in both design and intensities
        self.sample_names = [name for name in self.design.index.values if name in presented_samples]

        self.design = self.design.loc[self.sample_names, :]
        self.intensities = self.intensities.loc[:, self.sample_names]
        self.n_samples = len(self.sample_names)

        # remove samples that are not belong to any of target classes
        if self.target_classes:
            self.design = self.design.loc[self.design[self.target_classes].sum(axis=1) > 0, :]
            self.sample_names = list(self.design.index.values)
            self.intensities = self.intensities.loc[:, self.sample_names]
            self.n_samples = len(self.sample_names)
            logging.info(f"Client {self.cohort_name}: Samples are filtered based on target classes.")
            logging.info(f"Client {self.cohort_name}: {self.n_samples} samples are kept.")

    def process_tmt_files(self):
        """Validate and process the TMT files."""
        self.validate_tmt_files()

        if self.ref_type != "in_silico_reference":
            for tmt in self.tmt_names:
                if REF_SAMPLE not in self.design.loc[self.design[self.plex_column] == tmt, SAMPLE_TYPE].values:
                    logging.error(
                        f"Client {self.cohort_name}: No reference sample found in TMT-plex {tmt}. All samples will be excluded."
                    )
                    self.tmt_names.discard(tmt)

            self.counts = self.counts.loc[:, self.tmt_names]
            self.design = self.design.loc[self.design[self.plex_column].isin(self.tmt_names), :]
        
        logging.info(f"Client {self.cohort_name}: TMT data loaded successfully.")
        return True

    def validate_tmt_files(self):
        """
        Validates the TMT files.
        """
        if self.plex_column not in self.design.columns:
            logging.error(f"Client {self.cohort_name}: Design matrix does not contain '{self.plex_column}' column.")
            return
        logging.info(f"Client {self.cohort_name}: Using TMT-plex column '{self.plex_column}'.")

        if self.ref_type != "in_silico_reference":
            self.design[self.plex_column] = self.design[self.plex_column].apply(lambda x: str(x) + "_" + self.cohort_name)
            self.counts.rename(lambda x: str(x) + "_" + self.cohort_name, axis="columns", inplace=True)
            self.tmt_names = set(self.design[self.plex_column].values)
            if not set(self.counts.columns.values) == self.tmt_names:
                shared_tmt_plexes = set(self.counts.columns.values).intersection(self.tmt_names)
                logging.error(
                    f"Client {self.cohort_name}: Only {len(shared_tmt_plexes)} TMT-plexes are shared between design matrix and count table."
                )
                self.tmt_names = shared_tmt_plexes
        else:
            self.tmt_names = set(self.design[self.plex_column].values)

            # if any row has only one non-NA value inside TMT-plex, replace the value with NA
            # for plexes
            for tmt in self.tmt_names:
                tmt_samples = self.design.loc[self.design[self.plex_column] == tmt, :].index.values
                logging.info(f"Client {self.cohort_name}: Checking TMT-plex {tmt} for single non-NA values.")
                affected_rows = self.intensities.loc[:, tmt_samples].apply(lambda x: x.count() == 1, axis=1).sum()
                if affected_rows > 0:
                    logging.info(f"Client {self.cohort_name}: Rows will be affected: {self.intensities.loc[:, tmt_samples].apply(lambda x: x.count() == 1, axis=1).sum()}")
                    self.intensities.loc[:, tmt_samples] = self.intensities.loc[:, tmt_samples].apply(
                        lambda x: x if x.count() > 1 else np.nan, axis=1
                    )

        self.tmt_names = sorted(list(self.tmt_names))
        logging.info(f"Client {self.cohort_name}: Found {len(self.tmt_names)} TMT-plexes. Plexes: {self.tmt_names}")
        
        return True

    def validate_inputs(self, stored_features, variables):
        """
        Checks if protein_names match global protein group names.
        Important to ensure that client's protein names are in the global protein names. But it's not important that all global protein names are in the client's protein names.
        """
        # store variables for filtering
        self.variables = variables
        self.validate_protein_names(stored_features)
        self.validate_variables(variables)
        logging.info(
            f"Client {self.cohort_name}:\tValidated {len(self.sample_names)} samples and {len(self.prot_names)} proteins."
        )

    def validate_protein_names(self, stored_features):
        """
        Ensure that gene names are the same and are in the same order.
        """
        global_prots = set(stored_features)
        self_prots = set(self.prot_names)

        if self_prots > global_prots:
            logging.error(
                f"Client {self.cohort_name}:\tSome protein groups are not in the global list: {len(self_prots - global_prots)}"
            )
            raise ValueError(f"Client {self.cohort_name}:\tSome protein groups are not in the global list.")
        
        elif self_prots < global_prots:
            logging.info(
                f"Client {self.cohort_name}:\tSome protein groups are not in the client list: {len(global_prots - self_prots)} and will be added."
            )

        self_prots = global_prots

        # reorder genes
        self.prot_names = list(self_prots)
        self.intensities = self.intensities.reindex(self.prot_names)
        if self.use_counts:
            self.counts = self.counts.loc[self.prot_names]
        logging.info(f"Client {self.cohort_name}:\tProtein groups are validated.")

    def validate_variables(self, variables):
        """ensure that design matrix contains all variables"""

        self_variables = set(self.design.columns)

        if self_variables != set(variables):
            missing_variables = set(variables).difference(self_variables)
            if len(missing_variables) > 0:
                raise ValueError(
                    f"Client {self.cohort_name}: Some variables are missing in the design matrix: {missing_variables}"
                )

            extra_variables = self_variables.difference(set(variables))
            if len(extra_variables) > 0:
                logging.info(
                    f"Client {self.cohort_name}:\t{len(extra_variables)} columns are excluded from the design matrix: {extra_variables}"
                )

            # keep only necessary columns in the design matrix
            if self.experiment_type == EXPERIMENT_TYPE:
                if self.ref_type != "in_silico_reference":
                    self.design = self.design.loc[:, variables + [self.plex_column, SAMPLE_TYPE]]
                else:
                    self.design = self.design.loc[:, variables + [self.plex_column]]
            else:
                self.design = self.design.loc[:, variables]

        # find how many TMT detect each protein group
        if self.experiment_type == EXPERIMENT_TYPE and self.ref_type != "in_silico_reference":
            self.n_tmt_per_prot = self.n_tmt - self.counts.isna().sum(axis=1)

    def add_cohort_effects_to_design(self, cohorts, plex_covariate=False):
        """add covariates to model cohort effects."""
        logging.info(f"Client {self.cohort_name}:\tAdding cohort effects to the design matrix.")
        reference_col = cohorts[0]
        cohorts = cohorts[1:]

        if self.experiment_type == EXPERIMENT_TYPE and plex_covariate:
            logging.info(f"Client {self.cohort_name}:\tAdding TMT-plex as a covariate.")
            # use cohort as plex covariate list, check using self.plex_column
            for cohort in cohorts:
                # for each row in design check value of plex column,
                # if it is equal to cohort, set 1, otherwise 0
                self.design[cohort] = self.design[self.plex_column].apply(lambda x: 1 if x == cohort else 0)
    
        else:
            logging.info(f"Client {self.cohort_name}:\tAdding cohort as a covariate.")
            for cohort in cohorts:
                if self.cohort_name == cohort:
                    self.design[cohort] = 1
                else:
                    self.design[cohort] = 0

        # if reference column belongs to the client, set self.check_collinearity = True
        if reference_col == self.cohort_name:
            self.check_collinearity = True
            self.coll_samples = self.design.loc[self.design[reference_col] == 1, :].index.values
            logging.info(f"Client {self.cohort_name}:\tCollinearity will be checked.")
        elif self.experiment_type == EXPERIMENT_TYPE and plex_covariate:
            if reference_col in self.tmt_names:
                self.check_collinearity = True
                self.coll_samples = self.design.loc[self.design[self.plex_column] == reference_col, :].index.values
                logging.info(f"Client {self.cohort_name}:\tCollinearity will be checked.")
        logging.info(f"Client {self.cohort_name}:\tCohort effects are added to the design matrix.")

    ######### Filtering ##################
    def remove_single_pept_prots(self):
        """
        Remove from intensities and counts proteins supported by just a single peptide.
        """
        logging.info(f"Client {self.cohort_name}:\tProtein groups: {len(self.intensities.index)}")
        logging.info(f"Client {self.cohort_name}:\tProtein groups supported by a single peptide will be excluded.")
        proteins_passing_filter = self.counts[self.counts > 1].index
        self.intensities = self.intensities.loc[proteins_passing_filter]
        self.counts = self.counts[proteins_passing_filter]
        logging.info(f"Client {self.cohort_name}:\tProtein groups after filter: {len(self.intensities.index)}")

    def check_not_na(self, sample_type="sample"):
        """
        If TMT data: For each protein group in each TMT returns 1 if all samples of sample_type are Na.
        """
        intensities = self.intensities
        design = self.design

        if self.experiment_type == EXPERIMENT_TYPE and self.ref_type != "in_silico_reference":
            samples = design.loc[design[SAMPLE_TYPE] == sample_type, :].index.values
            intensities = intensities.loc[:, samples]
            design = design.loc[samples, :]
            not_detected_in_tmt = {}

            for tmt in set(design[self.plex_column].values):
                tmt_samples = design.loc[design[self.plex_column] == tmt, :].index.values
                na_fraction = intensities.loc[:, tmt_samples].isna().sum(axis=1) * 1.0 / len(tmt_samples)
                # mark 1 peptides not detected in TMT
                not_found = na_fraction
                not_found[not_found < 1] = 0
                not_detected_in_tmt[tmt] = not_found

            not_detected_in_tmt = pd.DataFrame.from_dict(not_detected_in_tmt)
            return not_detected_in_tmt

        na_count_in_variable = pd.DataFrame(index=self.intensities.index.values, columns=self.variables)
        samples_per_class = {}

        for variable in self.variables:
            samples = design.loc[design[variable] == 1, :].index.values
            tmp_intensities = intensities.loc[:, samples]
            # na_fraction = tmp_intensities.isna().sum(axis=1) * 1.0 / len(samples)
            na_count = tmp_intensities.isna().sum(axis=1)
            na_count_in_variable[variable] = na_count
            # add the number of samples per class
            samples_per_class[variable] = len(samples)

        return na_count_in_variable, samples_per_class

    def apply_filters(self, min_f=0.5, remove_single_peptide_prots=False):
        sample_type = SAMPLE_TYPE

        if remove_single_peptide_prots:
            self.remove_single_pept_prots()

        if self.experiment_type == EXPERIMENT_TYPE and self.ref_type != "in_silico_reference":
            # marks proteins not detected in refs
            # not_detected_in_refs = self.is_all_na(sample_type="ref")

            # marks proteins not detected in samples
            not_detected_in_samples = self.check_not_na(sample_type)
            not_detected = not_detected_in_samples  # + not_detected_in_refs
            not_detected[not_detected > 0] = 1
            not_detected = not_detected.sum(axis=1) / not_detected.shape[1]
            passed_prots = not_detected[not_detected <= min_f].index.values
            logging.info(
                f"Client {self.cohort_name}:\tProtein groups detected in less than {min_f} of all TMT-plexes will be excluded: {self.counts.shape[0] - len(passed_prots)}"
            )            
            return self.prot_names            
        
        else:           
            na_count_in_variable, samples_per_class = self.check_not_na(sample_type)
            logging.info(
                f"Client {self.cohort_name}:\tProtein groups detected in less than {min_f} of each target class will be excluded."
            )    
            return na_count_in_variable, samples_per_class

    def update_prot_names(self, passed_prots=None):
        if passed_prots is None:
            logging.error(f"Client {self.cohort_name}: No proteins passed the filters.")
            raise ValueError(f"Client {self.cohort_name}: No proteins passed the filters.")
        else:
            self.prot_names = passed_prots
            if self.use_counts:
                self.counts = self.counts[self.prot_names]
            self.intensities = self.intensities.loc[self.prot_names, self.design.index]

        if self.experiment_type == EXPERIMENT_TYPE and self.ref_type != "in_silico_reference":
            self.n_tmt_per_prot = self.counts.notna().sum(axis=1)

        logging.info(f'Samples in {self.cohort_name} data: {len(self.sample_names)}, protein groups: {len(self.prot_names)}')


    ######## Median Centering ###########
    def compute_medians(self):
        # computes and stores sample medians and returns their average
        self.sample_medians = self.intensities.median(skipna=True)
        return np.mean(self.sample_medians)

    def mean_median_centering(self, global_mean_median):
        # Based on weighted mean of samples medians
        self.intensities = self.intensities / self.sample_medians
        self.intensities = self.intensities * global_mean_median
        logging.info(f"Client {self.cohort_name}:\tMedian centering is applied.")

    # different version
    def median_centering(self, global_median):
        # R: sweep(log2(df.intensities+1), 2, sample_medians-global_median)
        self.intensities = self.intensities - self.sample_medians + global_median

    def count_leq_gt_samples(self, x):
        # leq - less or eqal x ;
        leq = self.intensities[self.intensities <= x].count().sum()
        # gt - grater than x
        gt = self.intensities[self.intensities > x].count().sum()
        return leq, gt

    ######### Row shift ##################
    def get_sum_refs(self):
        refs = list(self.design.loc[self.design[SAMPLE_TYPE] == REF_SAMPLE, :].index.values)
        return self.intensities.loc[:, refs].sum(axis=1)

    def apply_rs(self, avg_ref_global):
        norm_intensities = []
        for tmt_plex in self.tmt_names:
            d = self.design.loc[self.design[self.plex_column] == tmt_plex, :]
            tmt_samples = self.intensities.loc[:, d.index.values]
            refs = d.loc[d[SAMPLE_TYPE] == REF_SAMPLE, :].index.values
            avg_ref = self.intensities.loc[:, refs].mean(axis=1)
            norm_intensities.append(tmt_samples.subtract(avg_ref - avg_ref_global, axis=0))
        norm_intensities = pd.concat(norm_intensities, axis=1)
        self.intensities = norm_intensities.loc[:, self.design.index.values]

    ######### IRS with in silico reference #########

    def get_samples(self, pool):
        # get samples for a pool        
        filtered_samples = self.design[self.design[self.plex_column] == pool].index.values.astype(str).tolist()
        return filtered_samples

    def aggregate_data(self, data, method="average"):
        # Aggregate data based on the method
        if method == "average":
            return data.mean(axis=1, skipna=True)
        else:
            return data.sum(axis=1, skipna=True)

    def calculate_geometric_mean(self, irs):
        def geometric_mean(x):
            non_zero_values = x[x > 0]  # Exclude zeros from the computation
            if len(non_zero_values) == 0:
                return np.nan
            return np.exp(np.mean(np.log(non_zero_values)))

        geometric_means = irs.apply(geometric_mean, axis=1)
        return geometric_means

    def irsNorm_in_silico_single_center(self, aggregation_method="average"):
        # Initialize IRS table
        logging.info(f"Client {self.cohort_name}:\tIRS normalization with in silico reference.")
        irs = pd.DataFrame()

        # Process each pool to create IRS references
        pools = self.design[self.plex_column].unique()
        if len(pools) < 2:
            logging.info(f"Client {self.cohort_name}:\tOnly one pool found. IRS normalization is not applied.")
            return
        
        logging.info(f"Client {self.cohort_name}:\tFound {len(pools)} pools.")
        for pool in pools:
            samples = self.get_samples(pool)

            if not all(sample in self.intensities.columns for sample in samples):
                logging.error(f"Client {self.cohort_name}:\tSome samples for pool {pool} are missing in the data columns.")
                raise ValueError(f"Some samples for pool {pool} are missing in the data columns.")
            
            data = self.intensities[samples].copy()
            processed_data = self.aggregate_data(data, aggregation_method)
            irs[pool] = processed_data

        # Compute geometric mean across all IRS references
        irs_average = self.calculate_geometric_mean(irs)

        # Calculate scaling factors and normalize the self.intensities
        corrected_data = self.intensities.copy()
        
        for pool in pools:
            scaling_factor = irs_average / irs[pool]
            scaling_factor.replace([np.inf, -np.inf], np.nan, inplace=True)
            scaling_factor.fillna(1, inplace=True)
            pool_samples = self.get_samples(pool)

            if not all(sample in corrected_data.columns for sample in pool_samples):
                raise ValueError(f"Some samples for pool {pool} are missing in the corrected data columns.")

            corrected_data[pool_samples] = corrected_data[pool_samples].multiply(scaling_factor, axis=0)

        self.intensities = corrected_data
        logging.info(f"Client {self.cohort_name}:\tIRS normalization is applied.")

        return 

    ######### limma: preparation step #########
    def prepare_for_limma(self, stored_features):
        # remove unnecessary columns from the design matrix
        if self.experiment_type == EXPERIMENT_TYPE:
            if self.ref_type != "in_silico_reference":
                self.design = self.design.loc[self.design[SAMPLE_TYPE] == SAMPLE_TYPE, :]
                self.design = self.design.drop(columns=[self.plex_column, SAMPLE_TYPE])
            
            else:
                self.design = self.design.drop(columns=[self.plex_column])
                logging.info(f"Client {self.cohort_name}:\tPlex column is removed from the design matrix.")

            if not self.log_transformed:
                self.intensities = np.log2(self.intensities + 1)
                self.log_transformed = True
                logging.info(f"Client {self.cohort_name}:\tLog2(x+1) transformed intensities.")

        self.prot_names = stored_features
        self.sample_names = self.design.index.values
        self.intensities = self.intensities.loc[self.prot_names, self.sample_names]
        self.n_samples = len(self.sample_names)
        logging.info(f"Client {self.cohort_name}:\tPrepared for limma. Samples: {self.n_samples}, Proteins: {len(self.prot_names)}.")

    def get_mask(self):
        X = self.design.values
        Y = self.intensities.values
        n = Y.shape[0]  # genes
        k = self.design.shape[1]  # variables

        mask_X = np.zeros((n, k))

        for i in range(0, n):
            y = Y[i, :]
            ndxs = np.argwhere(np.isfinite(y)).reshape(-1)
            if len(ndxs) > 0:
                x = X[ndxs, :]
                column_variances = np.var(x, axis=0)
                zero_var_mask = column_variances == 0
                if np.any(zero_var_mask):
                    mask_X[i, zero_var_mask] = 1
            else:
                mask_X[i, :] = 1

        return mask_X

    def updated_mask(self, mask_glob):
        if self.check_collinearity:
            for i, protein in enumerate(self.prot_names):
                # Check if Y values all are NA for the given row
                y_selected = self.intensities.loc[protein, self.coll_samples]
                if np.all(np.isnan(y_selected)):
                    # Find indices of False values in mask_glob[i, len(self.target_classes):]
                    false_indices = np.where(mask_glob[i, len(self.target_classes):] == False)[0]
                    if false_indices.size > 0:
                        # Adjust index to consider the offset from len(self.target_classes)
                        last_false_index = false_indices[-1] + len(self.target_classes)
                        # Set the last False value to True
                        mask_glob[i, last_false_index] = True

            logging.info(f"Client {self.cohort_name}:\tCollinearity check and update completed.")

        mask_glob = mask_glob.astype(int)
        logging.info(f"Client {self.cohort_name}:\tCollinearity check completed.")
        return mask_glob

    ####### limma: linear regression #########
    def compute_XtX_XtY(self):
        X = self.design.values
        Y = self.intensities.values  # Y - intensities (proteins x samples)
        n = Y.shape[0]  # genes
        k = self.design.shape[1]  # variables

        self.XtX = np.zeros((n, k, k))
        self.XtY = np.zeros((n, k))

        # linear models for each row
        for i in range(0, n): 
            y = Y[i, :]
            # check NA in Y
            ndxs = np.argwhere(np.isfinite(y)).reshape(-1)
            if len(ndxs) > 0: 
                x = X[ndxs, :]
                y = y[ndxs]
                self.XtX[i, :, :] = x.T @ x
                self.XtY[i, :] = x.T @ y
                
        return self.XtX, self.XtY

    def compute_SSE_and_cov_coef(self, beta, mask_glob):
        X = self.design.values
        Y = self.intensities.values
        n = Y.shape[0]
        self.SSE = np.zeros(n)
        self.mu = np.empty(Y.shape)

        for i in range(0, n):
            # check NA in Y
            y = Y[i, :]
            # remove NA from Y and from x
            ndxs = np.argwhere(np.isfinite(y)).reshape(-1)
            x = X[ndxs, ]
            x = x[:, ~mask_glob[i]]
            y = y[ndxs]
            self.mu[i, ndxs] = x @ beta[i, ~mask_glob[i]]  # fitted Y
            self.SSE[i] = np.sum((y - self.mu[i, ndxs]) ** 2)  # local SSE
        
        self.cov_coef = X.T @ X
        
        return self.SSE, self.cov_coef

    def get_not_na(self):
        # number of non-NA samples per protein
        n_na = self.intensities.isna().astype(int).sum(axis=1)
        return self.n_samples - n_na

    def sum_intensities(self):
        return self.intensities.sum(axis=1)
    
    # #### DEqMS #####
    def get_min_count(self):
        return self.counts.min(axis=1)

