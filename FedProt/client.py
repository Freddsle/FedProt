import pandas as pd
import numpy as np
import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)

EXPERIMENT_TYPE = "TMT"
SAMPLE_TYPE = "sample_type"
TMT_PLEX = "TMT-plex"
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
    ):
        self.experiment_type = experiment_type
        self.tmt_names = None if self.experiment_type == EXPERIMENT_TYPE else None
        self.cohort_name = cohort_name
        self.intensities = None
        self.design = None
        self.prot_names = None
        self.sample_names = None

        self.variables = None

        # df for filters
        self.counts = None
        self.pep_counts = None

        if not self.open_dataset(intensities_file_path, 
                                 count_file_path, 
                                 annotation_file_path, 
                                 log_transformed=log_transformed):
            raise Exception("Failed to open dataset")

        self.XtX = None
        self.XtY = None
        self.SSE = None
        self.cov_coef = None
        self.fitted_logcounts = None
        self.mu = None

    def open_dataset(self, intensities_file_path, count_file_path, design_file_path, count_pep_file_path=None, log_transformed=False):
        """
        For LFQ-data:
        Reads data and design matrices and ensures that sample names are the same.
        Log2(x + 1) transforms intensities.

        For TMT-data:
        Reads data and design matrices and ensures that sample names are the same,
        each TMT-plex has at least one reference sample. Excludes proteins detected in less than 'min_f' plexes.
        Log2(x+1) transforms intensities.
        """
        self.read_files(intensities_file_path, count_file_path, design_file_path, count_pep_file_path)

        if not self.process_files(log_transformed=log_transformed):
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

    def process_files(self, log_transformed=False):
        """Process the loaded data based on experiment type."""
        if self.experiment_type == EXPERIMENT_TYPE:
            if not SAMPLE_TYPE in self.design.columns:
                logging.error(f"Client {self.cohort_name}: Design matrix does not contain '{SAMPLE_TYPE}' column.")
                return False
            if not self.process_tmt_files():
                return False

        # if intensities not log2 transformed
        # TODO: ADD filter ot parameter!
        if not log_transformed:
            self.intensities = np.log2(self.intensities + 1)        
            logging.info(f"Client {self.cohort_name}: Log2(x+1) transformed intensities.")
        else:
            logging.info(f"Client {self.cohort_name}: Intensities are already log2 transformed.")
        
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

    def process_tmt_files(self):
        """Validate and process the TMT files."""
        self.validate_tmt_files()

        for tmt in self.tmt_names:
            if REF_SAMPLE not in self.design.loc[self.design[TMT_PLEX] == tmt, SAMPLE_TYPE].values:
                logging.error(
                    f"Client {self.cohort_name}: No reference sample found in TMT-plex {tmt}. All samples will be excluded."
                )
                self.tmt_names.discard(tmt)

        self.counts = self.counts.loc[:, self.tmt_names]
        self.design = self.design.loc[self.design[TMT_PLEX].isin(self.tmt_names), :]

    def validate_tmt_files(self):
        """
        Validates the TMT files.
        """
        if TMT_PLEX not in self.design.columns:
            logging.error(f"Client {self.cohort_name}: Design matrix does not contain '{TMT_PLEX}' column.")
            return

        self.design[TMT_PLEX] = self.design[TMT_PLEX].apply(lambda x: str(x) + "_" + self.cohort_name)
        self.counts.rename(lambda x: str(x) + "_" + self.cohort_name, axis="columns", inplace=True)
        self.tmt_names = set(self.design[TMT_PLEX].values)

        if not set(self.counts.columns.values) == self.tmt_names:
            shared_tmt_plexes = set(self.counts.columns.values).intersection(self.tmt_names)
            logging.error(
                f"Client {self.cohort_name}: Only {len(shared_tmt_plexes)} TMT-plexes are shared between design matrix and count table."
            )
            self.tmt_names = shared_tmt_plexes
        self.tmt_names = sorted(list(self.tmt_names))

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
            f"Client {self.cohort_name}: Validated {len(self.sample_names)} samples and {len(self.prot_names)} proteins."
        )

    def validate_protein_names(self, stored_features):
        """
        Ensure that gene names are the same and are in the same order.
        """
        global_prots = set(stored_features)
        self_prots = set(self.prot_names).intersection(global_prots)
        #self_prots = set(self.prot_names)
        
        if len(self_prots) != len(set(self_prots)):
            logging.info("Client %s:\tDuplicate protein names found." % self.cohort_name)

        if self_prots != global_prots:
            extra_prots = self_prots.difference(global_prots)
            if len(extra_prots) > 0:
                logging.info(
                    "Client %s:\t%s protein groups absent in other datasets are dropped:"
                    % (self.cohort_name, len(extra_prots))
                )
            else:
                extra_prots = global_prots.difference(self_prots)
                logging.info("Client %s:\t%s protein groups not found." % (self.cohort_name, len(extra_prots)))

        if not self_prots.issubset(global_prots):
            extra_prots = self_prots.difference(global_prots)
            logging.info("Client %s:\t%s protein groups not found in Global proteins. Data load failed." % (self.cohort_name, len(extra_prots)))
            raise ValueError(
                f"Client {self.cohort_name}: Some proteins are missing in the global protein list: {extra_prots}"
            )    
        # reorder genes
        self.prot_names = list(self_prots)
        self.intensities = self.intensities.loc[self.prot_names, :]
        if self.use_counts:
            self.counts = self.counts.loc[self.prot_names, :]

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
                    "Client %s:\t%s columns are excluded from the design matrix:" % (self.cohort_name, len(extra_variables))
                )

            # keep only necessary columns in the design matrix
            if self.experiment_type == EXPERIMENT_TYPE:
                self.design = self.design.loc[:, variables + [TMT_PLEX, SAMPLE_TYPE]]
            else:
                self.design = self.design.loc[:, variables]

        # find how many TMT detect each protein group
        if self.experiment_type == EXPERIMENT_TYPE:
            self.n_tmt_per_prot = self.n_tmt - self.counts.isna().sum(axis=1)
        
    
    def add_cohort_effects_to_design(self, cohorts):
        """add covariates to model cohort effects."""
        for cohort in cohorts:
            if self.cohort_name == cohort:
                self.design[cohort] = 1
            else:
                self.design[cohort] = 0

    ######## Median Centering ###########
    def compute_medians(self):
        # computes and stores sample medians and returns their average
        self.sample_medians = self.intensities.median()
        return np.mean(self.sample_medians)

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
            d = self.design.loc[self.design[TMT_PLEX] == tmt_plex, :]
            tmt_samples = self.intensities.loc[:, d.index.values]
            refs = d.loc[d[SAMPLE_TYPE] == REF_SAMPLE, :].index.values
            avg_ref = self.intensities.loc[:, refs].mean(axis=1)
            norm_intensities.append(tmt_samples.subtract(avg_ref - avg_ref_global, axis=0))
        norm_intensities = pd.concat(norm_intensities, axis=1)
        self.intensities = norm_intensities.loc[:, self.design.index.values]

    ######### Filtering ##################
    def remove_single_pept_prots(self):
        """
        Replaces with NA all intensities for proteins supported by a single peptide
        """
        # TODO: need to be tested
        n_pepts = self.counts.shape[0]

        if self.experiment_type == EXPERIMENT_TYPE:
            for tmt in self.counts.columns.values:
                samples = self.design.loc[self.design[TMT_PLEX] == tmt, :].index.values
                failed_prots = self.counts[tmt][self.counts[tmt] <= 1].index.values
                # replace internsities with NA
                self.intensities.loc[failed_prots, samples] = np.NaN
                self.counts.loc[self.counts[tmt] <= 1] = np.NaN
        else:
            failed_prots = self.counts[self.counts <= 1].index.values
            self.intensities.loc[failed_prots, :] = np.NaN
            self.counts.loc[self.counts <= 1] = np.NaN

        # drop lines with all NA and update the number of TMT plexes per protein
        self.counts.dropna(axis=0, how="all", inplace=True)
        self.intensities = self.intensities.loc[self.counts.index, :]

        self.update_prot_names()

        logging.info(
            f"Client {self.cohort_name}:\tProteins supported by just a single peptide will be excluded: {n_pepts - self.counts.shape[0]}"
        )

    def is_all_na(self, sample_type="sample"):
        """
        If TMT data: For each protein group in each TMT returns 1 if all samples of sample_type are Na.
        """
        intensities = self.intensities
        design = self.design

        if self.experiment_type == EXPERIMENT_TYPE:
            samples = design.loc[design[SAMPLE_TYPE] == sample_type, :].index.values
            intensities = intensities.loc[:, samples]
            design = design.loc[samples, :]
            not_detected_in_tmt = {}

            for tmt in set(design[TMT_PLEX].values):
                tmt_samples = design.loc[design[TMT_PLEX] == tmt, :].index.values
                na_fraction = intensities.loc[:, tmt_samples].isna().sum(axis=1) * 1.0 / len(tmt_samples)
                # mark 1 peptides not detected in TMT
                not_found = na_fraction
                not_found[not_found < 1] = 0
                not_detected_in_tmt[tmt] = not_found

            not_detected_in_tmt = pd.DataFrame.from_dict(not_detected_in_tmt)
            return not_detected_in_tmt

        else:
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

    def apply_filters(self, min_f=0.5, remove_single_peptide_prots=False, sample_type="sample"):
        if remove_single_peptide_prots:
            self.remove_single_pept_prots()

        if self.experiment_type == EXPERIMENT_TYPE:
            # marks proteins not detected in refs
            # not_detected_in_refs = self.is_all_na(sample_type="ref")

            # marks proteins not detected in samples
            not_detected_in_samples = self.is_all_na(sample_type)
            not_detected = not_detected_in_samples  # + not_detected_in_refs
            not_detected[not_detected > 0] = 1
            not_detected = not_detected.sum(axis=1) / not_detected.shape[1]
            passed_prots = not_detected[not_detected <= min_f].index.values

            logging.info(
                f"Client {self.cohort_name}:\tProtein groups detected in less than {min_f} of all TMT-plexes will be excluded: {self.counts.shape[0] - len(passed_prots)}"
            )
            
            return self.prot_names
        
        else:           
            na_count_in_variable, samples_per_class = self.is_all_na(sample_type)
            logging.info(
                f"Client {self.cohort_name}:\tProtein groups detected in less than {min_f} of each target class will be excluded:"
            )    
            return na_count_in_variable, samples_per_class

    def update_prot_names(self, passed_prots=None):
        if passed_prots is None:
            logging.error(f"Client {self.cohort_name}: No proteins passed the filters.")
            raise ValueError(f"Client {self.cohort_name}: No proteins passed the filters.")
        else:
            self.prot_names = passed_prots
            if self.use_counts:
                self.counts = self.counts.loc[self.prot_names, :]
            self.intensities = self.intensities.loc[self.prot_names, self.design.index]

        if self.experiment_type == EXPERIMENT_TYPE:
            self.n_tmt_per_prot = self.counts.notna().sum(axis=1)

        logging.info(f'Samples in {self.cohort_name} data: {len(self.sample_names)}, protein groups: {len(self.prot_names)}')


    def prepare_for_limma(self, stored_features):
        # remove unnecessary columns from the design matrix
        if self.experiment_type == EXPERIMENT_TYPE:
            self.design = self.design.loc[self.design[SAMPLE_TYPE] == SAMPLE_TYPE, :]
            self.design = self.design.drop(columns=[TMT_PLEX, SAMPLE_TYPE])
        else:
            self.prot_names = stored_features
    
        self.sample_names = self.design.index.values
        self.intensities = self.intensities.loc[self.prot_names, self.sample_names]
        self.n_samples = len(self.sample_names)

    ####### limma: linear regression #########
    def compute_XtX_XtY(self):
        X = self.design.values
        Y = self.intensities.values  # Y - intensities (proteins x samples)
        n = Y.shape[0]  # genes
        k = self.design.shape[1]  # variables
        self.XtX = np.zeros((n, k, k))
        self.XtY = np.zeros((n, k))

        # linear models for each row
        for i in range(0, n):  #
            y = Y[i, :]
            # check NA in Y
            ndxs = np.argwhere(np.isfinite(y)).reshape(-1)
            if len(ndxs) != len(y):
                # remove NA from Y and from x
                x = X[ndxs, :]
                y = y[ndxs]
            else:
                x = X
            self.XtX[i, :, :] = x.T @ x
            self.XtY[i, :] = x.T @ y
        return self.XtX, self.XtY

    def compute_SSE_and_cov_coef(self, beta):
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
            x = X[ndxs, :]
            y = y[ndxs]
            self.mu[i, ndxs] = x @ beta[i, :]  # fitted Y
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

