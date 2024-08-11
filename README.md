# FedProt

Privacy-Preserving Federated Multi-center Differential Protein Abundance Analysis tool. 

https://arxiv.org/abs/2407.15220

![logo](https://github.com/Freddsle/FedProt/blob/main/logo2.jpg?raw=true)

It is a federated version of state of the art [DEqMS](https://pubmed.ncbi.nlm.nih.gov/32205417/) workflow. 

Current implementaion is available for DIA-LFQ and DDA-TMT MS data, as well as for any other data type that does not require additional preprocessing. 

Available normalization methods:
- log2 transformation;
- median normalization across all clients (for TMT data);
- IRS normalization inside each client (for TMT data).

## Content:

- [FedProt](#fedprot)
  - [Content:](#content)
- [Config Settings](#config-settings)
  - [Configuration File Description](#configuration-file-description)
- [Running the app](#running-the-app)
  - [Prerequisite](#prerequisite)
  - [Test data](#test-data)
  - [Run](#run)
  - [Output](#output)
- [FedProt states](#fedprot-states)
- [Evaluation](#evaluation)
  - [Run the evaluation:](#run-the-evaluation)
- [Repo structure:](#repo-structure)
- [Citation:](#citation)
  - [Preprint:](#preprint)
    - [BibTeX](#bibtex)

# Config Settings 

Here is an example config file. Optional parameters are marked.  
In data folder you can find example structure of client's data and config files.

```
fedprot:
  counts: protein_counts.tsv       
  design: design.tsv       
  intensities: protein_groups_matrix.tsv 

  sep: '\t' 
  
  use_smpc: true            # or false, control using SMPC
  
  max_na_rate: 0.8
  log_transformed: false

  experiment_type: 'TMT'            # default is "DIA"
  ref_type: 'in_silico_reference'   # optional, only if "TMT" specified
  plex_covariate: true              # optional, only if "TMT" specified
  plex_column: "Pool"               # optional, only if "TMT" specified

  use_median: true                  # default - false
  use_irs: true                     # default - false
  
  use_counts: true
  only_shared_proteins: false
  
  remove_single_pep_protein: true   # default - false  
  remove_single_value_design: true  # default - true, privacy filter
  TEST_MODE: false                  # optional, default - false. Use only for testing - for skip additional privacy protection filters.

  target_classes: [...]             # example - ["heathy", "FSGS"]
  covariates: []

  result_table: "DPE.csv"
```

## Configuration File Description

This configuration file is used to configure a federated proteomics analysis pipeline, particularly for data processing and analysis in mass spectrometry experiments. Below is a detailed explanation of each parameter in the `fedprot` configuration:

**Input Files**: 
- **`counts: protein_counts.tsv`**
  - Optional, use if "use_counts" set to true. Specifies the path to the protein group counts file, tab-separated values (.tsv) format, with two columns: one for protein groups (PG) and one for counts. 
- **`design: design.tsv`**
  - Path to the design file. The first column should contain sample names, and the columns representing target classes should be with boolean values (0 or 1).
- **`intensities: protein_groups_matrix.tsv`**
  - This file contains the intensity values for protein groups across samples, with rows representing protein groups and columns representing samples. Default is in .tsv format.

**Data Formatting**:  
- **`sep: '\t'`**
  - Specifies that all input files are tab-separated.

**Processing Options**:  
- **`use_smpc: true`**
  - Determines whether Secure Multi-Party Computation (SMPC) is used. Set to `true` to enable SMPC, or `false` to disable it.
- **`max_na_rate: 0.8`**
  - Specifies the maximum proportion of missing values (NA) allowed within each target class. A value of 0.8 allows up to 80% missing data per class.
- **`log_transformed: false`**
  - Indicates whether the data has already been log-transformed. Set to `true` if the data is already log-transformed; otherwise, set to `false`.

**Experiment Type Settings**:  
- **`experiment_type: 'TMT'`**
  - Specifies the type of experiment being conducted. The default is "DIA" for DIA and other data types, could be "TMT" (Tandem Mass Tagging) - specifically for TMT data.

- **`ref_type: 'in_silico_reference'`**
  - This option is applicable only if `experiment_type` is set to "TMT". It specifies the type of reference used, with "in_silico_reference" being the only tested option for now.

- **`plex_covariate: true`**
  - Relevant only for "TMT" experiments. Set to `true` if multiple TMT-plexes are present within a single client.

- **`plex_column: "Pool"`**
  - Specifies the name of the design file column containing TMT-plex information. Ensure that TMT-plexes names are not repeated between clients.

**Normalization Options**:  
- **`use_median: true`**
  - Enables median normalization if set to `true`. The default setting is `false`.
- **`use_irs: true`**
  - Enables Internal Reference Scaling (IRS) normalization if set to `true`. The default setting is `false`.

**Protein Filtering**:  
- **`remove_single_pep_protein: true`**
  - If set to `true`, proteins identified by only a single peptide will be removed, which improves data quality. The default setting is `false`.  
- **`remove_single_value_design: true`**
  - A privacy filter that transforms any single non-NA value in a design column subgroup to NA to protect privacy. The default setting is `true`.
- **`TEST_MODE: false`**
  - Optional setting for testing purposes. If set to `true`, additional privacy protection filters are skipped. The default setting is `false`.

**Analysis Options**:  
- **`use_counts: true`**
  - Determines whether protein group counts will be used in the analysis. Set to `true` to include counts.
- **`only_shared_proteins: false`**
  - If set to `true`, only proteins detected in all samples will be included in the analysis. The default setting is `false`.

**Target Classes and Covariates**:  
- **`target_classes: [...]`**
  - This field should contain a list of target classes, such as `["healthy", "FSGS"]`. These represent the groups under study.

- **`covariates: []`**
  - A list of covariates can be included here to adjust for in the analysis.

**Output**:  
- **`result_table: "DPE.csv"`**
  - Specifies the filename for the output results table, which will be saved as a .csv file.


# Running the app

## Prerequisite

To run FedProt app, Docker and FeatureCloud pip package should be installed:

```shell
pip install featurecloud
```

Start controller. 
```shell
# first, create and go the the dir, where test folder will be created
cd path/to/dir/with/test
featurecloud controller start --data-dir=PATH/TO/DATA/data/bacterial_data/balanced
```

Download (or build locally) the app:

```shell
# download
featurecloud app download featurecloud.ai/fedprot

# OR build
featurecloud app build featurecloud.ai/fedprot
```
## Test data

You can find example test data in:
- `data/bacterial_data/balanced` - clients lab_A,lab_B,lab_C,lab_D,lab_E,
- `data/TMT_data/01_smaller_lib_balanced_PG_MajorPG` - Center1,Center2,Center3
- and `data/simulated_data/balanced`, `data/simulated_data/mild_imbalanced`, `data/simulated_data/imbalanced` - clients lab1,lab2,lab3.


## Run

You can run FedProt as a standalone app in the FeatureCloud test-bed [FeatureCloud test-bed](https://featurecloud.ai/development/test), or you can also run the app using CLI:

```shell
featurecloud test start --controller-host=http://localhost:8000 --app-image=fedprot --query-interval=1 --client-dirs=lab_A,lab_B,lab_C,lab_D,lab_E
```

The results could be found in the featurecloud tests folder. 

You can use provided example data or you own data.

## Output

The results file contains logFC, p-values and adj.p-values or count-adjusted p-values nd adj.p-values. The result file is in the same format as DEqMS or limma result tables.

# FedProt states

The FedProt app states:

![states](https://github.com/Freddsle/FedProt/blob/main/states.png?raw=true)

Represent the client-side workflow, red represent the coordinator's workflow, and violet represent transitions that involve both the client and coordinator. 

More about it you can read in the [FedProt paper](https://arxiv.org/abs/2407.15220) methods.

# Evaluation

Required libraries to run FedProt evaluation code:
1. R:
   ```
   DEqMS
   tidyverse        # purrr, dplyr, ggplot2
   ggrepel
   gridExtra
   grid
   patchwork
   data.table
   foreach
   MetaVolcanoR
   RankProd
   invgamma
   RobNorm 
   ```
2. Python:
   ```
   pandas
   numpy
   scipy
   statsmodels
   matplotlib
   yaml
   upsetplot
   ```

## Run the evaluation:

You can more quickly familiarize yourself with how FedProt works by using the  `evaluation_utils/fedprot_prototype/fedprot_script.py` script.
Be aware that this version does not have SMPC and runs locally, only as an introduction and test.

The examples and evaluation is in `evaluation` folder. Evaluation was done using 5 datasets, two real-world: bacterial DIA-LFQ and human plasma DDA-TMT, and 3 simulated.

# Repo structure:

For real datasets - in `evaluation/TMT_data/` and `evaluation/bacterial/` data folders - code to run the analysis (central, FedProt, meta-analyses). The code for evaluation and plot figures based on the results are in `evaluation/aggregated_eval/` folder.

Analysis of 'Handling of batch effects' are in `evaluation/batch_effects_eval/` folder.

Code for the simulated data analysis and evaluation are in `evaluation/simulated/`. Only final aggregated results are present.

# Citation:

## Preprint:
### BibTeX
```
@misc{burankova2024privacypreservingmulticenterdifferentialprotein,
      title={Privacy-Preserving Multi-Center Differential Protein Abundance Analysis with FedProt}, 
      author={Yuliya Burankova and Miriam Abele and Mohammad Bakhtiari and Christine von Törne and Teresa Barth and Lisa Schweizer and Pieter Giesbertz and Johannes R. Schmidt and Stefan Kalkhof and Janina Müller-Deile and Peter A van Veelen and Yassene Mohammed and Elke Hammer and Lis Arend and Klaudia Adamowicz and Tanja Laske and Anne Hartebrodt and Tobias Frisch and Chen Meng and Julian Matschinske and Julian Späth and Richard Röttger and Veit Schwämmle and Stefanie M. Hauck and Stefan Lichtenthaler and Axel Imhof and Matthias Mann and Christina Ludwig and Bernhard Kuster and Jan Baumbach and Olga Zolotareva},
      year={2024},
      eprint={2407.15220},
      archivePrefix={arXiv},
      primaryClass={q-bio.QM},
      url={https://arxiv.org/abs/2407.15220}, 
}
```
