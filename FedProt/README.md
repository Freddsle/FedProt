# FedProt

Privacy-Preserving Federated Differential Protein Expression Analysis tool.

It is a federated version of state of the art [DEqMS](https://pubmed.ncbi.nlm.nih.gov/32205417/) workflow. 

Current implementaion is available for DIA MS data.

## Config 

```
fedprot:
  counts: protein_counts.tsv
  design: design.tsv
  intensities: protein_groups_matrix.tsv
  sep: '\t'
  
  use_smpc: true
  
  max_na_rate: 0.8
  log_transformed: false

    # label: "Condition"
  target_classes: [...] 
  covariates: []

  result_table: "DPE.csv"
```

## Running app
You can run FedProt as a standalone app in the FeatureCloud test-bed. You can also run the app using CLI:

```
featurecloud test start --app-image featurecloud.ai/fedprot --client-dirs './c1,./c2,./c3' --generic-dir './generic'
```

# Publication
...