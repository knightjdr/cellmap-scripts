## NMF

### Generate matrix from SAINT file

This script will grab preys that pass a specified FDR for at least one bait and get values across all baits. It will subtract control values and normalize prey counts across baits to 1.

Requires:
* SAINT file
* FDR, default 0.01

1. Run script
```
saintFormat.R saint.txt 0.01
```

2. Output
* sc_matrix.csv

### NMF

Requires:
* Python3
* matrix from previous step, sc-matrix.csv
* rank to start testing from
* total number of ranks to test

1. Run script with caffeinate
```
caffeinate node $HOME/cellmap-scripts/nmf/nmf.js sc-matrix.csv
```

2. Output to `results` folder with folder for each rank containing:
* basis.csv (prey information)
* scores.csv (bait information)

Notes:
* l1_ratio parameter – 1 is better for sparse data
* The init parameter is “nndsvd” which is better for sparsity
* See http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html for options
* can test for sparsity using matrix and the Matlab script `is-sparse.M`

### Processing NMF output for enriched GO terms

This has now been set up in batch to automatically follow in the script nmf.js, so this doesn’t not need to be done separately.

This will use the input prey list as background.

Requires:
* basis and scores matrices from NMF
* GO category to use, default GO:CC (one of GO:CC, GO:BP or GO:MF)

1. Run script
```
HOME/cellmap-scripts/nmf/rankProfile.R basis.csv scores.csv GO:CC
```

2. Output
* summary.xlsx (for inspecting rank information)
* rank_order.txt  (clustering order for ranks, probably won’t use)
* gene-localizations.txt (for each gene, lists its rank)
* terms_perrank.txt (list of GO terms and pvalues for each rank);
