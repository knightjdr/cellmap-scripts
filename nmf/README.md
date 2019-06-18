## NMF

### Generate matrix from SAINT file

This script will grab preys that pass a specified FDR for at least one bait and get values across all baits. It will subtract control values and normalize prey counts across baits to 1.

Requires:
* SAINT file
* FDR, default 0.01

1. Run script
```
"$CMSCRIPTS"/nmf/saintFormat.R saint.txt 0.01
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
caffeinate node "$CMSCRIPTS"/nmf/nmf.js sc-matrix.csv
```

2. Output to `Results` folder with folder for each rank containing:
* `basis.csv`: prey information
* `scores.csv`: bait information

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
"$CMSCRIPTS"/nmf/rankProfile.R basis.csv scores.csv GO:CC
```

2. Output
* `summary.xlsx`: for inspecting rank information
* `rank_order.txt`: clustering order for ranks, probably won’t use
* `gene-localizations.txt`: for each gene, lists its rank
* `terms_perrank.txt`: list of GO terms and pvalues for each rank
* `top-rank.txt`: max and median value in each rank

### Move files for assessment

Grab all attribute (terms_perrank.txt) and node files (gene-localizations.txt), rename them to match their parent folders and move them into a folder called “nmf_assessment” and a subfolder called “attributes” or “nodes” depending on their type.

1. Move into folder where NMF was run
```
cd "$CMSCRIPTS"/nmf/analysis
```

2. Run script
```
"$CMSCRIPTS"/nmf/nmf-move-files.pl
```

### Assess NMF terms

Take all terms within each rank (including children), and then for each prey gene it will see if the genes are being assigned to a rank with a previously known term.

Requires
* GO annotations, (goa_human.gaf without header)
* List of children for each GO term
* Map of GO IDs to terms
* GO namespace, default C, one of C, F or P

1. Move into assessment folder
```
cd "$CMSCRIPTS"/nmf/analysis/nmf_assessment
```

2. Run script
```
"$CMSCRIPTS"/nmf/nmf-assessor.pl -g go-children.txt -l goa_human_nohead.gaf -m go-map.txt -n C
```

3. Output
* console

#### Assign representative term for ranks

Use the summary.xlsx from the NMF step to determine the representative term to use for each NMF rank. See file in scripts folder called sample-nmf-terms.txt on how to record this information into a file.

1. Select a representative term. If multiple terms, separate them with comma.

2. Also need to specify the “displayname” to use for the term on the cell map network, as well as the GO identifier, synonyms and the IC content value. Synonyms can be found in go-basic.obo file. 

3. Name the file rank-summary.txt.

### Recovery of known localizations

Calculate the recovery of known localizations for each NMF score. This script can be run a second time with an explicit threshold to generate a file with a boolean indicating if a prey above (inclusive) the threshold or not.

Requires:
* GO annotations, (goa_human.gaf without header)
* List of children for each GO term
* Assigned NMF rank for each prey (gene-localizations.txt)
* Rank details (rank-summary.txt)
* GO namespace, default cc, one of bp, cc or mf
* Threshold below which preys are considered low confidence (default 0)

1. Run script
```
"$CMSCRIPTS"/nmf/nmf-recovery.pl -c go-children.txt -g gene-localizations.txt -l goa_human.gaf -n cc -r rank-summary.txt -t 0.15
```

2. Output
* `nmf-recovery.txt`:
  * for each nmf bin lists the total number of preys in that bin, how many have known localizations and the fraction
  * there are three bin assessments: increments of 0.01, 0.05 and 0.1
  * each bin is inclusive for the lower threshold and exclusive of the upper
  * anything with an NMF value over 1 is recorded as 1
* `prey-confidence.txt`:
  * gene-localizations.txt file
  * column for whether localization is known
  * column for whether score is confident or not

### Robustness box plot

After running cmgo nmf-robustness module, use the file `robustness.txt` with this script to generate boxplots

### Assess prey moonlighting

For each prey calculate whether it has a secondary localization within X% of its primary. X is set from 0-1 in 0.01 increments.

Requires:
* basis.csv

1. Run script
```
"$CMSCRIPTS"/nmf/nmf-moonlighting.pl -b basis.csv
```

2. Output
* moonlighting.txt

### List top three categories for each prey

For each prey output primary, secondary and tertiary rank, along with the score and the difference between the primary/secondary and primary/tertiary.

Requires:
* basis.csv matrix
* NMF summary file (rank-summary.txt)

1. Run script
```
"$CMSCRIPTS"/nmf/top3-ranks.pl -b basis.csv -r rank-summary.txt
```

2. Output
* `top3-ranks.txt`

### Subset prey matrix based on moonlighting

Given two groups of NMF ranks, subset the basis matrix based on preys that are moonlighting in both groups.

Requires:
* basis.csv
* list of ranks in first group
* list of ranks in second group
* moonlighting threshold, default 30(%)

1. Run script
```
"$CMSCRIPTS"/nmf/nmf-moonlighting-subset.pl -b basis.csv -poola 3,15 -poolb 13 -t 30
```

2. Output
* basis-moonlighting.csv

3. Run clustering script on this matrix
```
"$CMSCRIPTS"/nmf/nmf-moonlighting-clustering.R basis-moonlighting.csv
```

4. Output
* clustered_moonlighting_basis.csv, reordered basis matrix

5. Convert CSV to ProHits-viz file
```
HOME/cellmap-scripts/nmf/matrix-to-heatmap.pl -m clustered_moonlighting_basis.csv
```

6. Output
* matrix.tsv, can load in interactive viewer at ProHits-viz

### Generate heat map and Cytoscape network based on correlation

Cluster basis and scores matrices for viewing at Prohits-viz

Requires:
* basis.csv
* scores.csv
* correlation cutoff, default 0.7

1. Run clustering script
```
"$CMSCRIPTS"/nmf/nmf-clustering.R basis.csv scores.csv 0.7
```

2. Output
* clustered_basis.csv
* clustered_scores.csv
* cytoscape_nmf_correlation.tsv for importing to Cytoscape

3. Run script to generate files for prohits-viz
```
"$CMSCRIPTS"/nmf/nmf-heatmap.pl -b clustered_basis.csv -s clustered_scores.csv -r rank-summary.txt
```

4. Output
* pv_basis.tsv
* pv_scores.tsv
* rank-color.tsv for colouring Cytoscape nodes via color attribute (import as table, then for node color, select column and use passthrough). A transparency value will also be included based on the NMF value (NMF = 1, no transparency, NMF = 0, full transparency, with these values scaled from 0-255 for cytoscape).
* cytoscape_nmf.tsv. This can be used to build a network in Cytoscape. If two preys have similar NMF values (default is within 5%) in the same rank, create an edge between them.

### Create t-SNE map

Generate a coordinate file for Cytoscape from the NMF data using t-SNE. Either the basis matrix can be used directly (nmf-wrapper.m) or the Euclidean distance matrix can be calculated from the basis matrix and that used for t-SNE (nmf-wrapper_D.m).

1. Open Matlab and add t-SNE script folder to path: addpath(genpath('"$CMSCRIPTS"/nmf/tsne/'))

2. Open script `nmf/tsne-wrapper.m`

3. Specify path to nmf basis.csv file in script

4. Run script

5. Output
* tsne_nmf.txt. This is a file with X and Y coordinates for each gene.

### Create Cytoscape file from t-SNE map

Requires:
* tsne_nmf.txt, t-SNE map
* rank-color.tsv, colours to use for ranks
* rank-summary.txt, summary file manually created from NMF results

1. Run script
```
"$CMSCRIPTS"/nmf/tsne-to-cytoscape.pl -c rank-color.tsv -d rank-summary.txt -t tsne_nmf.txt
```

2. Output
* tsne_nmf_details.txt. Will have additional columns where x and y are multiplied by 30, 50 and 100. This allows me to pick a good scale for cytoscape.

3. Import `tsne_nmf_details` into cytoscape as network (will only import the first column as nodes)

4. Then import this same file as table

5. On style tab, click properties dropdown and make sure X and Y are visible. Set X, Y and color to respective columns with pass-through mapping

6. Save file as .cys and .cyjs (.cyjs will be used for u+q network on cell-map.org)

### Configure Cytoscape file for u+q on website

Takes .cyjs file from previous step and removes unneeded data attributes.

Requires:
* .cyjs from previous step

1. Run script
```
"$CMSCRIPTS"/nmf/trim-cyjs.pl -c network.cyjs
```

2. Output
* trimmed.cyjs
