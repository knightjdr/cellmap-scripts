## Data set assessment

### SAINT file summary

Calculate summary statistics for SAINT file (number of interactions, etc.)

Requires:
* SAINT file
* Bait-Gene map
* BFDR cutoff, default 0.01

The bait gene map should be formatted like so:

| bait  | gene  |
|-------|-------|
| bait1 | genea |
| bait2 | geneb |
| ...   | ...   |

1. Run script
```
$HOME/cellmap-scripts/assessment/saint-summary.pl -m bait-gene.txt -s saint.txt
```

### Reproducibility

Assess reproducibility of replicates (calculate r-squared for replicates using interactors that have passed the FDR cutoff). Assumes two replicates.

Requires:
* SAINT file
* FDR, default 0.01

1. Run script
```
$HOME/cellmap-scripts/assessment/repCheck.R saint.txt 0.01
```

2. Output:
* regression.txt (table with r-squared and adjusted r-squared)
* rsquare_plots.pdf (replicate 1 v replicate 2 plots for all baits)


2. Output:
* summary.txt

### GO term recovery

Takes the significant preys for each bait and gets the enriched terms. Enrichment for each bait will be relative to other preys in the dataset.

Requires:
* SAINT file
* FDR, default 0.01
* Number of top preys to use, default all

1. Run script. If the third argument is a number > 0, only that number of preys for each bait will be used for the enrichment.
```
$HOME/cellmap-scripts/assessment/bait-enrichment.R 0.01
```

2. Output:
* bait-enrichment.xlsx or bait-enrichment_topX.xlsx will have five tabs, one for all GO terms, and separate tabs for BP, CC, MF and Corum