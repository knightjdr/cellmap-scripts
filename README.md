# Cell map analysis scripts

Scripts used for analyzing data at the cell-map.org. In this guide, all scripts are assumed to be in $HOME/cellmap-scripts. Adjust path as needed. All perl scripts can be run with the -h flag for command line options.

## SAINT file processing

Scripts are in `saint-processing`

### Merge SAINT files

Requires:
* SAINT files
* bait list

1. Place all SAINT files (ending in .txt) into folder with no other files.
2. Create a file with a list of baits. Only needs one column but can have more.

| bait  | other |
|-------|-------|
| bait1 | a     |
| bait2 | b     |
| ...   | ...   |

3. Run saint_merge.pl in the directory with the SAINT files and specifiy path to list of baits. In this example, the assumption is the bait list is one directory above the folder with the SAINT files.
```
$HOME/cellmap-scripts/saint-processing/saint-merge.pl -b ../baits.txt
```

4. Outputs:
* merged.txt conating the merged SAINT files
* missing_merged.txt with a list of any baits in the bait list file that were not found
* skipped_merged.txt to report any baits in the SAINT files that are not being included

### Filter SAINT file

Requires:
* SAINT file
* contaminants list
* FDR cutoff, default 0.01
* minimum prey threshold, default 5

Contaminant list is in `data-files/contaminants.txt`. It was generated from:
* ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta 
* http://www.coxdocs.org/doku.php?id=maxquant:common:download_and_installation and grabbing contaminants.fasta
* supplemented with BirA tag, streptavidin, mCherry, GFP, bovine, casein and adenoviral proteins

Script will filter contaminants, decoys and any baits that do not have at least 5 preys passing the cutoff of 0.01 FDR. These thresholds can be adjusted with the -m and -f flags respectively.

1. Run script using contaminants.txt and SAINT file.
```
$HOME/cellmap-scripts/saint-processing/saint-filter.pl -c contaminants.txt -s saint.txt
```

2. Output:
* saint_filtered.txt
* filtered_baits.txt to list any baits removed

Note: Need to manually remove Human prey genes found with mouse baits

## Data set assessment

Scripts are in `assessment`

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

## Interaction assessment

Scripts are in `interaction-assessment`

### Recovered ineractions per FDR cutoff

This script will take all interaction pairs (BioGRID, Intact or merged). For each FDR probability it will calculate the percentage of recovered known interactions.

Requires:
* SAINT file
* interaction file from BioGRID, Intact or merged
* optional: a list of all baits to include in analysis (default will assume all baits)
* bait to gene name map

| bait  | gene  |
|-------|-------|
| bait1 | genea |
| bait2 | geneb |
| ...   | ...   |

1. Run script
```
$HOME/cellmap-scripts/interaction-assessment/complex-validation.pl -b interactions.txt -m bait-gene.txt -s saint.txt -t m
```

2.	Output:
* fraction-recovered.txt

### Bait assessment

Generate a bait summary using an FDR of 0.01 when comparing against BioGrid (or Intact or merged).

Requires
* SAINT file
* interaction file from BioGRID, Intact or merged
* FDR cutoff, default 0.01
* bait to gene name map

| bait  | gene  |
|-------|-------|
| bait1 | genea |
| bait2 | geneb |
| ...   | ...   |

1. Run script
```
$HOME/cellmap-scripts/interaction-assessment/bait-summary.pl -b interactions.txt -m bait-gene.txt -s saint.txt -t m
```

2.	Output:
* bait-level-recovery.txt
* edge-recovered.txt for plotting known bait-prey interaction partners on Cytoscape

### Bait assessment, top 25 preys

* SAINT file
* interaction file from BioGRID, Intact or merged
* FDR cutoff, default 0.01
* adjust spectral counts to length, default is yes
* subtract control average from spectral counts, default is yes
* number of top preys to use, default is 25
* bait to gene name map

| bait  | gene  |
|-------|-------|
| bait1 | genea |
| bait2 | geneb |
| ...   | ...   |

1. Run script
```
$HOME/cellmap-scripts/interaction-assessment/bait-summary-topX.pl -b interactions.txt -m bait-gene.txt -s saint.txt -t m
```

2.	Output:
* bait-level-recovery_top25_lengthAdjusted.txt

## Dataset formatting

Scripts are in `data-sets`

### BioGRID

1. Download BioGRID dataset
* https://thebiogrid.org/download.php)
* go to Current Release
* download file BIOGRID-ORGANISM-X.X.X.tab2.zip, unzip and grab the human file

### IntAct

1. Download IntAct dataset
* http://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip

2. Parse for human genes
```
$HOME/cellmap-scripts/data-sets/intact_parse.pl -i intact.txt
```

3. Output:
* parsed-intact.txt

4. Merge BioGRID and IntAct datasets
```
$HOME/cellmap-scripts/data-sets/merge-biogrid-intact.pl -b BIOGRID-ORGANISM-Homo_sapiens-X.X.X.tab2.txt -i parsed-intact.txt
```

5. Output:
* interactions.txt
