# Cell map analysis scripts

Scripts used for analyzing data at the cell-map.org. In this guide, all scripts are assumed to be in $HOME/cellmap-scripts. Adjust path as needed. All perl scripts can be run with the -h flag for command line options.

## Script order

For third-party datasets, see `data-sets` README for instructions.

1. saint-processing
* merge files
* filter files to remove contaminants and poorly performing baits

2. assessment
* summary statistics: no. preys, unique preys, etc.
* r-squared between replicates
* bait GO enrichment
* bait localization summary (from known information)
* baits recovered as preys

3. interaction-assessment
* known interaction recovery per FDR
* known and new interactions recovered per bait
