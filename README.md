# Cell map analysis scripts

Scripts used for analyzing data at the cell-map.org. In this guide, all scripts are assumed to be in $HOME/cellmap-scripts. Adjust path as needed. All perl scripts can be run with the -h flag for command line options.

## Script order

For third-party datasets, see `data-sets` README for instructions on getting and/or reformatting files.

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

4. prey-prey
* prey-prey correlation
* generate list of preys in network
* generate Cytoscape files
* prey-prey interaction recovery
* prey-prey complex recovery

5. cytoscape
* format correlation file for Cytoscape
* generate network for SAFE

6. safe
* create annotation matrix
* run SAFE with a collection of networks, annotation matrices and parameters
* assess results
* generate Cytoscape network coloured by SAFE domain

7. nmf
* create matrix for NMF
* run NMF in batch
* summarize NMF ranks by GO terms
* assess NMF ranks
* assess prey moonlighting
* generate heat map for viewing NMF matrices
* generate Cytoscape network based on correlation
* create t-SNE map from basis matrix
