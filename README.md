# Cell map analysis scripts

> Scripts are currently being migrated to GO and can be found
> at this repo: https://github.com/knightjdr/cmgo

Scripts used for analyzing data at the cell-map.org. All perl scripts can be run with the -h flag for command line options.

## Script order

For third-party datasets, see `data-sets` README for instructions on getting and/or reformatting files.

Setup up an ENV variable `CMSCRIPTS` pointing to this repo. Run any script as `"$CMSCRIPTS"/folder/script.pl`

1. saint-processing
* merge files
* filter files to remove contaminants and poorly performing baits

2. assessment
* summary statistics: no. preys, unique preys, etc.
* r-squared between replicates
* bait GO enrichment
* prey overlap between baits using the Jaccard index
* bait localization summary (from known information)
* baits recovered as preys
* overlap with BioPlex

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
* assess prey localizaiton at varying NMF thresholds
* assess prey moonlighting
* report top three localizations (ranks) per prey
* generate heat map for viewing NMF matrices
* generate Cytoscape network based on correlation
* create t-SNE map from basis matrix

8. localization assessments
* assess NMF, SAFE and published dataset localizations against GO
* compare NMF and SAFE for overlap
* compare two datasets
* calculate the recovery of known terms per bait and avgspec
* assess primary and secondary NMF localizations

9. enrichments
* enriched domains for NMF or SAFE compartments
* enriched motifs for NMF or SAFE compartments

10. organelle assessments
* assess bait overlap for ER
* assess bait overlap for Mito
* assess bait overlap for ER-mito

11. networks (for website)
* create NMF tSNE network
* create correlation networks for NMF or SAFE

12. hierarchy
* build a localization hierarchy for retrieving parent and child terms

13. Upload related
* summary of prey values for evaluating localizations