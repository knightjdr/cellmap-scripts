## Localization assessment

### Recovery of known terms

* assess the number of preys assigned to previously known terms
* can be done for NMF, SAFE, the HPA subcellular_localization.csv or any published dataset
* other published datasets must have three columns `gene|localization|GO term`
* all localization terms should be mapped to the closest GO ID
* when a gene is assigned multiple localizations, it quality tier is set using the least informative localization based on the information content

Requries:
* a list of child terms for each parent GO term (`go-children_cc.txt`)
* a file of known associated GO terms for each gene
  * `goa_human.gaf` (no header lines)
  * remove header
* information content for each GO term
* for NMF/SAFE need GO ID(s) corresponding to each rank/domain
  * `rank-summary.txt` for NMF
  * `domain-summary.txt` for SAFE
* localization for each gene
 * NMF: `gene-localizations.txt`
 * SAFE: `node_properties_annotation-highest.txt` (remove header)
* type of input file (n = NMF, s = SAFE, h = HPA, o = other)
* prefix for output file

#### NMF

 1. Run script
 ```
 "$CMSCRIPTS"/localization-assessment/localization-assessment.pl -h go-children.txt -i information-content.txt -l goa_human.gaf -g gene-localizations.txt -r rank-details.txt -t n -o 'nmf-'
 ```

 2.	Output:
 * `nmf-localization-assessment.txt`
 * `node_attributes_cytoscape.txt`: for each gene prints its domain/rank if known (only for SAFE/NMF)

 #### SAFE

 1. Run script
 ```
 "$CMSCRIPTS"/localization-assessment/localization-assessment.pl -h go-children.txt -i information-content.txt -l goa_human.gaf -g node_properties_annotation-highest.txt -r domain-details.txt -t s -o 'safe-'
 ```

 2.	Output:
 * `safe-localization-assessment.txt`
 * `node_attributes_cytoscape.txt`: for each gene prints its domain/rank if known (only for SAFE/NMF)

 #### HPA

 1. Run script
 ```
 "$CMSCRIPTS"/localization-assessment/localization-assessment.pl -h go-children.txt -i information-content.txt -l goa_human.gaf -g subcellular_location.tsv -t h -o 'hpa-'
 ```

 2.	Output:
 * `hpa-localization-assessment.txt`

 #### Other

 1. Run script
 ```
 "$CMSCRIPTS"/localization-assessment/localization-assessment.pl -h go-children.txt -i information-content.txt -l goa_human.gaf -g localizations.txt -t o -o 'other-'
 ```

 2.	Output:
 * `other-localization-assessment.txt`

 ### NMF v SAFE

 * Compute overlap between NMF and SAFE

 Requires:
* NMF rank details (`rank-summary.txt`)
*	SAFE domain details (`domain-summary.txt`)
* file with GO child terms (`go-children_cc.txt`)
* assigned localization per gene:
  * NMF: `gene-localizations.txt`
  * SAFE: `node_properties_annotation-highest.txt` with four header lines removed

1. Run script
```
"$CMSCRIPTS"/localization-assessment/nmf-v-safe.pl -g go-children.txt -n rank-summary.txt -na gene-localizations.txt -s domain-summary.txt -sa node_properties_annotation-highest.txt
```

2. Output:
* `nmf-v-safe.txt`: for each gene it lists the NMF and SAFE compartment and whether they overlap
* outputs to STDOUT whether NMF term matches SAFE (or is a child of), SAFE term matches NMF or either

### Evidence

* assess NMF and SAFE localizations based on number of baits each prey was found with and the average spectral counts it was found with those baits

Requires:
* file with GO child terms (`go-children_cc.txt`)
* known GO annotations (`goa_human.gaf` without header lines)
* GO ID(s) corresponding to each rank/domain
  * NMF: rank-summary.txt
  * SAFE: romain-summary.txt
* assigned localization per gene:
  * NMF: `gene-localizations.txt`
  * SAFE: `node_properties_annotation-highest.txt` with four header lines removed
* SAINT file
* file type (n: NMF or s: SAFE)

1. Run script
```
"$CMSCRIPTS"/localization-assessment/assessment-metrics.pl -h go-children.txt -l goa_human.gaf -r rank-summary.txt -s saint.txt -g gene-localizations.txt -t n
```

2. Output:
* `metrics.txt`
  * first group of rows: for bait numbers 1-15+, list the fraction of preys with known localizations
  * second group of rows: for avgspec values between bins, list the fraction of preys with known localizations. Note (this bins will be listed like 0-4, 5-9, but the upper bound is actually < 5 and < 10 respectively.
  * third group of rows: combines the above two values.

### NMF primary and secondary localizations

* will output the fraction of genes whose primary or secondary localization is previously known
* will also output the number of genes with no known localization apart from:
  * cell
  * cell part
  * cellular_component
  * intracellular
  * intracellular organelle
  * intracellular part
  * organelle
  * organelle part

Requires:
* basis.csv matrix
* file with GO child terms (`go-children_cc.txt`)
* known GO annotations (`goa_human.gaf` without header lines)
* NMF rank information (`rank-summary.txt`)
* threshold for secondary localization (default 0.75, i.e. the secondary localization must be within 0.75 percent of the primary)

1. Run script
```
"$CMSCRIPTS"/localization-assessment/nmf-assessment.pl -b basis.csv -h go-children.txt -l goa_human.gaf -r rank-summary.txt
```

2. Output STDOUT:
* number of genes with no known localization
* number of genes with a known primary localization
* number of genes with a known secondary localization
* number of genes with either a known primary or secondary localization