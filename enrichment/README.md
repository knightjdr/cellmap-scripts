## Enrichment

For each NMF rank/SAFE domain, I calculate a fold-enrichment, along with its p-value and then perform the Benjamini-Hochberg procedure to see what passes a 5% FDR.

### Domains

Requires:
* Pfam domains
* a list of genes with associated ranks/domains
 * NMF: gene-localizations.txt
 * SAFE: node_properites_annotation-highest.txt which lists the best domain for each gene. Need to remove top 4 lines from this file after SAFE and only keep the columns 1 + 3.
* rank/domain details (optional: if this file is included, each page on the spreadsheet will show the rank name to make life easier, but if absent the tab number will represent the rank)
  * NMF: rank-summary.txt
  * SAFE: domain-summary.txt respectively
* type of analysis (default n; n for NMF or s for SAFE)

1. Run script
```
"$CMSCRIPTS"/enrichment/domain-enrichment.pl -g gene-localizations.txt -n rank-summary.txt -u pfam-domains.txt -t n
```

2. Ouput
* domain.xls for manual inspection (domains that pass the 1% FDR). Open once with excel and save as .xlsx
* domain.txt for use in generating the network file for the Cellmap

### Things other than domains

Requires:
* File with enrichment terms: a file with two columns listing the gene name and then the term for enrichment (one term per row, so a gene with multiple terms will have multiple rows)
* a list of genes with associated ranks/domains
 * NMF: gene-localizations.txt
 * SAFE: node_properites_annotation-highest.txt which lists the best domain for each gene. Need to remove top 4 lines from this file after SAFE and only keep the columns 1 + 3.
* rank/domain details (optional: if this file is included, each page on the spreadsheet will show the rank name to make life easier, but if absent the tab number will represent the rank)
  * NMF: rank-summary.txt
  * SAFE: safe_summary.txt respectively
* type of analysis (default n; n for NMF or s for SAFE)

1. Run script
```
"$CMSCRIPTS"/enrichment/general-enrichment.pl -g gene-localizations.txt -n rank-summary.txt -u enrichment-terms.txt -t n
```

2. Ouput
* enrichment.xls for manual inspection (terms that pass the 1% FDR). Open once with excel and save as .xlsx
* enrichment.txt for use in generating the network file for the Cellmap
