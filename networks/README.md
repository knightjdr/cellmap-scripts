## Networks

### Summary files for ranks/domains

Generate ranks.json/domains.json file that contains all details about each rank/domain, including its listed name on the cellmap, matching go term etc, all GO terms with p-values and domains with p-values

Requires:
* summary file (either rank-summary.txt or domain-summary.txt)
* GO terms for eacn rank with (and p-values for NMF)
 * SAFE: attribute_properties_annotation-highest.txt (need to delete the first 4 rows of this file)
 * NMF: terms_perrank.txt
* Domain, motif and/or diesease file associated with each rank
* NMF will require a file with the order for ranks (rank_order.txt)
* SAFE will require a map of go terms to ids (go-map.txt), produced from GO formatting step.

1. Run script

NMF
```
"$CMSCRIPTS"/network/rank_profile.pl -g terms_perrank.txt -r rank-summary.txt -s rank_order.txt -t n -diseases diseases.txt -domains domains.txt -motifs motifs.txt
```

SAFE
```
"$CMSCRIPTS"/network/rank_profile.pl -g terms_perrank.txt -o go-map.txt -r rank-summary.txt -t  -diseases diseases.txt -domains domains.txt -motifs motifs.txt
```

2. Output
* `ranks.json`: NMF rank details
* `domains.json`: SAFE domain details
