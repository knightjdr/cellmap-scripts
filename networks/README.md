## Networks

### Summary files for ranks/domains

Generate ranks.json/domains.json file that contains all details about each rank/domain, including its listed name on the cellmap, matching go term etc, all GO terms with p-values and domains with p-values

Requires:
* summary file
  * rank-summary.txt
  * domain-summary.txt
* GO terms for eacn rank with (and p-values for NMF)
 * SAFE: attribute_properties_annotation-highest.txt (need to delete the first 4 rows of this file)
 * NMF: terms_perrank.txt
* Domain, motif and/or diesease file associated with each rank
* NMF will require a file with the order for ranks (rank_order.txt)
* SAFE will require a map of go terms to ids (go-map.txt), produced from GO formatting step.

1. Run script

NMF
```
"$CMSCRIPTS"/network/rank-profile.pl -g terms_perrank.txt -r rank-summary.txt -s rank_order.txt -t n -diseases diseases.txt -domains domains.txt -motifs motifs.txt
```

SAFE
```
"$CMSCRIPTS"/network/rank-profile.pl -g terms_perrank.txt -o go-map.txt -r rank-summary.txt -t  -diseases diseases.txt -domains domains.txt -motifs motifs.txt
```

2. Output
* `ranks.json`: NMF rank details
* `domains.json`: SAFE domain details

### NMF tSNE network

Requires:
* coordinates file from tsne (tsne_nmf.txt)
* file with rank assigned to each gene (gene-localizations.txt)
* ranks.json file
* SAINT file (to get list of baits)
* annotations matrix used for SAFE (make sure this matches GO namespace of network). This file indicates whether a GO term is known or not for each gene (use this as opposed to goa-human.gaf because it will contain all parents as well as the GO terms themselves)

1. Run script
```
"$CMSCRIPTS"/network/network-tsne.pl -a cc_annotations.txt -g gene-localization.txt -n tsne_nmf.txt -r ranks.json -s saint.txt -t n
```

2. Output
* nmf-tsne-network-legend.json
* nmf-tsne-network.json

### Correlation network

Requires:
* network in .cyjs format
* file with rank/domain assigned to each gene
  * NMF: gene-localizations.txt
  * SAFE: node_properites_annotation-highest.txt, need to remove top 4 lines
* ranks.json file
* SAINT file (to get list of baits)
* annotations matrix used for SAFE (make sure this matches GO namespace of network). This file indicates whether a GO term is known or not for each gene (use this as opposed to goa-human.gaf because it will contain all parents as well as the GO terms themselves)

1. Run script
```
"$CMSCRIPTS"/network/network-corr.pl -a cc_annotations.txt -g gene-localization.txt -n tsne_nmf.txt -r ranks.json -s saint.txt -t n
```

2. Output
* corr-network-legend.json
* corr-network.json
