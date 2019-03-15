## Information content

* calculate the information content for GO terms as -log(p), where p is the probability a gene has an annotation, i.e. the number of genes with the annotation divided by the total number of genes

Requires:
* GO hiearchy (go-basic.obo)
* GO annotations (goa_human.gaf without header)
* GO namespace (default C, one of C, F, P)

1. Run script
```
"$CMSCRIPTS"/information-content/index.js --annotations=goa_human.gaf --obo=go-basic.obo --namespace=C
```

2. Output
* `information-content.txt`
* console: tier bins

### IC bins

* create information content tier bins to split terms into the following four levels
* the total number of genes can be determined by the number of genes with the term `Cellular_component`

| tier | % of genes with term | no. genes | lower IC limit (inclusive) | upper IC limit (exclusive) |
|------|----------------------------|----------------------------|----------------------|-----------|
| 1 | 1% | | | |
| 2 | 2% | | | |
| 3 | 10% | | | |
| 4 | 25% | | | |
| 5 | 100% | | | |
