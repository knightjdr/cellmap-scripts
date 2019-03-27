## Hierarchy

Annotate an organelle hierarchy. The `hierarchy.json` file contains a hierarchy of parent and child terms. Terms were selected to coincide with GO whenever possible, the main exceptions being the large-scale groupings, such as “membrane-bound”, etc. For each term there is the complete GO name specified by “name” and possible a shorter display name. The script will add the GO ID and two Booleans indicating if the term was seen by NMF or SAFE.

Requries:
* GO term-id map (go-map.txt)
* hierarchy in json format
* NMF details (rank-details.txt)
* SAFE details (domain-details.txt)

1. Run script
```
"$CMSCRIPTS"/hierarchy/hierarchy.pl -g go-map.txt -h hierarchy.json -n rank-details.txt -s domain-details.txt
```

2. Output
* hierarchy_annotated.json
