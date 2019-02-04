## Cytoscape

Repeat all of the below steps for each correlation cutoff to test.

### Format correlation file for Cytoscape

Remove self-links and duplicate links in input file for Cytoscape.

Requires
* preys_cytoscape.txt (from correlation)

1. Run script
```
"$CMSCRIPTS"/cytoscape/format_preycorr.pl -c preys_cytoscape.txt
```

2. Output
* nodes.txt

### Generate network for SAFE

Requires
* gene list with all preys in network
* nodes.txt from previous step

1. Import as text gene_list.txt into Excel

2. Add “name” as column header

3. Duplicate this column and change new column’s header name to “ORF”

4. Save at gene_list-ORF.tsv

5. Launch Cytoscape

6. Import network from file: nodes.txt

7. Generate edge-weighted spring embedded layout for this file (use “none” for edge-weighting)

8. Import gene_list-ORF.tsv as table using “name” column as key

9. Move isolated network clusters as needed to make best use of space

10. Save session file .cys

11.	Also export as .cyjs file for coordinates (Select File-->Export --> Network and View... —> Select .cyjs format and name file)
