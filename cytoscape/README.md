## Cytoscape

* Cytoscape files for different correlation cutoffs are generated during the prey-prey step
* use those files for generating networks

### Generate network for SAFE

Requires
* gene list with all preys in network
* cytoscape file

1. Import gene list as text into Excel

2. Add “name” as column header

3. Duplicate this column and change new column’s header name to “ORF”

4. Save as names-orf.txt

5. Launch Cytoscape

6. Import network from cytoscape text file

7. Generate edge-weighted spring embedded layout for this file (use “none” for edge-weighting)

8. Import names-orf.txt as table using “name” column as key

9. Move isolated network clusters as needed to make best use of space

10. Save session file .cys

11.	Also export as .cyjs file for coordinates (Select File-->Export --> Network and View... —> Select .cyjs format and name file)
