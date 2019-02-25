## SAFE

For details on running SAFE, see

* [Systematic Functional Annotation and Visualization of Biological Networks](https://www.ncbi.nlm.nih.gov/pubmed/27237738)
* [Bitbucket](https://bitbucket.org/abarysh/safe/src)

### Make prey vs GO term matrix

Generate matrix of preys vs GO terms, with cells equal to the binary value (1 if prey has the term, 0 if not). This only needs to be done for the lowest correlation cutoff to test. Generate a separate (or merged) matrix for each GO namespace you want to test.

Requires:
* Gene Ontology (go-basic.obo)
* Gene annotations (goa_human.gaf without no header)
* list of genes
* GO namespace, default cc (one of all, bp, cc or mf)

1. Run script

```
"$CMSCRIPTS"/safe/go-table.pl -f go-basic.obo -g goa_human_nohead.gaf -l gene-list.txt -n cc
```

2. Output
* cc_annotations.txt

### Run SAFE

The SAFE iterator will run on each network and annotation matrix supplied, and test a range of enrichment radii.

Requires
* annotation matrices
* Cytoscape network

1. Open Matlab and specify path SAFE folder. Note, $CMSCRIPTS must be replaced by the path name.
```
addpath(genpath('"$CMSCRIPTS"/safe/'));
```

2. Create and move into a folder for running the analysis
```
mkdir "$CMSCRIPTS"/safe/analysis
cd "$CMSCRIPTS"/safe/analysis
```

3.	Copy safe.ini file from script folder to analysis folder
```
cp "$CMSCRIPTS"/safe/safe.ini .
```

4. Open `safe-iterator.m`

5. In script, set name of analysisFolder to full path of your folder
```
analysisFolder = '"$CMSCRIPTS"/safe/analysis';
```

6. Place all annotation matrices in `"$CMSCRIPTS"/safe/analysis/annotationMatrices`. Files should be named cc_something.txt, and Matlab will grab the part before the underscore for naming the results folder.

7. Place all networks to test in `"$CMSCRIPTS"/safe/analysis/networkFolder`. Networks should be named 0.7_something.cys and Matlab will grab the part before the underscore for naming the results folder.

8.	Specify an array of all radii to test
```
radii = [4 5 6 7 8 9 10]
```

9. Output
* Results, folder containing subfolder for each network + namespace combination

### Move files for assessment

Grab all attribute and node files, rename them to match their parent folders, delete the first four lines from each and move them into a folder called “safe_assessment” and a subfolder called “attributes” or “nodes” depending on their type.

1. Move to parent folder containing the “Results” folder generated in the previous step. Would be `"$CMSCRIPTS"/safe/analysis/` in this example.

2. Run script
```
"$CMSCRIPTS"/safe/safe-move-files.pl
```

### Assess SAFE terms

Take all terms within each domain (including children), and then for each node/gene in the network, it will see if the genes are being assigned to a domain with a previously known term.

Requires
* GO annotations, (goa_human.gaf without header)
* List of children for each GO term ()
* Map of GO IDs to terms
* GO namespace, default C, one of C, F or P

1. Move into folder `"$CMSCRIPTS"/safe/analysis/safe_assessment`

2. Run script
```
"$CMSCRIPTS"/safe/safe-assessor.pl -g go-children.txt -l goa_human_nohead.gaf -m go-map_cc.txt -n C
```

3. Output
* console

### After selecting best SAFE parameters

The below script assumes a SAFE parameter set and results have been selected. 

#### Perform GO enrichment on SAFE domains

The SAFE GO terms listed for each domain are not ordered and no enrichment scores are present. The following script will take all preys assigned to each domain and perform a g:Profiler enrichment to get an ordered list.

Requries
* node_properties_annotation-highest.txt
* GO namespace, one of GO:CC, GO:BP or GO:MF

1. Move into results folder for parameter set, e.g for a network with a correlation cutoff of 0.6, annotated with CC terms and a radius of 3.5
```
cd "$CMSCRIPTS"/safe/analysis/results/0.6cc_ccns_3.5r/
```

2. Run script
```
"$CMSCRIPTS"/safe/domainProfile.R node_properties_annotation-highest.txt GO:CC
```

2. Output
* `domain-summary.xlsx`
* `terms_perdomain.txt`: list of GO terms and pvalues for each domain

#### Assign representative term for domains

Use the domain-summary.xlsx from the previous step and the attribute_properties_annotation-highest.txt file to determine the representative term to use for each SAFE domain. See file in scripts folder called sample-safe-terms.txt on how to record this information into a file.

1. Select a representative term. If multiple terms, separate them with comma.

2. Also need to specify the “displayname” to use for the term on the cell map network, as well as the GO identifier, synonyms and the IC content value. Synonyms can be found in go-basic.obo file. 

3. Name the file domain-summary.txt.

#### Color Cytoscape network

Generate colors to use for SAFE nodes on cytoscape network for SAFE and color network.

Requires:
* node_properties_annotation-highest.txt (without the first four header lines)
* domain-summary.txt

1. Run script to generate node colours
```
"$CMSCRIPTS"/safe/safe-color.pl -n node_properties_annotation-highest.txt -s domain-summary.txt
```

2. Output
* domain-color.tsv

3.	Open cytoscape network and import domain-color.tsv as table for colouring nodes. Save as network_safe.cyjs for database (file will be used in u+q).


