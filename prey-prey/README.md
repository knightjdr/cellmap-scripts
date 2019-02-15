## Prey-Prey analysis

### Prey-Prey correlation

1. Go to https://prohits-viz.lunenfeld.ca

2. Choose correlation tool
* use BFDR for filtering
* abundance cutoff for prey correlation of 0
* correlation cutoff for cytoscape output: 0.4

3. Once complete, download results folder and get `preys_cytoscape.txt` and `preyvprey_df.tsv`

### Correlation filtering

Filter cytoscape file to generate input files with different correlation cutoffs. This will remove any preys not passing the cutoff as well as self links and duplicate interactions.

Requires
* prey_cytoscape.txt
* correlation cutoff

1. Run script
```
"$CMSCRIPTS"/prey-prey/filter-correlation.pl -c 0.4 -f preys_cytoscape.txt
```

2. Output
* preys_cytoscape_filtered_0.4.txt

### Generate gene list

Get list of preys passing correlation cutoff. This can be done just once for the lowest cutoff being tested and that file will work for all greater than it too.

Requires
* preys_cytoscape.txt

1. Run script
```
"$CMSCRIPTS"/prey-prey/get_unique.perl -l preys_cytoscape.txt
```

2. Output
* gene_list.txt

### Capturing known interactions

Determine what proportion of interactions are previously known, assuming a prey-prey passing a correlation cutoff is an interaction. For each correlation cutoff this script will report the fraction that are known.

Requires
* preyvprey_df.tsv from correlation tool (need to delete JSON bit)
* list of genes in correlation file or list of genes to check
* interaction file

1. Run script 
```
"$CMSCRIPTS"/prey-prey/prey_interactor_validation.pl -b interaction.txt -c preyvprey_df.tsv -g gene-list.txt -t m
```

2.	Output
* prey-interactors-recovered.txt

### Capturing known complexes

Determine what proportion of interactions are previously known to be part of the same complex, assuming a prey pair passing a correlation cutoff is an interaction. For each correlation cutoff this script will report the fraction that are known to be part of a complex. The script does not include self-interactions as it really skews the results.

Requires
* preyvprey_df.tsv from correlation tool (need to delete JSON bit)
* list of genes in correlation file or list of genes to check
* list of genes in protein complexes (from huMAP or Corum)
* file type: c for Corum (default) or h for huMAP

1. Run script
```
"$CMSCRIPTS"/prey-compex-validation.pl -b corum.txt -c preyvprey_df.tsv -g gene-list.txt -t c
```

2. Output

