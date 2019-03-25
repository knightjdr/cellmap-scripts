## Upload related

### Crap

For the upload and query function, we want to allow users to ignore frequently flyers in BioID. From the cell map data set, we identify proteins seen frequently that we do not think will be informative by grabbing all significant preys, then defining crap as those seen with more than 20 baits (of 192). 

Requires:
* SAINT file
* FDR (default 0.01)
* BAIT threshold (default 20)

1. Run script
```
"$CMSCRIPTS"/upload-related/define-crap.R saint.txt
```

2. Output:
* `crap-list.txt`: list of crap proteins
* `crap-details.txt`: shows the number of occurrences and the median AvgSpec after control subtraction

### Enriched GO terms per bait

Determine enriched GO terms per bait for localization tab on bait report

Requires:
* SAINT file
* Map of bait names to database id
  * export from cell map using script in upload folder
  * tsv file with two columns (id, name)
* FDR (default: 0.01)

1. Run script
```
"$CMSCRIPTS"/upload-related/go-enrich.R saint.txt map.txt
```

2. Output
* `go-enrichment.json`
* validate the json at jsonlint.com (by copy paste) because of proxy errors with gProfiler

### Prey averages

Need average of the spectral count for every prey in the dataset, both with and without control subtraction. This is for the u+q report to display an extra column on the heat map for the dataset average.

Requires:
* Saint file

1. Run script
```
"$CMSCRIPTS"/upload-related/prey-average -s saint.txt
```

2. Output
* `prey_averages.txt`: contains prey name, average spectral count across baits, control subtracted average count and the best FDR for the prey. Spectral counts are rounded to 3 decimals.

### Prey summary

Generate summaries for preys in terms of how many baits they were seen with, the average spec (note: this is not control subtracted) and their NMF rank ratio.

Requires:
* NMF basis matrix
* SAINT file

1. Run Script
```
"$CMSCRIPTS"/upload-related/prey-summarize.pl -b basis.csv -s saint.txt
```

2. Output
* `prey-summary.txt`: For each prey reports the number of baits it was seen with, the average spec across those baits and its NMF ratio
