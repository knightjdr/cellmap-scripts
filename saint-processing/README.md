## SAINT file processing

### Merge SAINT files

Requires:
* SAINT files
* bait list

1. Place all SAINT files (ending in .txt) into folder with no other files.
2. Create a file with a list of baits. Only needs one column but can have more.

| bait  | other |
|-------|-------|
| bait1 | a     |
| bait2 | b     |
| ...   | ...   |

3. Run saint_merge.pl in the directory with the SAINT files and specifiy path to list of baits. In this example, the assumption is the bait list is one directory above the folder with the SAINT files.
```
"$CMSCRIPTS"/saint-processing/saint-merge.pl -b ../baits.txt
```

4. Outputs:
* merged.txt conating the merged SAINT files
* missing_merged.txt with a list of any baits in the bait list file that were not found
* skipped_merged.txt to report any baits in the SAINT files that are not being included

### Filter SAINT file

Requires:
* SAINT file
* contaminants list
* FDR cutoff, default 0.01
* minimum prey threshold, default 5

Contaminant list is in this repo `data-files/contaminants.txt`. It was generated from:
* ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta 
* http://www.coxdocs.org/doku.php?id=maxquant:common:download_and_installation and grabbing contaminants.fasta
* supplemented with BirA tag, streptavidin, mCherry, GFP, bovine, casein and adenoviral proteins

Script will filter contaminants, decoys and any baits that do not have at least 5 preys passing the cutoff of 0.01 FDR. These thresholds can be adjusted with the -m and -f flags respectively.

1. Run script using contaminants.txt and SAINT file.
```
"$CMSCRIPTS"/saint-processing/saint-filter.pl -c contaminants.txt -s saint.txt
```

2. Output:
* saint_filtered.txt
* filtered_baits.txt to list any baits removed

Note: Need to manually remove Human prey genes found with mouse baits