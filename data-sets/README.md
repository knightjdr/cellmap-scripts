## Dataset formatting

### BioGRID

1. Download BioGRID dataset
* https://thebiogrid.org/download.php)
* go to Current Release
* download file BIOGRID-ORGANISM-X.X.X.tab2.zip, unzip and grab the human file

### IntAct

1. Download IntAct dataset
* http://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip

2. Parse for human genes
```
$HOME/cellmap-scripts/data-sets/intact_parse.pl -i intact.txt
```

3. Output:
* parsed-intact.txt

4. Merge BioGRID and IntAct datasets
```
$HOME/cellmap-scripts/data-sets/merge-biogrid-intact.pl -b BIOGRID-ORGANISM-Homo_sapiens-X.X.X.tab2.txt -i parsed-intact.txt
```

5. Output:
* interactions.txt