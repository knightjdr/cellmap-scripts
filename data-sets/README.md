## Dataset formatting

### BioGRID

1. Download BioGRID dataset
* https://thebiogrid.org/download.php
* go to Current Release
* download file BIOGRID-ORGANISM-X.X.X.tab2.zip, unzip and grab the human file

2. Reformat BioGRID database so that it can be used in Cytoscape to display if an interaction edge is known
```
"$CMSCRIPTS"/data-sets/biogrid_forCytoscape.pl -b BIOGRID-ORGANISM-Homo_sapiens-X.X.X.tab2.txt
```

3. Output
* console

### IntAct

Requires:
* IntAct file
* organism, default 9606

1. Download IntAct dataset
* http://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip

2. Parse for human genes
```
"$CMSCRIPTS"/data-sets/intact_parse.pl -i intact.txt
```

3. Output:
* human-intact.txt

4. Merge BioGRID and IntAct datasets
```
"$CMSCRIPTS"/data-sets/merge-biogrid-intact.pl -b BIOGRID-ORGANISM-Homo_sapiens-X.X.X.tab2.txt -i human-intact.txt
```

5. Output:
* interactions.txt

### Gene Ontology

1. Download gene annotations (goa_human.gaf.gz) from http://current.geneontology.org/annotations/

2. Create a file called goa_human_nohead.gaf which simply removes the header lines from the above file.

3. Download ontology (go-basic.obo) from http://www.geneontology.org/page/download-ontology

#### Create a file without HPA annotations

Requires:
* `goa_human.gaf` without header lines

1. Run script
```
"$CMSCRIPTS"/data-sets/go-hpa.pl
```

2. Output
* `goa_human_nohpa.gaf`

#### Get child terms

Requires:
* go-basic.obo
* GO namespace, default cc (one of bp, cc or mf)

1. Run script
```
"$CMSCRIPTS"/data-sets/go-children.pl -h go-basic.obo -n cc
```

2. Output
* `go-terms.txt` list of all terms in specified namespace
* `go-children.txt` lists all children for each term

#### Map terms to IDs

Requires:
* `go-basic.obo`
* GO namespace, default C (one of C, F or P)

1. Run script
```
"$CMSCRIPTS"/data-sets/go-term-mapping.pl -h go-basic.obo -n C
```

2. Output
* go-map.txt

### UniProt

1. Download database
* http://www.uniprot.org/downloads
* Want the .txt format for the Reviewed (Swiss-Prot) release

### Information content (using ITGOM - oboslete, don't use this anymore)

This is obsolete since they haven't updated since 2014.

Papers:
* http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4256219
* http://www.hindawi.com/journals/bmri/2013/292063/

Requires
* list of GO IDs (go-terms.txt) from Gene Ontology step above

1. Go to http://web.cbio.uct.ac.za/ITGOM/tools/itgom.php.

2. Parameters:
* Calculating: “Term Information Content”
* Select Family: “Topology-based”
* Select Approach: Zhang et al.
* Select file go-terms.txt (The webtool can only take 2000 terms, so may need to split the list up into pieces to use)

3. Download table and save as za-information-content.txt.

### Protein complexes

1. Download protein complex file
* http://proteincomplexes.org/download
* Click `Protein Complex Map (genenames)`. Will produce file called `genename_clusters.txt` and rename to `humap.txt`.

### Corum

1. Download complex list
* http://mips.helmholtz-muenchen.de/corum/#download
* Download complete complexes in .txt format. Will produce file called `allComplexes.txt` and rename to `corum.txt`.

# Pfam

1. Download human domains
* ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/
* want the file 9606.tsv (need to delete first two lines)

## motifs

Pfam motifs (e.g. disordered and low complexity regions) can only be retrieved programmatically. For example, information for the Uniprot ID Q9NRA8 can be retrieved in JSON format at this URL: http://pfam.xfam.org/protein/Q9NRA8/graphic.

Requires:
* UniProt file (uniprot_sprot.txt)
* Species (specifed in script, Homo sapiens by default)

1. Run script
```
"$CMSCRIPTS"/data-sets/pfam-format.pl -u uniprot_sprot.txt
```

2. Output
* pfam_motifs.txt (will contain only two columns: gene name and motif). This file may have duplicate genes, so open in vim and do a `:sort u`
* pfam_motif_complete.txt will have all motifs, along with start and end positions, length, etc

## domains

Using UniProt as a map grab domains from Pfam

Requires:
* UniProt file (uniprot_sprot.txt)
* Pfam file (9606.tsv)

1. Run script
```
"$CMSCRIPTS"/data-sets/domain-extraction.pl -p 9606.tsv -u uniprot_sprot.txt
```

2. Output
* pfam-domains.txt
