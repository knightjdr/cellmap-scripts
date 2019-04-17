## Organelle assessments

### ER

Take all ER baits and their preys and classifies the preys as
* cytoplasmic: if only seen with cytoplasmic baits
* lumenal: if only seen with lumenal baits
* transmembrane: if seen with at least one cytoplasmic and one lumenal

An optional list of genes can be included that will label matching preys with an *. This is helpful to see which of the lumenal preys were actually assigned to the lumenal compartment by NMF.

Requires:
* in the script, specify the lumenal and cytoplasmic baits
* SAINT file
* FDR cutoff
* list of genes to highlight (optional)

1. Run script
```
"$CMSCRIPTS"/organelle/er-assessment.pl -s saint.txt -f 0.01 -e list.txt
```

2. Output
* er-details.txt: for the cytosol, lumen and transmembrane it will list
  1. the number of baits the preys were seen with
  2. the number of preys
  3. the number of preys in the highlight list
  4. prey names (an asterisks indicates if it is in the highlight list)
