## Upload related

### Prey summary

Generate summaries for preys in terms of how many baits they were seen with, the average spec (note: this is not control subtracted) and their NMF rank ratio.

Requires:
* NMF basis matrix
* SAINT file

1. Run Script
```
"$CMSCRIPTS"/upload-related/prey-summarize -b basis.csv -s saint.txt
```

2. Output
* `prey-summary.txt`: For each prey reports the number of baits it was seen with, the average spec across those baits and its NMF ratio
