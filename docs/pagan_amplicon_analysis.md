PAGAN: amplicon analysis<a name="top"></a>
========================

   
Amplicon analysis is not much different from standard phylogenetic placement, see [here](pagan_phylogenetic_placement.md) for instructions on that. However, in amplicon analysis query sequence are generated with PCR and are thus expected to match a specific target region in the reference alignment. Typically only the region that the query sequences match and only the reference sequences closely-related to the query sequences are of interest. To facilitate amplicon analysis, PAGAN contains options to trim and prune the resulting alignment.

* * *

[Back to PAGAN front page](../README.md)

The extra options for amplicon analysis are bets explained with an example. The data used in the examples below can be found [here](../examples/amplicon_analysis).

The following files are included:

*   ```pagan_placement.cfg``` : config file that runs a basic analysis;
*   ```query.fas``` : amplicon sequences to place in the alignment;
*   ```reference.fas``` : reference alignment (these from Silva);
*   ```reference.tree``` : reference phylogeny.

The data can be analysed using the following command:  
 
```
pagan pagan_placement.cfg
```

The same analysis could be performed with the command:

```
pagan --ref-seqfile reference.fas --ref-treefile reference.tree \
 --queryfile query.fas --outfile pagan_placement --xml \ 
 --one-placement-only --trim-extended-alignment \
 --prune-extended-alignment --prune-keep-number 0 \
 --prune-keep-closest 
```

This produces several output files:  
 

*   ```pagan_placement.???``` contain the full data;
*   ```pagan_placement.trimmed.???``` contain the target region for all the sequences;
*   ```pagan_placement.pruned.???``` contain the target region for the query sequences and for N reference sequences (here N=0);
*   ```pagan_placement.pruned_closest.???``` contain the target region for  the query sequences and the very closest reference sequences.

Files ending with ```.tre``` are the alignment guide trees for each output; those ending 
with ```.fas``` and ```.xml``` are the alignments in FASTA and HSAML format.  
 

If the amplicon sequences are very closely related to the reference sequences, one can  
add option ```--terminal-nodes``` and thus enforce their placement next to the tips.

Options related to trimming and pruning are:  
 

*   ```--prune-extended-alignment``` : prune output;
*   ```--prune-keep-number arg``` : prune output and keep N most distantly related sequences (if = 0, keep only query sequences);
*   ```--prune-keep-threshold arg``` : prune output and remove sequences with distance below threshold;

* ```--prune-keep-closest``` : prune output and keep only closest reference sequences;

*   ```--trim-extended-alignment``` : trim reference alignment outside query target region;
*   ```--trim-keep-sites arg``` : trim distance around query sequences.

The ```.xml``` output files can be visualised with Wasabi.