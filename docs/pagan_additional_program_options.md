PAGAN: additional program options<a name="top"></a>
=================================

 
This part of the documentation has not been completed.  
  
 

[Back to PAGAN front page](../README.md)

Some of the options printed by ```./pagan --help``` relate to unfinished features and may not function properly.

#### Generic options[](#Generic_options)

*   Option ```--silent``` minimises the output (doesn’t quite make it silent, though).

*   Option ```--xml``` writes the output alignment also in HSAML format that can be imported to Wasabi.

*   Options ```--indel-rate``` can be used adjust the insertion and deletion rates. Although it would be possible to consider the two processes separately, this has not been implemented yet.

*   Options ```--gap-extension``` and ```--end-gap-extension``` define the gap extension probability for regular and terminal gaps. For meaningful results, the latter should be greater (and, for pair-end data, equal to ```--pair-read-gap-extension```).

*   Options ```--dna-kappa``` and ```--dna-rho``` affect the DNA substitution scoring matrix; base frequencies are estimated from the data.

*   Options ```--codons``` translates DNA sequences to codons and aligns them using the codon substitution model. (experimental)

*   Options ```--scale-branches```, ```--truncate-branches``` and ```--fixed-branches``` override the branch lengths defined in the guide tree. By default, long branches are truncated to make the scoring matrix more informative; this can be prohibited with ```--real-branches```.

*   Option ```--output-ancestors``` writes the parsimony-reconstructed ancestral sequences for the internal nodes of the tree. The tree indicating the nodes is written in ```outfile.anctree```.

*   Option ```--config-file file_name``` specifies a config file. If an option is specified both in a config file and as a command-line argument, the latter one overrides the former.

*   Option ```--config-log-file file_name``` specifies a log file where (non-default) options used for the analysis are written. The format is compatible with the option input and the file can be used a config file.

There are many parameters related to "insertion calling", the type and amount of phylogenetic information required to consider insertion-deletion as an insertion and thus prevent the later matching of those sites. These parameters are still experimental (although some of them are used and affect the resulting alignment) and will be described in detail later.

#### Phylogenetic placement options[](#Phylogenetic_placement_options)

Some options are only relevant for the placement of sequences into existing alignment.

*   Option ```--query-distance``` sets the expected distance between the query and the pseudo-parent node (against which the query is aligned) and thus affects the substitution scoring used in the alignment. Having the distance very short (default), the alignment is stringent and expects high similarity.

*   Queries with too few sites aligned against sites of reference sequences are discarded. (The stringency of the alignment is set using the option above.) Options ```--min-query-overlap``` and ```--min-query-identity``` set the required overlap and base identity for accepting the query.

Many options are either not important for basic use or self-explanatory (or both).

[back to top](#top)