PAGAN: New features
===================

[Back to PAGAN home.](../README.md)

### New features in the latest version of PAGAN

#### v.1.51

First release of PAGAN2. Now using built-in BLAST library (instead of Exonerate)  
and doing efficient alignment achoring that massively reduces the memory usage.  
Can now align closely related sequences that are up to megabases long on a regular  
laptop computer.

#### v.0.61

Fixed a serious bug in query placement of independent sequences: previously  
substrings (insertions) may have been left out of the resulting alignment.  
Also fixed the branch length computation for the new added sequences.  
PAGAN code is now hosted at github: changed the update check to reflect that.

#### v.0.60

Using now explicit data types for Exonerate. Previously DNA sequences with  
lots of ambiguous characters could be handled as proteins and search failed.

#### v.0.59

Fixing an issue with multiple placement of the same query. Adding a new option  
“–output-discarded-queries” for one-by-one placement. Enabling placement to newly  
added taxa now also for translated data. Speed-up for Exonerate search.

#### v.0.58

Adding ancestral sequence reconstruction for codon data when “–translate” is used.  
A backup system for cases where BppAncestors fails.

#### v.0.57

Faster reference build (mainly for aa/codons). Less verbose output. Exonerate hits  
trimmed to avoid overlapping hits.

#### v.0.56

By mistake, PAGAN pruned the files even if no option was selected. This is now fixed.

#### v.0.55

The program homepage has been re-written and PAGAN is now compatible with that. Several changes in placement heuristics and tree extension.

#### v.0.54

Cleaning up the code: many unused options removed. Xml output with nhx tags. Changed search order of external tools. Now using BppAncestors for ancestor reconstruction.  
 

#### v.0.53

TID tags now considered also for ORFs. Fixing the double naming of ORFs. Changed a bad default value for branch-skip-probability (not fixing the bug yet!).  
 

#### v.0.52

Consistent naming of nucleotide and protein output of ORF alignments. Improvements for alignment pruning. Option to trim extended alignment.  
 

#### v.0.51

Functions for alignment pruning (to reduce large datasets). Anchoring coverage threshold for skipping the DP alignment altogether. Options ‘find-best-orf’ and ‘compare-reverse’ also with placement.  
 

#### v.0.50

Improvements for pileup of 454 data. Removal of gaps in query sequences before alignment.  
 

#### v.0.49

New option to create a distance tree using ‘bpp\_dist’ package. Significant changes (and improvements) for pileup alignment of 454 reads.  
 

#### v.0.48

New option’test-every-terminal-node’.  
 

#### v.0.47

Compiles on Darwin and Cygwin.  
 

#### v.0.46

First integration of MAFFT and RAxML: PAGAN is now capable of computing a guide tree by first calling MAFFT and then RAxML.

#### v.0.45

Minor: code reorganised and recently created bugs fixed.

#### v.0.44

PAGAN can now anchor alignments with Exonerate or substring matching. As a result, significant speed ups in the alignment of long and/or very similar similar sequences.

#### v.0.43

Now properly functional multi-threading for standard multiple alignment: by default, PAGAN uses all the cores available. This can be changed with \--threads #.

#### v.0.42

First attempt to introduce multi-threading for standard multiple alignment.

#### prior v.0.42

See VERSION\_HISTORY in the source code