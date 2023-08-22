PAGAN: phylogenetic placement<a name="top"></a>
=============================

 

PAGAN extends existing alignments by aligning the new sequences to their phylogenetic positions in the reference alignment. This ensures that the new sequences are compared to the closest possible reference sequence and the alignment created, including the gaps opened for insertions and deletions, is as accurate as possible. Addition of new sequences does not affect the relative alignment of original sequences.

This page explains unguided placement. See [here](pagan_guided_placement.md) for guided placement.

* * *

*   [Placement of sequences to an existing alignment](#unguided-placement)
*   [Placement with strand search and ORF search](#placement-with-strand-search)
    *   [Example of unguided placement](#example-of-unguided-placement)
    *   [Example of ORF search and translated placement](#example-of-orf-search-and-translated-placement)

[Back to PAGAN front page](../README.md)

### Placement of sequences to an existing alignment

PAGAN can reconstruct ancestral sequence history for a reference alignment related by a phylogeny. Sequences can then be added to this reference alignment by aligning the new sequences against terminal or internal nodes and making room for the added sequences in their correct evolutionary positions among the reference sequences.

The minimal command to perform sequence placement is:

```
pagan --ref-seqfile ref_alignment_file --ref-treefile ref_tree_file --queryfile query_sequence_file
```

The reference alignment has to be in FASTA format and the reference tree in Newick or Newick eXtended (NHX) format; the query sequences that will be added in the reference alignment can be either in FASTA or FASTQ format. If the reference alignment consists of one sequence only, the reference guide tree is not required. Again, you can define your own output file with option ```--outfile```.

Some placement functionalities have multiple alternative options.  
 

By default, PAGAN places the query sequences anywhere in the reference phylogeny. This can be changed with options:  
 

*   default: place queries anywhere
*   ```--terminal-nodes```: place queries to terminal nodes
*   ```--internal-nodes```: place queries to internal nodes

By default, the placement is done one query sequence add time and after each placement the reference tree and alignment are updated; the next sequence can then be placed also next to the newly added sequence. This strategy is not ideal for short, incomplete sequences and the use of option ```--fragments``` is recommended. This option finds the target for all queries first and then adds them one target node at time, adding the query sequence always to the ancestor node of the preceding query sequence. This may result in an unrealistic, ladder-like tree but the strategy avoids alignment errors when the query sequences overlap only partially.

PAGAN uses Exonerate to speed up the target selection. By default, it first pre-selects a smaller set of targets and then does the final selection among these. (This makes a difference only with very large reference sets.) From these, PAGAN first finds the best targets with Exonerate ungapped alignment, then from these with Exonerate gapped alignment, and finally from those with PAGAN alignment. The number of sequences retained at each stage can be selected with a couple of shorthand options, or by defining the numbers precisely:

*   default: 10 / 3 / 1+ (retained: ungapped/gapped/PAGAN)
*   ```--fast-placement```: 5 / 1 / 1
*   ```--very-fast-placement```: 1 / 1 / 1
*   ```--own-placement```: – / – / 1+

By default, PAGAN places a query sequence to multiple locations if more than one location shares the best score. This can be changed with option ```--one-placement``` that places the query only once, at the deepest (or last) node among the possible placements.

### Placement with strand search and ORF search

If the correct strands of the reads are unknown, PAGAN can be run with option \--both-strands that performs both forward and reverse-complement alignment and chooses the better one:

```
pagan --ref-seqfile ref_alignment_file --ref-treefile ref_tree_file --queryfile query_sequence_file --both-strands
```

If the query sequences come from protein-coding sequences (e.g. contigs from an RNA-seq experiment), PAGAN can be run with option ```--find-orfs``` that searches for open reading frames in the query (both in forward and reverse-complement strands) and chooses the one giving the best alignment:

```
pagan --ref-seqfile ref_alignment_file --ref-treefile ref_tree_file --queryfile query_sequence_file --translate --find-orfs --min-orf-length 30
```

The reference alignment should be DNA (and translated with option ```--translate```) or the resulting alignment cannot be back-translated to DNA. In the example above, the minimum ORF length is set to be 30 aa (default 50); the length can also be defined as proportion of the full length. The output is written separately for translated (peptide) and back-translated (DNA) data.

It should be noticed that the ORF search does not correct for reading frames and may thus not work well with data coming e.g. from 454.

[back to top](#top)

#### Example of unguided placement

The following files are included in directory [examples/ngs_placement](../examples/ngs_placement):

*   ```reference_codon.fas``` : simulated reference alignment.
*   ```reference_tree.nwk``` : reference phylogeny (no tags).
*   ```input_ngs.fastq``` : simulated Illumina reads.
*   ```input_ngs_primates.fastq``` : a subset of the one above.

The placement of short query sequences is approximate and is only performed to allow for the accurate alignment of the data (e.g. for downstream analyses to build longer contigs of the overlapping sequences). If the aim is to infer the phylogenetic origin of the query sequences, one should re-analyse the alignments with an appropriate program.

Placement can be performed using e.g. the following command:

```
pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nwk --queryfile input_ngs_primates.fastq --outfile read_alignment --fast-placement --fragments
```

PAGAN aligns and outputs full sequences but it may sometimes discard complete reads if no good placement is found (by Exonerate) or the resulting alignments (by PAGAN) looks poor. A placement for some reads discarded by Exonerate may be found with PAGAN’s exhaustive search that is selected with option ```--exhaustive-placement```; the threshold for the read rejection at the later stages can be adjusted with options ```--overlap-minimum```, ```--overlap-identity``` and ```--overlap-identical-minimum```. See ```pagan --help``` for the default values.

One should note that PAGAN can place one query sequence to multiple locations if they score equally well. This makes sense in some analyses but may cause problems in others. This behaviour can be disabled with option ```--one-placement-only```.

#### Example of ORF search and translated placement

If the query sequences are rather distantly related to the reference sequences, a better alignment may be obtained by translating the sequences and then performing the alignment as protein sequences. If all the input are DNA, the resulting extended alignment can be back-translated to DNA. Translated alignment also maintains the codon structure and the reading frame.

Translated placement can be performed using e.g. the following command:

```
pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nwk --queryfile input_ngs_primates.fastq --outfile read_alignment --fast-placement --fragments --translate --find-orfs --min-orf-length 30
```

By replacing ```--translate --find-orfs``` with ```--both-strands```, one can extend alignment with query data where some sequences may be in reverse-complement orientation.

[back to top](#top)