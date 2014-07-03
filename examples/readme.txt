
The three sub-directories contain example data for guided sequence placement with PAGAN
(protein _placement and ngs_placement) and for NGS read pile-up (454_pileup).

==================
UNGUIDED PLACEMENT
==================

The idea of unguided placement is to extend existing alignments with new sequences. The
optimal location for the new sequences is searched with different heuristics.

The examples are based on simulated data. In these examples, we have reference 
alignments that include gene sequences from several mammals; the hypothetical gene has 
duplicated first in the ancestor of primates and rodents and then again in the ancestor
of rodents. Our new sequences come from a primitive primate, believed to belong 
phylogenetically between tarsier and lemur, and a rodent, belonging between guinea pig
and ground squirrel. 

When sequences are placed with PAGAN, it first tries first finds the best location for
the query sequence and then aligns it there. If sevaral target nodes score equally well,
the sequence is, by default, assigned to each one of them. Importantly, the relative
alignment of existing reference sequences is not changed. 

The examples below demonstrate that the selection of the target node for each query
sequence is not straightforward. Long sequences are placed well but shorter ones
may be placed to wrong nodes. 


================
GUIDED PLACEMENT
================

The idea of guided placement is to extend existing alignments with new sequences for 
which the (rough) phylogenetic location is known. The situation may be complicated by 
paralogous genes and uncertainty of the right target for the new sequences. For such 
cases, PAGAN allows indicating several potential targets and then chooses the best 
location by aligning the query sequences to each specified target.

The examples are based on simulated data. In these examples, we have reference 
alignments that include gene sequences from several mammals; the hypothetical gene has 
duplicated first in the ancestor of primates and rodents and then again in the ancestor
of rodents. Our new sequences come from a primitive primate, believed to belong 
phylogenetically between tarsier and lemur, and a rodent, belonging between guinea pig
and ground squirrel. Because of the gene duplications, the primate and rodent sequences 
have two and three potential target locations, respectively, in the phylogeny. These
are indicated in the reference phylogenies using special TID tags. The phylogenies and
the target locations for guided placement can be viewed and edited with the 
ArchaeopteryxPE software available on the PAGAN download site.

When sequences are placed with PAGAN, it first tries placing each query sequence to 
every matching target node and chooses the best node. If sevaral target nodes score 
equally well, the sequence is assigned to each one of them. PAGAN then starts aligning 
the sequences to targets, adding several query sequences to the same target using a 
progressive approach. The complete new subtree is then inserted back to the reference
alignment, adding space for the new insertions if necessary. Importantly, the relative
alignment of existing reference sequences is not changed. 

The examples below demonstrate that the selection of the target node for each query
sequence is not straightforward. Long sequences are placed correctly but shorter ones
may be placed to a wrong node. Errors are also more likely between closely-related 
paralogs than between distant ones.


=================
AMPLICON ANALYSIS
=================

Amplicon analysis is like standard phylogenetic placement except that it is helpful to
trim and prune the resulting alignment such that it contains only the target reagion
and a smaller number (or no) reference sequences. Running the example should explain
that nicely and more infromation is provided at the program homepage.
 

================
PILEUP ALIGNMENT
================

PAGAN can make "pileup" alignments by adding the sequences in the order their appear 
in the input file. This is not recommended for distantly-related sequences but can be 
useful for highly similar sequences that require alignment. This may be relevante.g. 
in the analysis of overlapping noisy reads from the same locus; if the reads have been 
generated on Roche 454 or Ion Torrent platform, PAGAN's ability to model the 
homopolymer errors is especially useful.

Like in any alignment also in a pileup alignment the consenus sequence (the alignment 
of sequences included so far) and the next sequence to be added should overlap. PAGAN 
can include sequences that do not overlap but the region in between has to be bridged 
by other sequences and thus the order of adding the sequences can be important. If the
sequences are from the same species, the option to reconstruct the consensus sequence
(see below) can also be useful.



==================
protein_placement:
==================

The following files are included:

- reference_aa.fas   : simulated reference alignment.
- reference_tree.nhx : reference phylogeny with some nodes tagged for placement.
- reference_tree.nwk : reference phylogeny without tags.
- input_aa_full.fas  : simulated amino-acid sequences.
- input_aa_frags.fas : a subset of the one above, broken into fragments.

The data can be analysed using the following commands:

UNGUIDED PLACEMENT:

pagan --ref-seqfile reference_aa.fas --ref-treefile reference_tree.nwk \
 --queryfile input_aa_full.fas --outfile aa_full_alignment

pagan --ref-seqfile reference_aa.fas --ref-treefile reference_tree.nwk \
 --queryfile input_aa_frags.fas --outfile aa_frags_alignment --fragments


GUIDED PLACEMENT:

pagan --ref-seqfile reference_aa.fas --ref-treefile reference_tree.nhx \
 --queryfile input_aa_full.fas --outfile aa_full_alignment --guided

pagan --ref-seqfile reference_aa.fas --ref-treefile reference_tree.nhx \
 --queryfile input_aa_frags.fas --outfile aa_frags_alignment --guided \
 --fragments

The resulting alignments will be written to files aa_[full|frags]_alignment.fas.



==============
ngs_placement:
==============

The following files are included:

- reference_codon.fas : simulated reference alignment.
- reference_tree.nhx  : reference phylogeny with some nodes tagged for placement.
- reference_tree.nwk  : reference phylogeny without tags.
- input_ngs.fastq     : simulated Illumina reads.
- input_ngs_primates.fastq : a subset of the one above.

The data can be analysed using the following command:

UNGUIDED PLACEMENT:

pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nwk \
 --queryfile input_ngs.fastq --outfile read_alignment --fast --fragments


GUIDED PLACEMENT:

pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nhx \
 --queryfile input_ngs.fastq --outfile read_alignment --fast --guided \
 --fragments

The resulting alignments will be written to file read_alignment.fas.

If we add the option '--config-log-file':

pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nwk \
 --queryfile input_ngs.fastq --outfile read_alignment --fast --fragments \
 --config-log-file simple.cfg

the arguments used for the analysis will be written to the file 'simple.cfg'. This file 
works as a log file of the specific settings used for the analysis and allows repeating 
the same analysis with the command:

pagan simple.cfg  (or 'pagan --config-file simple.cfg')

The use of config files simplifies complex commands. For example, the following:
 
pagan --ref-seqfile reference_codon.fas --ref-treefile reference_tree.nwk \
 --queryfile input_ngs.fastq --build-contigs --use-consensus --show-contig-ancestor \
 --consensus-minimum 3 --outfile read_alignment --fast --fragments \
 --config-log-file contigs.cfg

specifies config file 'contigs.cfg'. As the arguments given on the command line 
override those given in a config file, the following command repeats a similar 
analysis for another data set.

pagan contigs.cfg --queryfile input_ngs_primates.fastq --outfile primate_alignment

You may check the content of files 'read_alignment_contigs.fas' and 
'primate_alignment_contigs.fas'. See 'pagan --help' and the program documentation 
for additional options.



=================
amplicon analysis
=================

The following files are included:

- pagan_placement.cfg : config file that runs a basic analysis.
- query.fas           : amplicon sequences to place in the alignment.
- reference.fas       : reference alignment (these from Silva).
- reference.tree      : reference phylogeny.

The data can be analysed using the following command:

pagan pagan_placement.cfg

This produces several output files. 

Files named "pagan_placement.???" contain the full data; "pagan_placement.trimmed.???" 
contain the target region for all the sequences; "pagan_placement.pruned.???" contain 
the target region for the query sequences and for N reference sequences (here N=0); and 
finally those named "pagan_placement.pruned_closest.???" contain the target region for 
the query sequences and the very closest reference sequences. 

Files ending with ".tre" are the alignment guide trees for each output; those ending 
with ".fas" and ".xml" are the alignments in FASTA and HSAML format. 
  
The same analysis could be performed with the command:

pagan --ref-seqfile reference.fas --ref-treefile reference.tree --queryfile query.fas \
 --outfile pagan_placement --one-placement-only --trim-extended-alignment \
 --prune-keep-number 0 --prune-keep-closest --xml

If the amplicon sequences are very closely related to the reference sequences, one can 
add option '--terminal-nodes' and thus enforce their placement next to the tips.



===========
454_pileup:
===========

The following files are included:

- 454_reads.fas : simulated overlapping 454 reads
- 454_reads_reversed.fas : simulated overlapping 454 reads, half of them reversed

The data can be analysed using the following command:

pagan --pileup-alignment --use-consensus --454 --queryfile 454_reads.fas \
 --outfile 454_reads_pagan --config-log-file 454.cfg

The resulting alignment will be written to file 454_reads_pagan.fas and the the 
arguments used to '454.cfg'.

A similar analysis could now be repeated for another input file with the command:

pagan 454.cfg --queryfile more_reads.fas --outfile another_output_pagan



==================
strand/ORF_search:
==================

If the strand of the reads is unknown, PAGAN can be run with an option that performs
both forward and reverse-complement alignment and chooses the better one:

pagan --pileup-alignment --use-consensus --454 --queryfile 454_reads_reversed.fas \
 --outfile 454_reads_reversed_pagan --both-strands

If the strand of the reads is unknown and the reads come from protein-coding sequences
(e.g. from RNA-seq experiment), PAGAN can be run with an option that searches for open
reading frames in the query (both in forward and reverse-complement strands) and 
chooses the one giving the best alignment:

pagan --pileup-alignment --ref-seqfile human.fas --queryfile input_ngs_primates.fas \
 --translate --find-orfs

The reference sequence should also be DNA or the resulting alignment cannot be 
back-translated. This option does not correct for reading frames and may not work well
with data coming e.g. from 454.

