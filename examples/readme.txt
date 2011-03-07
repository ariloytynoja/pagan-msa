
This directory contains example data for NGS read placement using PAGAN. 

The following files are included:

- reference_alignment.fas : a simulated reference alignment.
- reference_tree.nhx : a phylogeny relating the reference sequences with one node tagged for 
  read placement.
- illumina_reads.fastq : simulated Illumina reads.

The data can be analysed using the following command:

pagan --ref-seqfile reference_alignment.fas --ref-treefile reference_tree.nhx 
 --readsfile illumina_reads.fastq --outfile read_alignment

The resulting alignments will be written to files read_alignment.fas and read_alignment.xml 
in FASTA and HSAML formats.

If we add the option '--config-log-file':

pagan --ref-seqfile reference_alignment.fas --ref-treefile reference_tree.nhx 
 --readsfile illumina_reads.fastq --outfile read_alignment --config-log-file pagan_alignment.cfg

the arguments used for the analysis will be written to the file 'pagan_alignment.cfg'. 
This file works as a reminder of the specific settings used for the analysis and allows 
repeating the same analysis with the command:

pagan --config-file pagan_alignment.cfg

