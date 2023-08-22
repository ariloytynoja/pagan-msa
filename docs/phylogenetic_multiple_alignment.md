# PAGAN: phylogenetic multiple alignment <a name="phylogenetic-multiple-alignment"></a>


*   [Phylogenetic multiple alignment](#phylogenetic-multiple-alignment)
    *   [Input and output formats](#input-and-output-formats)
    *   [Codon alignment and translated alignment](#codon-alignment-and-translated-alignment)
    *   [Alignment anchoring](#alignment-anchoring)
    *   [An example analysis with PAGAN](#an-example-analysis-with-pagan)
*   [Inference of ancestral sequences for existing alignments](#inference-of-ancestral-sequences-for-existing-alignments)

[Back to PAGAN front page](../README.md)

PAGAN is based on a progressive algorithm that aligns sequences according to a guide tree. It can compute a guide tree using external programs (MAFFT and BppDist or RAxML) but also allows the user to provide a **rooted** binary tree relating the sequences. The leaf names in the tree and the sequence names (until the first space) in the sequence file have to match exactly. Alignment is only performed for the parts of the guide phylogeny that have sequences associated; the unnecessary branches and sequences are pruned/dropped out.

The PAGAN binary packages contain helper applications to compute guide trees, anchor alignments and reconstruct ancestral character states. When using one of those, the minimal command to perform the alignment is:

```
./pagan --seqfile sequence_file 
```

Using a user-specified tree, the minimal command to perform the alignment is:

```
./pagan --seqfile sequence_file --treefile tree_file 
```

The resulting alignment will be written in file outfile.fas. If you want to use another file name, you can specify that with option ```--outfile```:

```
./pagan --seqfile sequence_file [--treefile tree_file] --outfile another_name
```

PAGAN will automatically add a suffix for the corresponding output format (by default, ```.fas``` for FASTA).

#### Input and output formats

The sequence input file has to be in FASTA format and the guide tree in Newick tree format, with branch lengths as substitutions per site. By default, the resulting alignment will be written in FASTA format. Other output formats (nexus, paml, phylipi, phylips or raxml) can specified with option ```--outformat format_name```. With option ```--xml```, PAGAN writes an additional file, containing both the alignment and the guide phylogeny, in the XML-based [HSAML format](https://github.com/ariloytynoja/prank-msa/blob/master/docs/hsaml_format.md).

PAGAN writes several temporary files to communicate with the helper tools. By default, these files are written to ```/tmp``` directory (or current directory if ```/tmp``` does not exist) and are deleted afterwards. The directory can be changed with option ```--temp-folder directory_path```.

#### Codon alignment and translated alignment

PAGAN supports the alignment of nucleotide, amino-acid and codon sequences as well as translated alignment of nucleotide sequences. The type of data is detected automatically and either DNA or protein model is used. With option ```--codons```, PAGAN can align protein-coding DNA sequences using the codon substitution model. Alternatively, with option ```--translate```, PAGAN translates the DNA sequences to proteins, aligns them as proteins and writes the resulting alignment as proteins and back-translated DNA. Codon and translated alignment are performed in the first reading frame and no frame correction is performed.

#### Alignment anchoring

PAGAN can significantly speed up the alignment by anchoring the alignments with program Exonerate or by finding long, shared prefixes. Anchoring with Exonerate is used by default while prefix anchoring can chosen with option ```--use-prefix-anchors```. Anchoring can be disabled with option ```--no-anchors```; see ```pagan --help``` for related options.

#### An example analysis with PAGAN

As an example of combining several options, the following command:

```
./pagan --seqfile cDNA_unaligned.fas --outfile cDNA_aligned  --outformat paml --codons  --config-log-file cDNA_alignment.cfg
```

would 
 1. read in protein-coding DNA sequences from the file ```cDNA_unaligned.fas``` 
 2. translate the sequences to proteins and perform an initial alignment with MAFFT
 3. infer a guidetree for this protein alignment with BppDist
 4. treating the sequences as codon sequences, perform the PAGAN alignment using anchors inferred with Exonerate on protein-translated sequences
 5. write the output alignment to the file ```cDNA_aligned.phy``` in a format compatible with PAML and the BppDist-estimated guidetree to the file ```cDNA_aligned.tre``` in Newick format. 
 
The last option tells PAGAN to write the options used to the file ```cDNA_alignment.cfg```: this file can be stored as a record of the analysis performed but it also fully specifies the options for PAGAN and allows repeating the same analysis with the command:

```
./pagan cDNA_alignment.cfg
```

Parameters can overridden and added by specifying them after the config file. For example, a similar analysis as above but with a RAxML-estimated guidetree and output in the Nexus format would be performed with command:

```
./pagan cDNA_alignment.cfg --raxml-tree --outformat nexus
```

[back to top](#phylogenetic-multiple-alignment)

* * *
<a name="ancestral-sequences"></a>
### Inference of ancestral sequences for existing alignments 


In phylogenetic placement, PAGAN infers evolutionary history for a set of aligned sequences. This functionality can be used without placement, too. Command:

```
./pagan --ref-seqfile alignment_file --ref-treefile tree_file --xml
```

writes the aligned sequences from ```alignment_file``` to files ```outfile.fas``` and ```outfile.xml```, thus providing a tool to convert FASTA files to HSAML format; the HSAML — and the inferred ancestral sequences included — can be visualised with Wasabi. The ancestral sequences can also be outputted in FASTA format using option ```--output-ancestors```. Command:

```
./pagan --ref-seqfile alignment_file --ref-treefile tree_file --output-ancestors
```

includes the ML-reconstructed internal nodes in the FASTA output, providing a efficient tool to infer gap structure for ancestral sequences based on existing alignments. The tree indicating the ancestral nodes is written in outfile.anctree.