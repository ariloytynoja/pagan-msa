**See [https://ariloytynoja.github.io/pagan-msa](https://ariloytynoja.github.io/pagan-msa).**


### Note!

This is the original version of the program and the **much improved PAGAN2** can be found at [https://github.com/ariloytynoja/pagan2-msa](https://github.com/ariloytynoja/pagan2-msa).  

---

PAGAN
=====


PAGAN is a general-purpose method for the alignment of sequence graphs. PAGAN is based on the phylogeny-aware progressive alignment algorithm and uses graphs to describe the uncertainty in the presence of characters at certain sequence positions. However, graphs also allow describing the uncertainty in input sequences and modelling e.g. homopolymer errors in Roche 454 reads, or representing inferred ancestral sequences against which other sequences can then be aligned. PAGAN is still under development and will hopefully evolve to an easy-to-use, general-purpose method for phylogenetic sequence alignment.

As the graph representation has features that make PAGAN especially powerful for phylogenetic placement of sequences into existing alignments, the functionality necessary for that was implemented first. The method and its uses for alignment extension are described in  [http://bioinformatics.oxfordjournals.org/content/28/13/1684.full](http://bioinformatics.oxfordjournals.org/content/28/13/1684.full).


* * *

![](docs/data/pagan_logo.png)

*   **[Download the original PAGAN](binaries/)**

Pre-compiled binaries for Linux (and OSX).

*   **[Download the fancy new PAGAN2](https://github.com/ariloytynoja/pagan2-msa)**

Source and pre-compiled binary for Linux.

*   **[Installing PAGAN](docs/pagan_installation.md)**

Instructions for installation and using of PAGAN on different operating systems.

*   **[Using PAGAN](#using-pagan)**

Instructions for the use of PAGAN with example commands and test data.

*   **[Building PAGAN2 with Docker](https://github.com/ariloytynoja/pagan2-msa)**

Instructions for building PAGAN2 with Docker and using the PAGAN2 Docker container.

*   **[New features](docs/pagan_new_features.md)**

Brief description of the new features in the latest versions of the program


* * *

### Using PAGAN

At the simplest, PAGAN can be run with command:

```
pagan --seqfile input_file
```

where ```input_file``` contains sequences in FASTA format.

PAGAN2 has **much** improved anchoring and memory handling and often this runs much faster:

```
pagan2 --seqfile input_file
```

See the links below for the description of the most central program options.

*   [Phylogenetic multiple alignment](docs/phylogenetic_multiple_alignment.md)
*   [Inference of ancestral sequences for existing alignments](docs/phylogenetic_multiple_alignment.md)
*   [Phylogenetic placement of sequences into an existing alignment](docs/pagan_phylogenetic_placement.md)
*   [Guided placement with assignment of queries to specific nodes](docs/pagan_guided_placement.md)
*   [Phylogenetic placement of amplicon (e.g. 18S) data](docs/pagan_amplicon_analysis.md)
*   [Pileup alignment without a guide tree](docs/pagan_pileup_alignment.md)
*   [Additional program options](docs/pagan_additional_program_options.md)


* * *

### PAGAN command-line options

PAGAN is a command-line program. It can be used by (a) specifying a list of options (command-line arguments) when executing the program, or (b) creating a configuration file with the options and specifying that when executing the program. The configuration file does not need to be created from scratch as PAGAN can output the options specified for an analysis in a file. This file can then be edited if necessary and specified as the configuration file for another analysis. Alternatively, the file can be considered as a record of a particular analysis with a full description of options and parameters used.

A list of the most important program options is outputted if no arguments are provided:

```
./pagan
```

and a more complete list is given with the option ```--help```:

```
./pagan --help
```

In general, the option names start with ```--``` and the option name and value (if any) are separated by a space. The configuration file makes an exception and can be specified without the option name:

```
./pagan option_file
```

Also this one can be given in the standard format and the following command is equivalent:

```
./pagan --config-file option_file
```

Configuration files contain option names and values separated by ```=``` sign, one option per row. Rows starting with a hash sign # are comments and ignored. Thus, if the content of file config.cfg is:

```
# this is an uninformative comment
ref-seqfile = reference_alignment.fas
ref-treefile = reference_tree.nhx
queryfile = illumina_reads.fastq
outfile = read_alignment
xml = 1
```

the command:

```
./pagan config.cfg
```

(or .```/pagan --config-file config.cfg```)

is equivalent to:

```
./pagan --ref-seqfile reference_alignment.fas --ref-treefile reference_tree.nhx --queryfile illumina_reads.fastq --outfile read_alignment --xml
```

By adding the option ```--config-log-file config.cfg``` in the command above, PAGAN creates a config file that is equivalent with the one above (with some more comments). Config files can, of course, be written or extended manually using the same format. One should note, however, that also boolean options need a value assigned, such as ```xml = 1``` in the example above. If a boolean option is not wanted, it should not specified in the config file (or it should be commented out with a hash sign) as setting an option e.g. ```0``` or ```false``` does not disable it.

Options in a config file are overridden by re-defining them on command line. Thus, the command:

```
./pagan config.cfg --outfile another_name
```

is the same as the one above except that the results will be placed to a file with another name.


[back to top](#pagan)
