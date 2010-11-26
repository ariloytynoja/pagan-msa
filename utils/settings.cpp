#include <cstdio>
#include <iostream>
#include "settings.h"

namespace po = boost::program_options;

using namespace std;
using namespace ppa;

Settings::Settings()
{

}

int Settings::read_command_line_arguments(int argc, char *argv[])
{
    version = 0.12;
    date = "24 Nov, 2010";

    boost::program_options::options_description minimal("Minimal progressive alignment options");
    minimal.add_options()
        ("seqfile", po::value<string>(), "sequence infile (FASTA)")
        ("treefile", po::value<string>(), "tree file")
        ("outfile", po::value<string>(), "sequence outfile")
    ;

    boost::program_options::options_description generic("Generic options");
    generic.add_options()
        ("help", "display help message")
        ("silent","minimal output")
        ("noise", po::value<int>(), "output noise level")
        ("output-ancestors", "include ancestors in outfile")
        ("full-probability", "compute full probability")
        ("output-graph","output ancestral graph")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("no-terminal-edges", "expect terminal missing data")
    ;

    boost::program_options::options_description reads_alignment("Basic reads alignment options");
    reads_alignment.add_options()
        ("ref-seqfile", po::value<string>(), "reference alignment file (FASTA)")
        ("ref-treefile", po::value<string>(), "reference tree file (NH/NHX)")
        ("readsfile", po::value<string>(), "reads file (FASTA/FASTQ)")
        ("outfile", po::value<string>(), "sequence outfile")
        ("pair-end","connect paired reads")
        ("454", "correct homopolymer error")
        ("test-every-node","test every node for each read")
    ;

    boost::program_options::options_description reads_alignment3("Overlapping pair reads options");
    reads_alignment3.add_options()
        ("overlap-pair-end","merge overlapping paired reads")
        ("overlap-minimum", po::value<int>()->default_value(15), "minimum overlapping sites")
        ("overlap-identity", po::value<float>()->default_value(0.75), "minimum identity at overlap")
        ("overlap-identical-minimum", po::value<int>()->default_value(10), "minimum identical overlapping sites")
//        ("overlap-no-truncate","do not truncate suspicious overlap")
        ("overlap-merge-file", po::value<string>(), "output file for merged reads")
        ("overlap-merge-only","only merge overlapping paired reads")
    ;

    boost::program_options::options_description reads_alignment2("Additional reads alignment options");
    reads_alignment2.add_options()
        ("no-fastq", "do not use Q-scores")
        ("qscore-minimum", po::value<int>()->default_value(10), "threshold to mask low Q-score sites")
        ("trim-read-ends", "trim read ends with low Q-scores")
        ("trim-mean-qscore", po::value<int>()->default_value(15), "sliding window trimming threshold")
        ("trim-window-width", po::value<int>()->default_value(5), "sliding window width")
        ("minimum-trimmed-length", po::value<int>()->default_value(20), "minimum trimmed read length")
        ("rank-reads-for-nodes","rank reads within nodes for alignment")
        ("discard-overlapping-identical-reads", "discard embedded identical reads")
        ("discard-overlapping-reads", "discard embedded reads")
        ("discard-pairwise-overlapping-reads", "discard embedded reads (pairwise alignment)")
        ("align-reads-at-root", "ignore tags and align reads at root")
        ("align-bad-reads-at-root", "align non-matching reads at root")
        ("placement-only", "compute read placement only")
        ("placement-file", po::value<string>(), "read placement file")
        ("output-nhx-tree", "output tree with NHX TID tags")
        ("allow-skip-low-qscore", "allow skipping low scoring bases")
        ("pair-read-gap-extension", po::value<float>(), "read spacer extension probability")
        ("min-reads-overlap", po::value<float>()->default_value(0.1), "overlap threshold for read and reference")
        ("min-reads-identity", po::value<float>()->default_value(0.3), "identity threshold for aligned sites")
        ("reads-distance", po::value<float>()->default_value(0.01), "evolutionary distance from pseudo-root")
    ;

    boost::program_options::options_description graph("Graph options");
    graph.add_options()
        ("weight-sampled-edges", "use posterior scores to weight sampled edges")
        ("no-weight-transform", "no weight transform for sampled edges")
        ("cuberoot-weight-transform", "cuberoot weight transform for sampled edges")
    ;

    boost::program_options::options_description model("DNA model options");
    model.add_options()
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("gap-extension", po::value<float>(), "gap extension probability")
        ("end-gap-extension", po::value<float>(), "terminal gap extension probability")
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
    ;

    boost::program_options::options_description tree_edit("Tree manipulation options");
    tree_edit.add_options()
        ("scale-branches", po::value<float>(), "scale tree branches")
        ("truncate-branches", po::value<float>()->default_value(0.1), "truncate tree branches")
        ("real-branches", "use real tree branch lengths")
        ("fixed-branches", po::value<float>(), "use fixed-length tree branches")
    ;

    boost::program_options::options_description alignment("Alignment model options");
    alignment.add_options()
        ("any-skips-confirm-insertion", po::value<int>(), "#skips to confirm as insertion")
        ("match-skips-confirm-insertion", po::value<int>(), "#skips from match sites to confirm as insertion")
        ("branch-length-confirm-insertion", po::value<float>(), "total branch length skipped to confirm as insertion")
        ("branch-skip-weight-per-distance", po::value<float>(), "weighted (by branch length unit) probability for site(s) being skipped over and later matched (>default<)")
        ("branch-skip-penalty-per-branch", po::value<float>(), "fixed probability for site(s) being skipped over and later matched")
    ;

    boost::program_options::options_description output("Graph output options");
    output.add_options()
        ("mpost-graph-file", po::value<string>(), "sequence graph file for metapost")
        ("output-alignment-graphs", "include aligned graphs")
        ("output-leaf-graphs", "include terminal sequences")
        ("mpost-posterior-plot-file", po::value<string>(), "posterior plot file for metapost")
        ("plot-slope-up", "plot viterbi path climbing up")
    ;

    boost::program_options::options_description debug("Debugging options");
    debug.add_options()
        ("check-valid-graphs", "check that fwd and bwd edges are identical")
        ("help-all", "display help-all message")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("no-log-odds", "do not use log-odds substitutions scores")
    ;

    boost::program_options::options_description broken("Broken options");
    broken.add_options()
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
        ("cds-seqfile", po::value<string>(), "reference alignment file (FASTA)")
        ("cds-treefile", po::value<string>(), "reference tree file (NH/NHX)")
//        ("terminal-gap-cost-divider", po::value<float>()->default_value(2), "divider for terminal gap cost")
//        ("unused-edges-not-transferred", "mimic PRANK-F behaviour")
    ;

    full_desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(debug).add(broken);
    desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(tree_edit).add(alignment);
    max_desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(graph).add(tree_edit).add(alignment).add(output);
    min_desc.add(minimal).add(reads_alignment);

    po::store(po::parse_command_line(argc, argv, full_desc), vm);

    if(Settings::is("cds-seqfile") || Settings::is("cds-treefile"))
    {
        cout<<"\nThe program options '--cds-seqfile' and '--cds-treefile' have been renamed as '--ref-seqfile' and '--ref-treefile'.\n"
                "Please edit your command argument line accordingly. Exiting.\n\n";
        exit(0);

    }
    po::notify(vm);

    if(is("noise"))
        noise = get("noise").as<int>();


    if (vm.count("help")) {
        this->help();
        return 1;
    }

    if (vm.count("help-all")) {
        this->help_all();
        return 1;
    }


    return 0;
}

void Settings::print_msg()
{
    cout<<"\nPAGAN v. "<<version<<" ("<<date<<"). (C) 2010 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    cout<<" This is a development version and may contain bugs. Contact the author\n before using the program for any serious analysis.\n";
}

void Settings::help()
{
    this->print_msg();
    cout<< desc << "\n";
    exit(0);
}

void Settings::help_all()
{
    this->print_msg();
    cout<< max_desc << "\n";
    exit(0);
}

void Settings::info()
{
    this->print_msg();
    cout << min_desc << "\n";
    cout<<"Use option --help for more information.\n\n";
    exit(0);
}

int     Settings::noise             = 0;
float   Settings::resize_factor     = 1.5;
