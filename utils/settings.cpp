
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
    version = 0.007;

    boost::program_options::options_description minimal("Minimal options");
    minimal.add_options()
        ("seqfile", po::value<string>(), "sequence infile")
        ("treefile", po::value<string>(), "tree file")
        ("outfile", po::value<string>(), "sequence outfile")
    ;

    boost::program_options::options_description generic("Generic options");
    generic.add_options()
        ("help", "display help message")
        ("noise", po::value<int>(), "output noise level")
        ("output-ancestors", "include ancestors in outfile")
        ("full-probability", "compute full probability")
        ("output-graph","output ancestral graph")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
    ;

    boost::program_options::options_description reads_alignment("Reads alignment options");
    reads_alignment.add_options()
        ("cds-seqfile", po::value<string>(), "cds alignment infile")
        ("cds-treefile", po::value<string>(), "cds tree file")
        ("readsfile", po::value<string>(), "reads infile")
        ("min-reads-overlap", po::value<float>()->default_value(0.10), "read sequences' minimum overlap with reference alignment sites")
        ("reads-distance", po::value<float>()->default_value(0.01), "read sequences' evolutionary distance from root")
        ("454", "correct homopolymer error")
        ("pair-end","connect paired reads")
        ("qscore-minimum", po::value<int>()->default_value(10), "read sequences' minimum Q-score to be included")
        ("trim-read-ends", "trim read ends with low Q-scores")
        ("trim-mean-qscore", po::value<int>()->default_value(15), "sliding window average Q-score to be clipped")
        ("trim-window-width", po::value<int>()->default_value(5), "sliding window width for trimming")
        ("minimum-trimmed-length", po::value<int>()->default_value(20), "minimum trimmed read length")
        ("rank-reads-for-nodes","find best ranking for reads within nodes")
        ("discard-overlapping-identical-reads", "discard (in target) fully overlapping identical reads")
        ("discard-overlapping-reads", "discard (in target) fully overlapping reads")
        ("discard-pairwise-overlapping-reads", "discard (pairwise) fully overlapping reads")
        ("align-reads-at-root", "ignore tags and align reads at root")
        ("align-bad-reads-at-root", "align non-matching reads at root")
        ("output-nhx-tree", "output tree with NHX TID tags")
        ("allow-skip-low-qscore", "allow skipping low scoring bases")
        ("pair-read-gap-extension", po::value<float>(), "pair read middle gap extension probability")
        ("no-fastq", "do not use Q-scores")
    ;

    boost::program_options::options_description graph("Graph options");
    graph.add_options()
        ("weight-sampled-edges", "use posterior probabilities to weight sampled edges")
        ("no-weight-transform", "no weight transform for sampled edges")
        ("cuberoot-weight-transform", "cuberoot weight transform for sampled edges")
        ("no-terminal-edges", "assume terminal missing data")
    ;

    boost::program_options::options_description model("DNA model options");
    model.add_options()
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("gap-extension", po::value<float>(), "gap extension probability")
        ("end-gap-extension", po::value<float>(), "terminal gap extension probability")
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("no-log-odds", "do not use log-odds substitutions scores")
    ;

    boost::program_options::options_description tree_edit("Tree manipulation options");
    tree_edit.add_options()
        ("scale-branches", po::value<float>(), "scale tree branches")
        ("truncate-branches", po::value<float>()->default_value(0.10), "truncate tree branches")
        ("real-branches", "use real tree branch lengths")
        ("fixed-branches", po::value<float>(), "use fixed-length tree branches")
    ;

    boost::program_options::options_description alignment("Alignment model options");
    alignment.add_options()
        ("any-skips-confirm-insertion", po::value<int>(), "number of skips to confirm site as insertion")
        ("match-skips-confirm-insertion", po::value<int>(), "number of skips from flanking matches to confirm site as insertion")
        ("branch-length-confirm-insertion", po::value<float>(), "total branch length skipped to confirm site as insertion")
        ("branch-skip-weight-per-distance", po::value<float>(), "weighted (by branch length unit) probability for site(s) being skipped over and later matched (>default<)")
        ("branch-skip-penalty-per-branch", po::value<float>(), "fixed probability for site(s) being skipped over and later matched")
//        ("terminal-gap-cost-divider", po::value<float>()->default_value(2), "divider for terminal gap cost")
    ;
//        ("unused-edges-not-transferred", "mimic PRANK-F behaviour")

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
    ;

    desc.add(minimal).add(generic).add(reads_alignment).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(debug);
    min_desc.add(minimal);

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(is("noise"))
        noise = get("noise").as<int>();


    if (vm.count("help")) {
        this->help();
        return 1;
    }


    return 0;
}

void Settings::help()
{
        cout<<"\nPagan v. "<<version<<": input arguments:\n";
        cout << desc << "\n";
        exit(0);
}

void Settings::info()
{
        cout<<"\nPagan v. "<<version<<": minimal input arguments:\n";
        cout << min_desc << "\n";
        cout<<"use option --help for more information\n\n";
        exit(0);
}

int     Settings::noise             = 0;
float   Settings::resize_factor     = 1.5;
