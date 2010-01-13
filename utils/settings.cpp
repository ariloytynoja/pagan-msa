
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
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
        ("help", "display help message")
        ("seqfile", po::value<string>(), "sequence infile")
        ("treefile", po::value<string>(), "tree file")
        ("outfile", po::value<string>(), "sequence outfile")
        ("output-ancestors", "include ancestors in outfile")
        ("full-probability", "compute full probability")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("noise", po::value<int>(), "output noise level")
    ;

    boost::program_options::options_description graph("Graph options");
    graph.add_options()
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
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
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("use-log-odds", "use log-odds substitutions scores")
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
        ("terminal-gap-cost-divider", po::value<float>()->default_value(10), "divider for terminal gap cost")
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

    desc.add(generic).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(debug);

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(is("noise"))
        noise = get("noise").as<int>();


    if (vm.count("help")) {
        this->info();
        return 1;
    }


    return 0;
}

void Settings::info()
{
        float version = 0.001;

        cout<<"\nPapaya v. "<<version<<": input arguments:\n";
        cout << desc << "\n";
}

int     Settings::noise             = 0;
float   Settings::resize_factor     = 1.5;
