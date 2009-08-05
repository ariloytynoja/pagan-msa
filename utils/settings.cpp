
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
        ("noise", po::value<int>(), "output noise level")
    ;

    boost::program_options::options_description model("DNA model options");
    model.add_options()
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("gap-extension", po::value<float>(), "gap extension probability")
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
    ;

    boost::program_options::options_description alignment("Alignment model options");
    alignment.add_options()
        ("any-skips-confirm-insertion", po::value<int>(), "number of skips to confirm site as insertion")
        ("match-skips-confirm-insertion", po::value<int>(), "number of skips from flanking matches to confirm site as insertion")
        ("branch-length-confirm-insertion", po::value<float>(), "total branch length skipped to confirm site as insertion")
        ("branch-skip-weight-per-distance", po::value<float>(), "weighted (by branch length unit) probability for site(s) being skipped over and later matched (>default<)")
        ("branch-skip-penalty-per-branch", po::value<float>(), "fixed probability for site(s) being skipped over and later matched")
    ;
//        ("unused-edges-not-transferred", "mimic PRANK-F behaviour")

    boost::program_options::options_description graphs("Graph output options");
    graphs.add_options()
        ("mpost-graphfile", po::value<string>(), "sequence graphfile for metapost")
        ("output-alignment-graphs", "include aligned graphs")
        ("output-leaf-graphs", "include terminal sequences")
    ;

    boost::program_options::options_description debug("Debugging options");
    debug.add_options()
        ("check-valid-graphs", "check that fwd and bwd edges are identical")
    ;

    desc.add(generic).add(model).add(alignment).add(graphs).add(debug);

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
