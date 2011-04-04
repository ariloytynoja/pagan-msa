/***************************************************************************
 *   Copyright (C) 2010 by Ari Loytynoja                                   *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include "utils/settings.h"
#include "utils/check_version.h"

namespace po = boost::program_options;

using namespace std;
using namespace ppa;

Settings::Settings()
{

}

int Settings::read_command_line_arguments(int argc, char *argv[])
{
    version = 0.24;
    date = "1 Apr, 2011";

    boost::program_options::options_description minimal("Minimal progressive alignment options",100);
    minimal.add_options()
        ("seqfile", po::value<string>(), "sequence infile (FASTA)")
        ("treefile", po::value<string>(), "tree file")
    ;

    boost::program_options::options_description generic("Generic options",100);
    generic.add_options()
        ("help", "display help message")
        ("outfile", po::value<string>(), "sequence outfile")
        ("output-ancestors", "include ancestors in outfile")
        ("full-probability", "compute full probability")
        ("output-graph","output ancestral graph")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("no-terminal-edges", "expect terminal missing data")
        ("noise", po::value<int>(), "output noise level")
        ("silent","minimal output")
        ("no-xml","write FASTA alignment only")
        ("config-file",po::value<string>(),"config file with additional arguments")
        ("config-log-file",po::value<string>(),"log file for given arguments")
        ("version","show program version and check for updates")
    ;

    boost::program_options::options_description reads_alignment("Basic reads alignment options",100);
    reads_alignment.add_options()
        ("ref-seqfile", po::value<string>(), "reference alignment file (FASTA)")
        ("ref-treefile", po::value<string>(), "reference tree file (NH/NHX)")
        ("readsfile", po::value<string>(), "reads file (FASTA/FASTQ)")
        ("pair-end","connect paired reads")
        ("454", "correct homopolymer error")
        ("test-every-node","test every node for each read")
    ;

    boost::program_options::options_description reads_alignment3("Overlapping pair reads options",100);
    reads_alignment3.add_options()
        ("overlap-pair-end","merge overlapping paired reads")
        ("overlap-minimum", po::value<int>()->default_value(15), "minimum overlapping sites")
        ("overlap-identity", po::value<float>()->default_value(0.75), "minimum identity at overlap")
        ("overlap-identical-minimum", po::value<int>()->default_value(10), "minimum identical overlapping sites")
        ("overlap-merge-file", po::value<string>(), "output file for merged reads")
        ("overlap-merge-only","only merge overlapping paired reads")
        ("trim-before-merge", "trim read ends with low Q-scores before merging")
    ;

    boost::program_options::options_description reads_alignment2("Additional reads alignment options",100);
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
        ("reads-distance", po::value<float>()->default_value(0.1), "evolutionary distance from pseudo-root")
        ("perfect-reference", "assume perfect reference alignment")
    ;

    boost::program_options::options_description graph("Graph options",100);
    graph.add_options()
        ("weight-sampled-edges", "use posterior scores to weight sampled edges")
        ("no-weight-transform", "no weight transform for sampled edges")
        ("cuberoot-weight-transform", "cuberoot weight transform for sampled edges")
    ;

    boost::program_options::options_description model("DNA/Protein model options",100);
    model.add_options()
        ("indel-rate", po::value<float>(), "insertion-deletion rate")
        ("gap-extension", po::value<float>(), "gap extension probability")
        ("end-gap-extension", po::value<float>(), "terminal gap extension probability")
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
        ("codons", "translate and align codons")
        ("use-aa-groups", "reconstruct amino-acid parsimony with 51 groups")
    ;

    boost::program_options::options_description tree_edit("Tree manipulation options",100);
    tree_edit.add_options()
        ("scale-branches", po::value<float>(), "scale tree branches")
        ("truncate-branches", po::value<float>()->default_value(0.1), "truncate tree branches")
        ("real-branches", "use real tree branch lengths")
        ("fixed-branches", po::value<float>(), "fixed length for tree branches")
        ("min-branch-length", po::value<float>(), "minimum length for tree branches")
    ;

    boost::program_options::options_description alignment("Alignment model options",100);
    alignment.add_options()
        ("any-skips-confirm-insertion", po::value<int>(), "#skips to confirm as insertion")
        ("match-skips-confirm-insertion", po::value<int>(), "#skips from match sites to confirm as insertion")
        ("branch-length-confirm-insertion", po::value<float>(), "total branch length skipped to confirm as insertion")
        ("branch-skip-weight-per-distance", po::value<float>(), "weighted (by branch length unit) probability for site(s) being skipped over and later matched (>default<)")
        ("branch-skip-penalty-per-branch", po::value<float>(), "fixed probability for site(s) being skipped over and later matched")
    ;

    boost::program_options::options_description output("Graph output options",100);
    output.add_options()
        ("mpost-graph-file", po::value<string>(), "sequence graph file for metapost")
        ("output-alignment-graphs", "include aligned graphs")
        ("output-leaf-graphs", "include terminal sequences")
        ("mpost-posterior-plot-file", po::value<string>(), "posterior plot file for metapost")
        ("plot-slope-up", "plot viterbi path climbing up")
    ;

    boost::program_options::options_description debug("Debugging and testing options",100);
    debug.add_options()
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("check-valid-graphs", "check that fwd and bwd edges are identical")
        ("full-help", "display full-help message")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("no-log-odds", "do not use log-odds substitutions scores")
        ("time", "track time (debugging)")
    ;

    boost::program_options::options_description broken("Broken options",100);
    broken.add_options()
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
        ("cds-seqfile", po::value<string>(), "reference alignment file (FASTA)")
        ("cds-treefile", po::value<string>(), "reference tree file (NH/NHX)")
    ;

    po::positional_options_description pd;
    pd.add("config-file", 1);

    full_desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(debug).add(broken);
    desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(tree_edit).add(alignment);
    max_desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment3).add(reads_alignment2).add(model).add(graph).add(tree_edit).add(alignment).add(output);
    min_desc.add(minimal).add(reads_alignment);


//    po::store(po::parse_command_line(argc, argv, full_desc), vm);
    po::store(po::command_line_parser(argc, argv).options(full_desc).positional(pd).run(), vm);

    if(Settings::is("config-file"))
    {
        cout << endl<< "Reading command line options from file '" << Settings::get("config-file").as<string>()<<"'."<<endl;

        ifstream cfg(Settings::get("config-file").as<string>().c_str());
        po::store(po::parse_config_file(cfg, full_desc), vm);
    }

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

    if (vm.count("full-help")) {
        this->help_all();
        return 1;
    }

    if (vm.count("version")) {
        this->check_version();
        return 1;
    }



    if(Settings::is("config-log-file"))
    {
        cout << endl<< "Writing command line options to file '" << Settings::get("config-log-file").as<string>()<<"'."<<endl;

        ofstream log_out(Settings::get("config-log-file").as<string>().c_str());
        time_t s_time;
        time( &s_time );
        log_out <<this->print_log_msg()<< "#\n# Analysis started: " << asctime( localtime( &s_time ) );
        log_out<<"# Command line arguments:"<<endl<<endl;

        po::parsed_options opts = parse_command_line(argc, argv, full_desc);

        typedef vector< po::basic_option<char> > vec_opt;

        for(vec_opt::iterator iter = opts.options.begin(); iter != opts.options.end(); ++iter)
        {
            po::basic_option<char>& option = *iter;

            if(option.string_key == "config-log-file" || option.string_key == "config-file")
                continue;

            stringstream ss_opt;
            typedef vector< basic_string<char> > vec_string;

            for(vec_string::iterator s_iter = option.value.begin(); s_iter != option.value.end(); ++s_iter)
            {
                    ss_opt << *s_iter;
            }
            if(ss_opt.str().length()>0)
                log_out<<option.string_key << " = " << ss_opt.str()<<endl;
            else
                log_out<<option.string_key << " = 1"<<endl;
        }

        if(Settings::is("config-file"))
        {
            log_out<<"\n# Additional arguments from file '"<<Settings::get("config-file").as<string>()<<"':"<<endl<<endl;

            ifstream cfg(Settings::get("config-file").as<string>().c_str());
            opts = po::parse_config_file(cfg, full_desc);

            for(vec_opt::iterator iter = opts.options.begin(); iter != opts.options.end(); ++iter)
            {
                po::basic_option<char>& option = *iter;

                if(option.string_key == "config-log-file" || option.string_key == "config-file")
                    continue;

                stringstream ss_opt;
                typedef vector< basic_string<char> > vec_string;

                for(vec_string::iterator s_iter = option.value.begin(); s_iter != option.value.end(); ++s_iter)
                {
                        ss_opt << *s_iter;
                }
                if(ss_opt.str().length()>0)
                    log_out<<option.string_key << " = " << ss_opt.str()<<endl;
                else
                    log_out<<option.string_key << " = 1"<<endl;
            }
        }
        log_out<<endl;
    }

    return 0;
}

void Settings::print_msg()
{
    cout<<"\nPAGAN v."<<version<<" ("<<date<<"). (C) 2010-2011 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    cout<<" This is a development version and may contain bugs. Contact the author before using\n the program for any serious analysis.\n";
}

string Settings::print_log_msg()
{
    stringstream tmp;
    tmp<<"\n# PAGAN v."<<version<<" ("<<date<<"). (C) 2010-2011 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    tmp<<"# This is a development version and may contain bugs. Contact the author before using\n# the program for any serious analysis.\n";
    return tmp.str();
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

void Settings::info_noexit()
{
    this->print_msg();
    cout << min_desc << "\n";
    cout<<"Use option --help for more information.\n\n";
}

void Settings::check_version()
{
    Check_version cv(version);
    exit(0);
}

int     Settings::noise             = 0;
float   Settings::resize_factor     = 1.5;
