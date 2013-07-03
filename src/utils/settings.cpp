/***************************************************************************
 *   Copyright (C) 2010-2012 by Ari Loytynoja                              *
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
#include "utils/log_output.h"

namespace po = boost::program_options;

using namespace std;
using namespace ppa;

Settings::Settings(){}

int Settings::read_command_line_arguments(int argc, char *argv[])
{
    version = 0.50;
    date = "15 May, 2013";

    boost::program_options::options_description minimal("Minimal progressive alignment options",100);
    minimal.add_options()
        ("seqfile", po::value<string>(), "sequence infile (FASTA)")
        ("treefile", po::value<string>(), "tree file")
    ;

    boost::program_options::options_description generic("Generic options",100);
    generic.add_options()
        ("outfile", po::value<string>(), "alignment outfile")
        ("outformat", po::value<string>(), "alignment format [fasta,nexus,paml,phylipi,phylips,raxml]")
        ("raxml-tree","use RAxML for guide tree computation [default BppDist]")
        ("translate", "translate DNA input to protein")
        ("mt-translate", "translate mtDNA input to protein")
        ("output-ancestors", "include ancestors in outfile")
        ("config-file",po::value<string>(),"config file with additional arguments")
        ("config-log-file",po::value<string>(),"log file for given arguments")
        ("no-terminal-edges", "assume terminal gaps as missing data")
        ("silent","minimal output")
        ("noise", po::value<int>(), "output noise level")
        ("log-output-file",po::value<string>(),"output to file instead of stdout")
        ("xml","output also XML alignment")
        ("output-nhx-tree", "output alignment tree (with NHX TID tags)")
        ("temp-folder",po::value<string>(),"non-standard place for temp files")
        ("threads",po::value<int>(),"number of threads")
    ;
    boost::program_options::options_description help_update("Help and updates",100);
    help_update.add_options()
        ("help", "display all program options")
        ("version","show program version and check for updates")
    ;

    boost::program_options::options_description reads_alignment("Basic alignment extension options",100);
    reads_alignment.add_options()
        ("ref-seqfile", po::value<string>(), "reference alignment file (FASTA)")
        ("ref-treefile", po::value<string>(), "reference tree file (NH/NHX)")
        ("queryfile", po::value<string>(), "query file (FASTA/FASTQ)")
        ("454", "correct homopolymer error (DNA)")
        ("homopolymer", "correct homopolymer error (more agressively)")
        ("use-consensus", "use consensus for query ancestors")
        ("build-contigs", "build contigs of query clusters")
        ("test-every-node","test every node for each query")
        ("fast-placement","use Exonerate to quickly assign queries to nodes")
        ("very-fast-placement","shorthand for fast heuristic settings")
//        ("tagged-search","search over tagged nodes")
        ("upwards-search","stepwise search from root")
        ("prune-extended-alignment","remove closely related sequences")
        ("prune-keep-number",po::value<int>()->default_value(15),"keep N most distantly related sequences (default)")
        ("prune-keep-threshold",po::value<float>(),"remove sequences with distance below threshold")
    ;

    boost::program_options::options_description reads_alignment2("Additional alignment extension options",100);
    reads_alignment2.add_options()
        ("pair-end","connect paired reads (FASTQ)")
        ("pacbio","correct for missing data in PacBio reads (DNA)")
        ("show-contig-ancestor", "fill contig gaps with ancestral sequence")
        ("consensus-minimum", po::value<int>()->default_value(5), "threshold for inclusion in contig")
        ("consensus-minimum-proportion", po::value<float>()->default_value(0.5,"0.5"), "threshold for inclusion in contig")
        ("test-every-internal-node","test every internal node for each query")
        ("test-every-terminal-node","test every terminal node for each query")
        ("one-placement-only", "place only once despite equally good hits")
        ("placement-only", "compute query placement only")
        ("placement-file", po::value<string>(), "query placement file")
        ("query-distance", po::value<float>()->default_value(0.1,"0.1"), "evolutionary distance from pseudo-root")
        ("min-query-overlap", po::value<float>()->default_value(0.5,"0.5"), "overlap threshold for query and reference")
        ("overlap-with-reference","require overlap with reference")
        ("min-query-identity", po::value<float>()->default_value(0.5,"0.5"), "identity threshold for aligned sites")
        ("pair-read-gap-extension", po::value<float>(), "paired read spacer extension probability (DNA)")
    ;

    boost::program_options::options_description reads_alignment3("Overlapping paired reads options (FASTQ)",100);
    reads_alignment3.add_options()
        ("overlap-pair-end","merge overlapping paired reads")
        ("overlap-minimum", po::value<int>()->default_value(15), "minimum overlapping sites")
        ("overlap-identity", po::value<float>()->default_value(0.75,"0.75"), "minimum identity at overlap")
        ("overlap-identical-minimum", po::value<int>()->default_value(10), "minimum identical overlapping sites")
        ("overlap-merge-file", po::value<string>(), "output file for merged reads")
        ("overlap-merge-only","only merge overlapping paired reads")
        ("trim-before-merge", "trim read ends with low Q-scores before merging")
    ;

    boost::program_options::options_description reads_alignment4("Trimming and quality options (FASTQ)",100);
    reads_alignment4.add_options()
        ("no-fastq", "do not use Q-scores")
        ("trim-read-ends", "trim read ends with low Q-scores")
        ("trim-mean-qscore", po::value<int>()->default_value(15), "sliding window trimming threshold")
        ("trim-window-width", po::value<int>()->default_value(5), "sliding window width")
        ("qscore-minimum", po::value<int>()->default_value(10), "threshold to mask low Q-score sites")
        ("minimum-trimmed-length", po::value<int>()->default_value(20), "minimum trimmed read length")
        ("perfect-reference", "assume perfect reference alignment")
    ;

    boost::program_options::options_description anchoring("Anchoring options",100);
    anchoring.add_options()
        ("use-anchors","use Exonerate to anchor alignment")
        ("exonerate-hit-length", po::value<int>()->default_value(30), "Exonerate hit length for anchor (default)")
        ("exonerate-hit-score", po::value<int>(), "Exonerate hit score for anchor")
        ("anchors-offset", po::value<int>()->default_value(15), "anchors offset for alignment")
        ("use-prefix-anchors","use prefix approach to anchor alignment")
        ("prefix-hit-length", po::value<int>()->default_value(30), "prefix hit length for anchor")
        ("keep-temp-files","keep temporary files (mainly for debugging)")
    ;

    boost::program_options::options_description exonerate("Exonerate options",100);
    exonerate.add_options()
        ("exhaustive-placement","if Exonrate fails, use PAGAN to place the query")
        ("use-exonerate-local","use Exonerate local to map queries to nodes")
        ("exonerate-local-keep-best",po::value<int>()->default_value(5),"keep best # of local matches")
        ("exonerate-local-keep-above",po::value<float>(),"keep local matches above #% of the best score")
        ("use-exonerate-gapped","use Exonerate gapped to map queries to nodes")
        ("exonerate-gapped-keep-best",po::value<int>()->default_value(1),"keep best # of gapped matches")
        ("exonerate-gapped-keep-above",po::value<float>(),"keep gapped matches above #% of the best score")
        ("keep-despite-exonerate-fails", "keep queries that Exonerate fails to align")
    ;

    boost::program_options::options_description pileup("Pileup alignment options",100);
    pileup.add_options()
        ("pileup-alignment","make pileup alignment")
        ("compare-reverse","test also reverse-complement and keep better (DNA)")
        ("find-best-orf", "translate and use best ORF (DNA)")
        ("min-orf-coverage", po::value<float>()->default_value(0.95,"0.95"),"minimum ORF coverage to be considered (DNA)")
        ("query-cluster-attempts", po::value<int>()->default_value(1),"attempts to find overlap")
        ("find-cluster-reference","find optimal cluster reference")
        ("inlude-parent-in-contig", "include also ancestral parent in contigs")
        ("use-duplicate-weights", "use NumDuplicates=# to weight consensus counts")
        ("no-read-ordering","do not order reads by duplicate number")
        ("output-consensus", "output contig consensus alone")
    ;

    boost::program_options::options_description obscure("Additional obscure options",100);
    obscure.add_options()
        ("rank-reads-for-nodes","rank reads within nodes for alignment")
        ("discard-overlapping-identical-reads", "discard embedded identical reads")
        ("discard-overlapping-reads", "discard embedded reads")
        ("discard-pairwise-overlapping-reads", "discard embedded reads (pairwise alignment)")
        ("align-reads-at-root", "ignore tags and align reads at root")
        ("align-bad-reads-at-root", "align non-matching reads at root")
        ("use-identity-score", "choose target based on identity score")
        ("use-target-normalised-score", "choose target based on target-normalised substitution score")
        ("anchoring-coverage-threshold",po::value<float>()->default_value(1.0,"1.00"),"anchoring coverage threshold for skipping")
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
        ("truncate-branches", po::value<float>()->default_value(0.2,"0.2"), "truncate tree branches")
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
        ("keep-all-edges","nothing of those -- keep everything forever")
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
        ("full-probability", "compute full probability")
        ("output-graph","output ancestral graph")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("check-valid-graphs", "check that fwd and bwd edges are identical")
        ("full-help", "display full-help message")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("no-log-odds", "do not use log-odds substitutions scores")
        ("time", "track time (debugging)")
        ("recompute-reference-alignment-model", "recompute reference alignment model")
        ("no-score-scaling","no subsistitution score scaling")
        ("plot-anchors-for-R","plot for R")
    ;

    boost::program_options::options_description broken("Broken options",100);
    broken.add_options()
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
        ("readsfile", po::value<string>(), "reads file (FASTA/FASTQ)")
        ("reads-pileup","make pileup alignment")
        ("cluster-pileup","pileup clustered reads (DNA)")
    ;

    po::positional_options_description pd;
    pd.add("config-file", 1);

    full_desc.add(minimal).add(generic).add(anchoring).add(reads_alignment).add(reads_alignment2).add(reads_alignment3).add(reads_alignment4).add(exonerate).add(pileup).add(obscure).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(debug).add(broken).add(help_update);
    desc.add(minimal).add(generic).add(anchoring).add(reads_alignment).add(reads_alignment2).add(reads_alignment3).add(reads_alignment4).add(exonerate).add(pileup).add(model).add(tree_edit).add(alignment).add(help_update);
    max_desc.add(minimal).add(generic).add(anchoring).add(reads_alignment).add(reads_alignment2).add(reads_alignment3).add(reads_alignment4).add(exonerate).add(pileup).add(obscure).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(help_update);
    min_desc.add(minimal).add(reads_alignment).add(help_update);


//    po::store(po::parse_command_line(argc, argv, full_desc), vm);
    po::store(po::command_line_parser(argc, argv).options(full_desc).positional(pd).run(), vm);

    if(Settings::is("config-file"))
    {
        stringstream ss;
        ss<<"\nReading command line options from file '" << Settings::get("config-file").as<string>()<<"'.\n";
        Log_output::write_out(ss.str(),0);

        ifstream cfg(Settings::get("config-file").as<string>().c_str());

        if (!cfg)
        {
            Settings::info_noexit();

            stringstream ss;
            ss<<"Config file '" << Settings::get("config-file").as<string>()<<"' not found. Exiting.\n\n";
            Log_output::write_out(ss.str(),0);

            exit(1);
        }


        po::store(po::parse_config_file(cfg, full_desc), vm);
    }

    if(is("readsfile"))
    {
        Settings::info_noexit();
        stringstream ss;
        ss<<"\nThe program option '--readsfile' has been renamed as '--queryfile'. Please edit your command argument. Exiting.\n\n";
        Log_output::write_out(ss.str(),0);
        exit(1);
    }
    if(is("reads-pileup"))
    {
        Settings::info_noexit();
        stringstream ss;
        ss<<"\nThe program option '--reads-pileup' has been renamed as '--pileup-alignment'. Please edit your command argument. Exiting.\n\n";
        Log_output::write_out(ss.str(),0);
        exit(1);
    }

    po::notify(vm);

    if(is("noise"))
        noise = get("noise").as<int>();

    if(is("silent"))
        noise = -1;

    if(is("exonerate-local-keep-best") && get("exonerate-local-keep-best").as<int>() > 0 )
        exonerate_local_keep_best = get("exonerate-local-keep-best").as<int>();

    if(is("exonerate-gapped-keep-best") && get("exonerate-gapped-keep-best").as<int>() > 0 )
        exonerate_gapped_keep_best = get("exonerate-gapped-keep-best").as<int>();

    if(is("very-fast-placement"))
    {
        exonerate_local_keep_best = 1;
        exonerate_gapped_keep_best = 1;

        if( is("use-exonerate-local") || /*is("exonerate-local-keep-best") ||*/ is("exonerate-local-keep-above") ||
            is("use-exonerate-gapped") || /*is("exonerate-gapped-keep-best") ||*/ is("exonerate-gapped-keep-above") ||
            is("exhaustive-placement") || is("keep-despite-exonerate-fails") )
        {
            Log_output::write_out("Please disable Exonerate-related options if using '--very-fast-placement'. Exiting.\n",0);
            exit(1);
        }
    }

    if( (is("very-fast-placement") || is("fast-placement")) && !(is("test-every-node") || is("test-every-internal-node") || is("test-every-terminal-node") ) )
    {
        Log_output::write_out("\nWarning: When using option '--(very-)fast-placement', either option '--test-every-[internal-|terminal-]node'\nor a reference tree with TID tags should be used. "
                              "If the latter is true, this warning can be ignored.\n",0);
    }

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
        stringstream ss;
        ss << endl<< "Writing command line options to file '" << Settings::get("config-log-file").as<string>()<<"'."<<endl;
        Log_output::write_out(ss.str(),1);

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
    stringstream ss;
    ss<<"\nPAGAN v."<<version<<" ("<<date<<"). (C) 2010-2012 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    ss<<" This program is provided \"as-is\", with NO WARRANTY whatsoever; this is a development version\n and may contain bugs.\n";
    Log_output::write_out(ss.str(),0);
}

string Settings::print_log_msg()
{
    stringstream tmp;
    tmp<<"\n# PAGAN v."<<version<<" ("<<date<<"). (C) 2010-2012 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    tmp<<"# This program is provided \"as-is\", with NO WARRANTY whatsoever; this is a development version\n# and may contain bugs.\n";
    return tmp.str();
}

void Settings::help()
{
    this->print_msg();
    stringstream ss;
    ss<< desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::help_all()
{
    this->print_msg();
    stringstream ss;
    ss<< max_desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::info()
{
    this->print_msg();
    stringstream ss;
    ss<< min_desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::info_noexit()
{
    this->print_msg();
    stringstream ss;
    ss<< min_desc << "\n\n";
    Log_output::write_out(ss.str(),0);
}

void Settings::check_version()
{
    Check_version cv(version);
    exit(0);
}

int   Settings::noise             = 0;
float Settings::resize_factor     = 1.5;

int   Settings::exonerate_local_keep_best = 1;
int   Settings::exonerate_gapped_keep_best = 1;
