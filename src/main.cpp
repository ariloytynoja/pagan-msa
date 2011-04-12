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

#include <string>
#include <vector>
#include <ctime>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/newick_reader.h"
#include "utils/fasta_reader.h"
#include "utils/xml_writer.h"
#include "utils/model_factory.h"
#include "utils/evol_model.h"
#include "main/node.h"
#include "main/reads_aligner.h"

using namespace std;

using namespace ppa;

int main(int argc, char *argv[])
{
    clock_t t_start=clock();

    /*
    / Input: command line arguments & data
   */

    // Read the arguments
    try
    {
        int rv = Settings_handle::st.read_command_line_arguments(argc, argv);
    }
    catch ( const boost::program_options::error& e ) {
            Settings_handle::st.info_noexit();
            cout<<"Error in command line arguments: "<<e.what()<<"."<<endl<<endl;
            exit(0);
    }

    srand(time(0));

    clock_t analysis_start_time=clock();

    /***********************************************************************/
    /*  Overlapping paired-end read merge only                             */
    /***********************************************************************/

    if(Settings_handle::st.is("overlap-merge-only"))
    {

        if(!Settings_handle::st.is("readsfile"))
        {
            cout<<"No reads file given. Exiting.\n\n";
            exit(0);
        }

        Reads_aligner ra;
        ra.merge_reads_only();

        exit(0);
    }


    /***********************************************************************/
    /*  Read the sequences                                                 */
    /***********************************************************************/

    bool reference_alignment = false;
    vector<Fasta_entry> sequences;
    Fasta_reader fr;

    if(Settings_handle::st.is("seqfile"))
    {
        string seqfile =  Settings_handle::st.get("seqfile").as<string>();
        cout<<"Data file: "<<seqfile<<endl;

        fr.read(seqfile, sequences, true);
    }
    else if(Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("ref-seqfile").as<string>();
        cout<<"Reference alignment file: "<<seqfile<<endl;

        fr.read(seqfile, sequences, true);
        reference_alignment = true;
    }
    else
    {
        cout<<endl<<"Error: No sequence file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }



    /***********************************************************************/
    /*  Read the guidetree                                                 */
    /***********************************************************************/

    Node *root;
    bool tree_ok = false;

    if(Settings_handle::st.is("treefile"))
    {
        string treefile =  Settings_handle::st.get("treefile").as<string>();
        cout<<"Tree file: "<<treefile<<endl;

        Newick_reader nr;
        string tree = nr.read_tree(treefile);
        root = nr.parenthesis_to_tree(tree);

        tree_ok = true;
    }
    else if(Settings_handle::st.is("ref-treefile"))
    {
        string treefile =  Settings_handle::st.get("ref-treefile").as<string>();
        cout<<"Reference tree file: "<<treefile<<endl;

        Newick_reader nr;
        string tree = nr.read_tree(treefile);
        root = nr.parenthesis_to_tree(tree);

        tree_ok = true;
    }
    else if(reference_alignment && sequences.size()==1)
    {
        int data_type = fr.check_sequence_data_type(&sequences);

        root = new Node();
        root->set_name(sequences.at(0).name);
        root->add_name_comment( sequences.at(0).comment );
        root->set_distance_to_parent(0);

        root->add_sequence( sequences.at(0), data_type);

        tree_ok = true;
    }
    else
    {
        cout<<endl<<"Error: No tree file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }

    if(!Settings_handle::st.is("silent"))
    {
        Settings_handle::st.print_msg();
        time_t s_time;
        time( &s_time );
        cout <<endl<< "The analysis started: " << asctime( localtime( &s_time ) )<<endl;
    }


    /***********************************************************************/
    /*  Check that input is fine                                           */
    /***********************************************************************/

    // Check that the guidetree and sequences match

    vector<Node*> leaf_nodes;
    root->get_leaf_nodes(&leaf_nodes);

    bool tree_branches_ok = fr.check_sequence_names(&sequences,&leaf_nodes);
    if(!tree_branches_ok)
    {
        cout<<"Attempting to prune the tree.\n";
        root->prune_tree();
        cout<<"Pruning done. New tree:"<<endl;
        cout<<root->print_tree()<<"\n\n";

    }


    leaf_nodes.clear();
    root->get_leaf_nodes(&leaf_nodes);

    // Check that the sequences are fine; also computes the base frequencies.

    int data_type = fr.check_sequence_data_type(&sequences);

    if(Settings_handle::st.is("time"))
        cout <<"Time main::input: "<<double(clock()-t_start)/CLOCKS_PER_SEC << endl;


    Model_factory mf(data_type);

    if(!fr.check_alphabet(&sequences,data_type))
        cout<<"\nWarning: Illegal characters in input sequences removed!"<<endl;

    if(data_type==Model_factory::dna && Settings_handle::st.is("codons"))
    {
        // Create a codon alignment model using K&G.
        if(Settings::noise>2)
            cout<<"creating a codon model\n";
        mf.codon_model(&Settings_handle::st); // does it need the handle????
    }
    else if(data_type==Model_factory::dna)
    {
        // Create a DNA alignment model using empirical base frequencies.
        if(Settings::noise>2)
            cout<<"creating a DNA model\n";
        float *dna_pi = fr.base_frequencies();
        mf.dna_model(dna_pi,&Settings_handle::st);
    }
    else if(data_type==Model_factory::protein)
    {
        // Create a protein alignment model using WAG.
        if(Settings::noise>2)
            cout<<"creating a protein model\n";
        mf.protein_model(&Settings_handle::st); // does it need the handle????
    }


    if(Settings_handle::st.is("time"))
        cout <<"Time main::model: "<<double(clock()-t_start)/CLOCKS_PER_SEC << endl;


    /***********************************************************************/
    /*  Place the sequences to nodes and align them!                       */
    /***********************************************************************/

    fr.place_sequences_to_nodes(&sequences,&leaf_nodes,reference_alignment,data_type);

    int count = 1;
    root->name_internal_nodes(&count);

    root->start_alignment(&mf);

    if(Settings_handle::st.is("time"))
        cout <<"Time main::align: "<<double(clock()-t_start)/CLOCKS_PER_SEC << endl;

    // If reads sequences, add them to the alignment

    if( Settings_handle::st.is("readsfile") )
    {
        Reads_aligner ra;
        ra.align(root,&mf,count);

        root = ra.get_global_root();

        if(Settings_handle::st.is("time"))
            cout <<"Time main::reads_align: "<<double(clock()-t_start)/CLOCKS_PER_SEC << endl;
    }



    /***********************************************************************/
    /*  Collect the results and output them                                */
    /***********************************************************************/

//    root->show_seqs();

    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,Settings_handle::st.is("output-ancestors"));

    // Save results in output file
    //
    string outfile =  "outfile";

    if(Settings_handle::st.is("outfile"))
        outfile =  Settings_handle::st.get("outfile").as<string>();


    if(Settings_handle::st.is("no-xml"))
        cout<<"Alignment files: "<<outfile<<".fas"<<endl;
    else
        cout<<"Alignment files: "<<outfile<<".fas, "<<outfile<<".xml"<<endl;

    fr.set_chars_by_line(70);
    fr.write(outfile, aligned_sequences, true);

    count = 1;
    root->set_name_ids(&count);

    if(!Settings_handle::st.is("no-xml"))
    {
        Xml_writer xw;
        xw.write(outfile, root, aligned_sequences, true);
    }

    if(Settings_handle::st.is("output-ancestors"))
    {
        fr.write_anctree(outfile, root);
    }

    if(Settings_handle::st.is("output-nhx-tree"))
    {
        root->write_nhx_tree(outfile);
    }

    if(Settings::noise>1 )
    {
        if( Settings_handle::st.is("scale-branches") ||
            Settings_handle::st.is("truncate-branches") ||
            Settings_handle::st.is("fixed-branches") )
        {
            cout << "modified guide tree: " << root->print_tree()<<endl;
        }
    }

    if(Settings_handle::st.is("output-graph"))
    {
        fr.write_graph(outfile, root, true);
    }

    if(Settings_handle::st.is("mpost-graph-file")){
        root->write_sequence_graphs();
    }

    if(Settings_handle::st.is("time"))
        cout <<"Time main::exit: "<<double(clock()-t_start)/CLOCKS_PER_SEC << endl;

    if(!Settings_handle::st.is("silent"))
    {
        time_t s_time;
        time( &s_time );
        cout <<endl<< "The analysis finished: " << asctime( localtime( &s_time ) );
        cout<<"Total time used: "<<double_t(clock()-analysis_start_time)/CLOCKS_PER_SEC<<" sec."<<endl<<endl;
    }

    delete root;

    return 1;
}
