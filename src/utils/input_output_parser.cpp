/***************************************************************************
 *   Copyright (C) 2012 by Ari Loytynoja                                   *
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

#include "input_output_parser.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"
#include "utils/newick_reader.h"
#include "utils/fasta_reader.h"
#include "utils/xml_writer.h"
#include "utils/model_factory.h"
#include "utils/evol_model.h"
#include "utils/optimal_reference.h"
#include "main/node.h"
#include "main/reads_aligner.h"
#include "utils/substring_hit.h"
#include <iostream>
#include <vector>

using namespace ppa;
using namespace std;


Input_output_parser::Input_output_parser()
{
}

/************************************************************************************/

void Input_output_parser::parse_input_sequences(vector<Fasta_entry> *sequences, bool *reference_alignment)
{

    /***********************************************************************/
    /*  Overlapping paired-end read merge only                             */
    /***********************************************************************/

    if (Settings_handle::st.is("overlap-merge-only")) {

        if (!Settings_handle::st.is("queryfile")) {
            Log_output::write_out("No reads file given. Exiting.\n\n",0);
            exit(1);
        }

        Reads_aligner ra;
        ra.merge_reads_only();

        exit(0);
    }


    /***********************************************************************/
    /*  Read the sequences                                                 */
    /***********************************************************************/

    Fasta_reader fr;

    if(Settings_handle::st.is("seqfile"))
    {
        string seqfile =  Settings_handle::st.get("seqfile").as<string>();
        Log_output::write_out("Sequence file: "+seqfile+"\n",1);

        try
        {
            fr.read(seqfile, *sequences, true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the sequence file '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }
    }
    else if(Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("ref-seqfile").as<string>();
        Log_output::write_out("Reference alignment file: "+seqfile+"\n",1);

        try
        {
            fr.read(seqfile, *sequences, true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the reference alignment  file '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        *reference_alignment = true;
    }
    else if(Settings_handle::st.is("queryfile") && Settings_handle::st.is("find-cluster-reference") && !Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("queryfile").as<string>();
        Log_output::write_out("Optimal reference sequence from: "+seqfile+"\n",1);

        try
        {
            fr.read(seqfile, *sequences, true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the queryfile '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        Optimal_reference ore;
        ore.find_optimal_reference(sequences);

        *reference_alignment = true;

    }
    else if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
    {
        string seqfile =  Settings_handle::st.get("queryfile").as<string>();
        Log_output::write_out("Reference sequence from: "+seqfile+"\n",1);

        try
        {
            fr.read(seqfile, *sequences, true);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the queryfile '"+seqfile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        vector<Fasta_entry>::iterator it = sequences->begin();
        it++;
        for(;it!=sequences->end();)
            sequences->erase(it);

       *reference_alignment = true;
    }
    else
    {
        Log_output::write_out("\nError: No sequence file defined.\n",0);
        Settings_handle::st.info();

        exit(1);
    }
}

/************************************************************************************/

Node *Input_output_parser::parse_input_tree(vector<Fasta_entry> *sequences,bool reference_alignment)
{
    /***********************************************************************/
    /*  Read the guidetree                                                 */
    /***********************************************************************/

    Fasta_reader fr;
    Node *root;
    if(Settings_handle::st.is("treefile"))
    {
        string treefile =  Settings_handle::st.get("treefile").as<string>();
        Log_output::write_out("Tree file: "+treefile+"\n",1);

        Newick_reader nr;
        string tree;
        try
        {
            tree = nr.read_tree(treefile);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the guide tree file '"+treefile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        root = nr.parenthesis_to_tree(tree);
    }
    else if(Settings_handle::st.is("ref-treefile"))
    {
        string treefile =  Settings_handle::st.get("ref-treefile").as<string>();
        Log_output::write_out("Reference tree file: "+treefile+"\n",1);

        Newick_reader nr;
        string tree;
        try
        {
            tree = nr.read_tree(treefile);
        }
        catch (ppa::IOException& e) {
            Log_output::write_out("Error reading the reference tree file '"+treefile+"'.\nExiting.\n\n",0);
            exit(1);
        }

        root = nr.parenthesis_to_tree(tree);
    }
    else if(reference_alignment && sequences->size()==1)
    {
        int data_type = fr.check_sequence_data_type(sequences);

        root = new Node();
        root->set_name(sequences->at(0).name);
        root->add_name_comment( sequences->at(0).comment );
        root->set_distance_to_parent(0);

        root->add_sequence( sequences->at(0), data_type);
    }
    else
    {
        Log_output::write_out("\nError: No tree file defined.\n",0);
        Settings_handle::st.info();

        exit(1);
    }

    return root;
}

/************************************************************************************/

void Input_output_parser::match_sequences_and_tree(std::vector<Fasta_entry> *sequences,Node *root,bool reference_alignment,int *data_type)
{
    /***********************************************************************/
    /*  Check that input is fine and place the sequences to nodes          */
    /***********************************************************************/

    Fasta_reader fr;

    // Check that the guidetree and sequences match

    vector<Node*> leaf_nodes;
    root->get_leaf_nodes(&leaf_nodes);


    bool tree_branches_ok = fr.check_sequence_names(sequences,&leaf_nodes);
    if(!tree_branches_ok)
    {
        Log_output::write_out("Attempting to prune the tree: ",2);
        Log_output::flush();

        root->prune_tree();
        Log_output::write_out("pruning done.\n",2);
        Log_output::write_out("New tree:\n"+root->print_tree()+"\n\n",3);
    }

    leaf_nodes.clear();
    root->get_leaf_nodes(&leaf_nodes);


    // Check the data type

    *data_type = fr.check_sequence_data_type(sequences);

    if(!fr.check_alphabet(sequences,*data_type))
        Log_output::write_out(" Warning: Illegal characters in input sequences removed!\n",2);


    //  Place the sequences to nodes

    fr.place_sequences_to_nodes(sequences,&leaf_nodes,reference_alignment,*data_type);


}

/************************************************************************************/

void Input_output_parser::define_alignment_model(Model_factory *mf,int data_type)
{
    /***********************************************************************/
    /*  Define the alignment model                                         */
    /***********************************************************************/

    Fasta_reader fr;

    if(data_type==Model_factory::dna && Settings_handle::st.is("codons"))
    {
        // Create a codon alignment model using KHG.
        Log_output::write_out("Model_factory: creating a codon model\n",3);
        mf->codon_model(&Settings_handle::st); // does it need the handle????
    }
    else if(data_type==Model_factory::dna)
    {
        // Create a DNA alignment model using empirical base frequencies.
        Log_output::write_out("Model_factory: creating a DNA model\n",3);
        float *dna_pi = fr.base_frequencies();
        mf->dna_model(dna_pi,&Settings_handle::st);
    }
    else if(data_type==Model_factory::protein)
    {
        // Create a protein alignment model using WAG.
        Log_output::write_out("Model_factory: creating a protein model\n",3);
        mf->protein_model(&Settings_handle::st); // does it need the handle????
    }

}

/************************************************************************************/

void Input_output_parser::output_aligned_sequences(std::vector<Fasta_entry> *sequences,Node *root)
{

    /***********************************************************************/
    /*  Collect the results and output them                                */
    /***********************************************************************/

//    root->show_seqs();

    Fasta_reader fr;

    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,Settings_handle::st.is("output-ancestors"));

    Log_output::clean_output();


    // See if any sequences were placed
    //
    if(Settings_handle::st.is("queryfile") && (int)sequences->size() == (int)aligned_sequences.size())
    {
        Log_output::write_out("Failed to extend the alignment. No output created.\n",0);
    }
    else
    {

        // Save results in output file
        //
        string outfile =  "outfile";

        if(Settings_handle::st.is("outfile"))
            outfile =  Settings_handle::st.get("outfile").as<string>();

        string format = "fasta";
        if(Settings_handle::st.is("outformat"))
            format = Settings_handle::st.get("outformat").as<string>();

        if(Settings_handle::st.is("xml"))
            Log_output::write_out("Alignment files: "+outfile+fr.get_format_suffix(format)+", "+outfile+".xml\n",0);
        else
            Log_output::write_out("Alignment file: "+outfile+fr.get_format_suffix(format)+"\n",0);

        fr.set_chars_by_line(70);
        fr.write(outfile, aligned_sequences, format, true);

        int count = 1;
        root->set_name_ids(&count);

        if(Settings_handle::st.is("xml"))
        {
            Xml_writer xw;
            xw.write(outfile, root, aligned_sequences, true);
        }

        if(Settings_handle::st.is("build-contigs"))
        {
            vector<Fasta_entry> contigs;
            root->reconstruct_contigs(&contigs,false);

            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append("_contigs");
            Log_output::write_out("Contig file: "+outfile+".fas\n",1);

            fr.set_chars_by_line(70);
            fr.write(outfile, contigs, "fasta", true);
        }

        if(Settings_handle::st.is("output-consensus"))
        {
            vector<Fasta_entry> contigs;
            root->reconstruct_contigs(&contigs,false,true);

            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append("_consensus");
            Log_output::write_out("Consensus file: "+outfile+".fas\n",1);

            fr.remove_gap_only_columns(&contigs);

            fr.set_chars_by_line(70);
            fr.write(outfile, contigs, "fasta", true);
        }

        if(Settings_handle::st.is("translate") || Settings_handle::st.is("mt-translate") || Settings_handle::st.is("find-best-orf"))
        {
            string outfile =  "outfile";
            if(Settings_handle::st.is("outfile"))
                outfile =  Settings_handle::st.get("outfile").as<string>();

            outfile.append(".dna");
            Log_output::write_out("Back-translated alignment file: "+outfile+".fas\n",0);

            fr.set_chars_by_line(70);
            fr.write_dna(outfile, aligned_sequences, *sequences,root,true,Fasta_reader::plain_alignment);

            if(Settings_handle::st.is("build-contigs"))
            {
                string outfile =  "outfile";
                if(Settings_handle::st.is("outfile"))
                    outfile =  Settings_handle::st.get("outfile").as<string>();

                outfile.append("_contigs.dna");
                Log_output::write_out("Back-translated contig file: "+outfile+".fas\n",1);

                fr.set_chars_by_line(70);
                fr.write_dna(outfile, aligned_sequences, *sequences,root,true,Fasta_reader::contig_alignment);
            }
            if(Settings_handle::st.is("output-consensus"))
            {
                string outfile =  "outfile";
                if(Settings_handle::st.is("outfile"))
                    outfile =  Settings_handle::st.get("outfile").as<string>();

                outfile.append("_consensus.dna");
                Log_output::write_out("Back-translated consensus file: "+outfile+".fas\n",1);

                fr.set_chars_by_line(70);
                fr.write_dna(outfile, aligned_sequences, *sequences,root,true,Fasta_reader::consensus_only);
            }
        }

        if(Settings_handle::st.is("output-ancestors"))
        {
            fr.write_anctree(outfile, root);
        }

        if(Settings_handle::st.is("output-nhx-tree"))
        {
            root->write_nhx_tree(outfile);
        }

        if( Settings_handle::st.is("scale-branches") ||
             Settings_handle::st.is("truncate-branches") ||
              Settings_handle::st.is("fixed-branches") )
        {
            Log_output::write_out("Modified guide tree: " +root->print_tree()+"\n",2);
        }

        if(Settings_handle::st.is("output-graph"))
        {
            fr.write_graph(outfile, root, true);
        }

        if(Settings_handle::st.is("mpost-graph-file")){
            root->write_sequence_graphs();
        }
    }

}
