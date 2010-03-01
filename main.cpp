#include <string>
#include <vector>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/newick_reader.h"
#include "utils/node.h"
#include "utils/fasta_reader.h"
#include "utils/xml_writer.h"
#include "utils/model_factory.h"
#include "utils/evol_model.h"
#include "main/reads_alignment.h"

using namespace std;

using namespace ppa;

int main(int argc, char *argv[])
{

    /*
    / Input: command line arguments & data
   */

    // Read the arguments
    int rv = Settings_handle::st.read_command_line_arguments(argc, argv);
    
    srand(time(0));
    
    // Read the sequences
    bool gapped_seqs = false;
    Fasta_reader fr;
    vector<Fasta_entry> sequences;
    if(Settings_handle::st.is("seqfile")){

        string seqfile =  Settings_handle::st.get("seqfile").as<string>();
        cout<<"Data file: "<<seqfile<<endl;

        fr.read(seqfile, sequences, true);
    }
    else if(Settings_handle::st.is("cds-seqfile")){

        string seqfile =  Settings_handle::st.get("cds-seqfile").as<string>();
        cout<<"CDS alignment file: "<<seqfile<<endl;

        fr.read(seqfile, sequences, true);
        gapped_seqs = true;
    }
    else
    {
        cout<<endl<<"Error: No sequence file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }

    // Read the guidetree

    Node *root;
    bool tree_ok = false;
    if(Settings_handle::st.is("treefile")){
        string treefile =  Settings_handle::st.get("treefile").as<string>();
        cout<<"Tree file: "<<treefile<<endl;

        Newick_reader nr;
        string tree = nr.read_tree(treefile);
        root = nr.parenthesis_to_tree(tree);

        tree_ok = true;
    }
    else if(Settings_handle::st.is("cds-treefile")){
        string treefile =  Settings_handle::st.get("cds-treefile").as<string>();
        cout<<"CDS tree file: "<<treefile<<endl;

        Newick_reader nr;
        string tree = nr.read_tree(treefile);
        root = nr.parenthesis_to_tree(tree);

        tree_ok = true;
    }
    else if(gapped_seqs && sequences.size()==1)
    {
        root = new Node();
        root->set_name(sequences.at(0).name);
        root->add_name_comment( sequences.at(0).comment );

        string full_char_alphabet = "ACGTRYMKWSBDHVN";

        root->add_sequence( sequences.at(0), full_char_alphabet);

        tree_ok = true;
    }
    else
    {
        cout<<endl<<"Error: No tree file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }

    /*
    / Check that input is fine.
   */

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

    // Get the sequence data type
    //
    int data_type = fr.check_sequence_data_type(sequences);


    // Check that the sequences are fine; also computes the base frequencies.

    Model_factory mf(data_type);

    if(!fr.check_alphabet(mf.get_char_alphabet(),mf.get_full_char_alphabet(),sequences))
        cout<<"\nWarning: Illegal characters in input sequences removed!"<<endl;

    if(data_type==Model_factory::dna)
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

    // Place the sequences to nodes
    // and align them!

    fr.place_sequences_to_nodes(&sequences,&leaf_nodes,mf.get_full_char_alphabet(),gapped_seqs);

    int count = 1;
    root->name_internal_nodes(&count);

    root->start_alignment(&mf);

    // If reads sequences, add them to the alignment

    if( Settings_handle::st.is("readsfile") )
    {
        Reads_alignment ra;
        ra.align(root,&mf,count);

        root = ra.get_global_root();
    }

    // Collect results.
    //
    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,Settings_handle::st.is("output-ancestors"));

    // Save results in output file
    //
    string outfile =  "outfile";

    if(Settings_handle::st.is("outfile")){
        outfile =  Settings_handle::st.get("outfile").as<string>();
    }

    cout<<"Alignment files: "<<outfile<<".fas, "<<outfile<<".xml"<<endl;

    fr.set_chars_by_line(70);
    fr.write(outfile, aligned_sequences, true);

    if(Settings_handle::st.is("output-graph"))
        fr.write_graph(outfile, root, true);

    count = 1;
    root->set_name_ids(&count);

    Xml_writer xw;
    xw.write(outfile, root, aligned_sequences, true);

    if(Settings::noise>1 )
    {
        if( Settings_handle::st.is("scale-branches") ||
            Settings_handle::st.is("truncate-branches") ||
            Settings_handle::st.is("fixed-branches") )
        {
            cout << "modified guide tree: " << root->print_tree()<<endl;
        }
    }

    if(Settings_handle::st.is("mpost-graph-file")){
        root->write_sequence_graphs();
    }

    delete root;
}
