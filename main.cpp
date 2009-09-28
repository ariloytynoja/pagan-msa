#include <string>
#include <vector>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/newick_reader.h"
#include "utils/node.h"
#include "utils/fasta_reader.h"
#include "utils/xml_writer.h"
#include "utils/model_factory.h"
#include "utils/dna_model.h"

using namespace std;

using namespace ppa;

int main(int argc, char *argv[])
{

    /*
    / Input: command line arguments & data
   */

    // Read the arguments

//    Settings st;
//    Settings_handle handle;
    int rv = Settings_handle::st.read_command_line_arguments(argc, argv);

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
    else
    {
        cout<<endl<<"Error: No tree file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }


    // Read the sequences

    Fasta_reader fr;
    vector<Fasta_entry> sequences;
    bool seqs_ok = false;
    if(Settings_handle::st.is("seqfile")){

        string seqfile =  Settings_handle::st.get("seqfile").as<string>();
        cout<<"Data file: "<<seqfile<<endl;

        fr.read(seqfile, sequences, true);

        seqs_ok = true;
    }
    else
    {
        cout<<endl<<"Error: No sequence file defined."<<endl;
        Settings_handle::st.info();

        exit(0);
    }


    /*
    / Check that input is fine.
   */

    // Check that the guidetree and sequences match

    vector<Node*> leaf_nodes;
    root->get_leaf_nodes(&leaf_nodes);
    fr.check_sequence_names(&sequences,&leaf_nodes,&Settings_handle::st);


    // Check that the sequences are fine; also computes the base frequencies.

    Model_factory mf;

    if(!fr.check_alphabet(mf.get_dna_alphabet(),mf.get_full_dna_alphabet(),sequences))
        cout<<"\nWarning: Illegal characters in input sequences removed!"<<endl;

    float *dna_pi = fr.base_frequencies();


    // Create a DNA alignment model using empirical base frequencies.

    mf.dna_model(dna_pi,&Settings_handle::st);


    // Place the sequences to nodes
    // and align them!

    fr.place_sequences_to_nodes(&sequences,&leaf_nodes,mf.get_full_dna_alphabet());

    int count = 1;
    root->name_internal_nodes(&count);

    root->start_alignment(&mf);


    // Collect results.

    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,Settings_handle::st.is("output-ancestors"));


    if(1) {
        if(Settings_handle::st.is("outfile")){
            string outfile =  Settings_handle::st.get("outfile").as<string>();
            cout<<"Alignment files: "<<outfile<<".fas, "<<outfile<<".xml"<<endl;

            Fasta_reader fr;
            fr.set_chars_by_line(70);
            fr.write(outfile, aligned_sequences, true);


            int count = 1;
            root->set_name_ids(&count);

            Xml_writer xw;
            xw.write(outfile, root, aligned_sequences, true);

        }
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

    if(Settings_handle::st.is("mpost-graph-file")){
        root->write_sequence_graphs();
    }

}
/*

  Make a factory that gioves models.
  Both DNA and codon.
  Many?
*/
