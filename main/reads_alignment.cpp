#include "reads_alignment.h"
#include <sstream>

using namespace std;
using namespace ppa;

Reads_alignment::Reads_alignment(){}

void Reads_alignment::align(Node *root, Model_factory *mf, int count)
{
    string file = Settings_handle::st.get("readsfile").as<string>();

    Fasta_reader fr;
    vector<Fasta_entry> reads;
    cout<<"Reads data file: "<<file<<endl;

    fr.read(file, reads, true);

    if(!fr.check_alphabet(mf->get_char_alphabet(),mf->get_full_char_alphabet(),reads))
        cout<<"\nWarning: Illegal characters in input reads sequences removed!"<<endl;

    double r_dist = Settings_handle::st.get("reads-distance").as<float>();
    global_root = root;
    for(int i=0;i<reads.size();i++)
    {
        cout<<"Aligning read ">>i+1<<"/"<<reads.size()<<endl;

        Node * node = new Node();

        stringstream ss;
        ss<<"#"<<count<<"#";
        node->set_name(ss.str());

        global_root->set_distance_to_parent(0.001);
        node->add_left_child(global_root);

        Node * reads_node = new Node();
        reads_node->set_distance_to_parent(r_dist);
        reads_node->set_name(reads.at(i).name);
        reads_node->add_name_comment( reads.at(i).comment );
        reads_node->add_sequence( reads.at(i).sequence, mf->get_full_char_alphabet());
        node->add_right_child(reads_node);

        node->align_sequences_this_node(mf,true);

        count++;
        global_root = node;
    }
}
