#include "utils/optimal_reference.h"
#include "utils/fasta_reader.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/regex.hpp>

using namespace std;
using namespace ppa;

Optimal_reference::Optimal_reference()
{

}

void Optimal_reference::find_optimal_reference(vector<Fasta_entry> *sequences)
{

    int status = system("exonerate  >/dev/null");
    if( WEXITSTATUS(status) != 1 )
    {
        cout<<"The executable for exonerate not found! Using the first read as the reference!";

        vector<Fasta_entry>::iterator it = sequences->begin();
        it++;
        for(;it!=sequences->end();)
            sequences->erase(it);

        return;
    }

    int r = rand();

    int best_score = -1;
    int most_hits = -1;
    Fasta_entry best_read;
    Fasta_entry most_read;

    vector<Fasta_entry>::iterator fit = sequences->begin();
    for(;fit != sequences->end(); fit++)
    {

        int total_score = 0;
        set<string> total_hits;

        this->write_exonerate_input(&(*fit),sequences,r);

        stringstream command;
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";

        FILE *fpipe;
        if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
        {
            perror("Problems with exonerate pipe.\nExiting.\n");
            exit(1);
        }

        char line[256];
        while ( fgets( line, sizeof line, fpipe))
        {
            hit h;
            bool valid = split_sugar_string(string(line),&h);

            if(valid)
            {
                total_score += h.score;
                total_hits.insert(h.target);
            }
        }
        pclose(fpipe);

        int hits_here = total_hits.size();

        if(total_score > best_score)
        {
            best_score = total_score;
            best_read = *fit;
        }

        if(hits_here > most_hits)
        {
            most_hits = total_hits.size();
            most_read = *fit;
        }


        this->delete_files(r);
    }

    fit = sequences->begin();
    for(;fit!=sequences->end();)
        sequences->erase(fit);

    sequences->push_back(best_read);
//    sequences->push_back(most_read);

//    cout<<endl<<best_read.name<<" "<<best_score<<endl;
//    cout<<most_read.name<<" "<<most_hits<<endl;

}

bool Optimal_reference::find_best_exonerate_hit(Fasta_entry *q, vector<Fasta_entry> *sequences, Fasta_entry *best_read)
{

    int r = rand();

    int best_score = -1;
    string best_name;
    bool best_reversed = false;

    this->write_exonerate_input(q,sequences,r);

//    cout<<"query: "<<q->name<<" "<<q->sequence<<endl;
    stringstream command;
    command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -Q dna -T dna -m affine:local -E 2>&1";

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        perror("Problems with exonerate pipe.\nExiting.\n");
        exit(1);
    }

    char line[256];
    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {
//            cout<<"hit: "<<h.target<<" "<<h.score<<endl;
            if(h.score > best_score)
            {
                best_score = h.score;
                best_name  = h.target;
                best_reversed = (h.q_strand != h.t_strand);

//                cout<<line<<endl;
//                cout<<"new best: best_name "<<best_name<<endl;
            }
        }
    }
    pclose(fpipe);

    if(best_score<0)
        return false;

    this->delete_files(r);

    vector<Fasta_entry>::iterator fit = sequences->begin();
    for(;fit!=sequences->end();fit++)
    {
//        cout<<"search: fit "<<fit->name<<", best "<<best_name<<endl;
        if(fit->name == best_name)
        {
            *best_read = *fit;
            best_read->reversed = best_reversed;

            cout<<"Best "<<best_read->name<<" ("<<best_read->reversed<<") "<<"\n";

//           sequences->erase(fit);

            break;
        }
    }

    return true;
}

void Optimal_reference::find_best_hit(Fasta_entry *query, vector<Fasta_entry> *sequences, Model_factory *mf, Fasta_entry *best_read)
{

    int r = rand();

    float best_score = -1;

    for(int i=0;i<(int)sequences->size();i++)
    {

        if(query->name == sequences->at(i).name)
            continue;

        Node root_node;
        this->copy_node_details(&root_node,query);
        root_node.set_distance_to_parent(0.001);


        Node node;
        node.set_name("#0#");

        Node reads_node;
        this->copy_node_details(&reads_node,&sequences->at(i));

        node.add_left_child(&root_node);
        node.add_right_child(&reads_node);

        node.align_sequences_this_node(mf,true,false);

        float read_identity = this->read_alignment_identity(&node, sequences->at(i).name, root_node.get_name());

        Node node_rc;
        node_rc.set_name("#0#");

        Node reads_node_rc;
        this->copy_node_details(&reads_node_rc,&sequences->at(i),true);

        node_rc.add_left_child(&root_node);
        node_rc.add_right_child(&reads_node_rc);


        node_rc.align_sequences_this_node(mf,true,false);

        float read_identity_rc = this->read_alignment_identity(&node_rc, sequences->at(i).name, root_node.get_name());

//        if(!Settings_handle::st.is("silent"))
//            cout<<sequences->at(i).name<<": forward identity: "<<read_identity<<"; backward identity: "<<read_identity_rc<<endl;


        if(read_identity > read_identity_rc && read_identity > best_score)
        {
            *best_read = sequences->at(i);
            best_read->reversed = false;
            best_score = read_identity;
        }
        else if(read_identity < read_identity_rc && read_identity_rc > best_score)
        {
            *best_read = sequences->at(i);
            best_read->reversed = true;
            best_score = read_identity_rc;
        }

        node.has_left_child(false);
        node.has_right_child(false);

        node_rc.has_left_child(false);
        node_rc.has_right_child(false);
    }

    cout<<best_read->name<<" "<<best_score<<" ("<<best_read->reversed<<")"<<endl;


//    vector<Fasta_entry>::iterator fit = sequences->begin();
//    for(;fit!=sequences->end();fit++)
//    {
////        cout<<"search: fit "<<fit->name<<", best "<<best_read->name<<endl;
//        if(fit->name == best_read->name)
//        {
////            cout<<"delete "<<best_read->name<<"\n";
//            sequences->erase(fit);

//            break;
//        }
//    }

}

void Optimal_reference::align(Node *root, Model_factory *mf, int count)
{

    int status = system("exonerate  >/dev/null");
    if( WEXITSTATUS(status) != 1 )
    {
        cout<<"The executable for exonerate not found! Exiting!";
        exit(0);
    }


    string file = Settings_handle::st.get("readsfile").as<string>();

    Fasta_reader fr;
    vector<Fasta_entry> reads;
    cout<<"Reads data file: "<<file<<endl;

    try
    {
        fr.read(file, reads, true);
    }
    catch (ppa::IOException& e) {
        cout<<"Error reading the reads file '"<<file<<"'.\nExiting.\n\n";
        exit(1);
    }


    int data_type = fr.check_sequence_data_type(&reads);

    if(!fr.check_alphabet(&reads,data_type))
        if(!Settings_handle::st.is("silent"))
            cout<<"\nWarning: Illegal characters in input reads sequences removed!"<<endl;


    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;


    // Trim read ends
    //
    if(Settings_handle::st.is("trim-read-ends"))
    {
        fr.trim_fastq_reads(&reads);

        vector<Fasta_entry>::iterator fit1 = reads.begin();

        for(;fit1 != reads.end();fit1++)
        {
            stringstream trimming;
            trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end;
            fit1->comment += trimming.str();
        }
    }



    global_root = root;

    int remaining = reads.size();

    bool use_exonerate = false;

    while(remaining>0)
    {

        cout<<"remaining: "<<remaining<<endl;

        Fasta_entry best_read;

        Fasta_entry root_sequence;
        global_root->get_node_sequence(&root_sequence);


        bool exonarate_hit = false;
        if(use_exonerate)
            exonarate_hit = this->find_best_exonerate_hit(&root_sequence,&reads,&best_read);

        if( exonarate_hit )
            cout<<"exonerate hit: "<<best_read.name<<endl;
        else
        {
            this->find_best_hit(&root_sequence,&reads,mf,&best_read);
            cout<<"pagan hit: "<<best_read.name<<endl;
            use_exonerate = true;
        }

        Node * node = new Node();

        stringstream ss;
        ss<<"#"<<count<<"#";
        node->set_name(ss.str());

        global_root->set_distance_to_parent(0.001);
        node->add_left_child(global_root);


        Node * reads_node = new Node();
        this->copy_node_details(reads_node,&best_read,best_read.reversed);
        node->add_right_child(reads_node);


        node->align_sequences_this_node(mf,true,false);

        float read_overlap;
        float read_identity;
        this->read_alignment_scores(node, best_read.name, global_root->get_name(), &read_overlap, &read_identity);

        if(!Settings_handle::st.is("silent"))
            cout<<"overlap: "<<read_overlap<<", identity "<<read_identity<<endl;

        float min_overlap = Settings_handle::st.get("min-reads-overlap").as<float>();
        float min_identity = Settings_handle::st.get("min-reads-identity").as<float>();

        if(read_overlap < min_overlap || read_identity < min_identity)
        {

            node->has_left_child(false);
            delete node;

        } else {

            count++;
            global_root = node;

        }

        vector<Fasta_entry>::iterator fit = reads.begin();
        for(;fit!=reads.end();fit++)
        {
            if(fit->name == best_read.name)
            {
                cout<<"delete "<<best_read.name<<"\n";
                reads.erase(fit);

                break;
            }
        }

        remaining = reads.size();
    }

}

void Optimal_reference::copy_node_details(Node *reads_node,Fasta_entry *read,bool turn_revcomp)
{
    double r_dist = Settings_handle::st.get("reads-distance").as<float>();

    reads_node->set_distance_to_parent(r_dist);
    reads_node->set_name(read->name);
    reads_node->add_name_comment(read->comment);
    reads_node->add_sequence( *read, read->data_type, false, true, turn_revcomp);
    reads_node->get_sequence()->is_read_sequence(true);

}

float Optimal_reference::read_alignment_overlap(Node * node, string read_name, string ref_node_name)
{
    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;

    for( int j=0; j < node_sequence->sites_length(); j++ )
    {
        bool read_has_site = node->has_site_at_alignment_column(j,read_name);
        bool any_other_has_site = node->any_other_has_site_at_alignment_column(j,read_name);

        if(read_has_site)
            read_length++;

        if(read_has_site && any_other_has_site)
        {
            aligned++;
        }

    }

    if(Settings::noise>0)
        cout<<"  aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"]"<<endl;

    return (float)aligned/(float)read_length;
}

float Optimal_reference::read_alignment_identity(Node * node, string read_name, string ref_node_name)
{
    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;
    int matched = 0;

    for( int j=0; j < node_sequence->sites_length(); j++ )
    {
        bool read_has_site = node->has_site_at_alignment_column(j,read_name);
        bool any_other_has_site = node->any_other_has_site_at_alignment_column(j,read_name);

        if(read_has_site)
            read_length++;

        if(read_has_site && any_other_has_site)
        {
            aligned++;

            int state_read = node->get_state_at_alignment_column(j,read_name);
            int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);
            if(state_read == state_ref)
                matched++;
        }


    }

    if(Settings::noise>0)
        cout<<"  matched positions "<<(float)matched/(float)read_length<<" ["<<matched<<"/"<<read_length<<"]"<<endl;

    return (float)matched/(float)read_length;
}

void Optimal_reference::read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity)
{
    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;
    int matched = 0;

    for( int j=0; j < node_sequence->sites_length(); j++ )
    {
        bool read_has_site = node->has_site_at_alignment_column(j,read_name);
        bool any_other_has_site = node->any_other_has_site_at_alignment_column(j,read_name);

        if(read_has_site)
            read_length++;

        if(read_has_site && any_other_has_site)
        {
            aligned++;

            int state_read = node->get_state_at_alignment_column(j,read_name);
            int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);
            if(state_read == state_ref)
                matched++;
        }


    }

    if(Settings::noise>0)
    {
        cout<<"  aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"]"<<endl;
        cout<<"  matched positions "<<(float)matched/(float)aligned<<" ["<<matched<<"/"<<aligned<<"]"<<endl;
    }

    *overlap  = (float)aligned/(float)read_length;
    *identity = (float)matched/(float)aligned;

}

bool Optimal_reference::split_sugar_string(const string& row,hit *h)
{

    const boost::regex pattern("sugar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->target     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );
    }

    return valid;
}


void Optimal_reference::write_exonerate_input(Fasta_entry *q, vector<Fasta_entry> *seqs, int r)
{

    // create exonerate input
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    q_output<<">"<<q->name<<endl<<q->sequence<<endl;
    q_output.close();

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    vector<Fasta_entry>::iterator it = seqs->begin();
    for(;it!=seqs->end();it++)
    {
        if(it->name != q->name)
            t_output<<">"<<it->name<<endl<<it->sequence<<endl;
    }
    t_output.close();

    return;
}

void Optimal_reference::write_exonerate_input(string q, vector<Fasta_entry> *seqs, int r)
{

    // create exonerate input
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    q_output<<">query"<<endl<<q<<endl;
    q_output.close();

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    vector<Fasta_entry>::iterator it = seqs->begin();
    for(;it!=seqs->end();it++)
    {
        t_output<<">"<<it->name<<endl<<it->sequence<<endl;
    }
    t_output.close();

    return;
}

void Optimal_reference::delete_files(int r)
{
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    if( remove( q_name.str().c_str() ) != 0 )
       perror( "Error deleting file" );
    if( remove( t_name.str().c_str() ) != 0 )
       perror( "Error deleting file" );
}
