#include "reads_alignment.h"
#include <sstream>
#include <fstream>

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


    if( Settings_handle::st.is("pair-end") )
    {
        this->find_paired_reads( &reads );
    }
    else
    {
        this->add_trimming_comment( &reads );
    }


    if(Settings_handle::st.is("align-reads-at-root"))
    {
        if(Settings_handle::st.is("discard-pairwise-overlapping-reads"))
        {
            this->remove_overlapping_reads(&reads, mf);
        }

        string ref_root_name = root->get_name();

        global_root = root;
        for(int i=0;i<(int)reads.size();i++)
        {
            cout<<"read "<<i+1<<"/"<<reads.size()<<"; ";

            Node * node = new Node();

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            global_root->set_distance_to_parent(0.001);
            node->add_left_child(global_root);

            Node * reads_node = new Node();
            this->copy_node_details(reads_node,&reads.at(i), mf->get_full_char_alphabet());

            node->add_right_child(reads_node);

            node->set_nhx_tid(node->get_left_child()->get_nhx_tid());
            node->get_right_child()->set_nhx_tid(node->get_left_child()->get_nhx_tid());

            node->align_sequences_this_node(mf,true);

            // check if the alignment significantly overlaps with the reference alignment
            //
            bool read_overlaps = this->read_alignment_overlaps(node, reads.at(i).name, ref_root_name);

            if(read_overlaps)
            {
                count++;
                global_root = node;
            }
            // else delete the node; do not use the read
            else
            {
                node->has_left_child(false);
                delete node;
            }

        }
    }
    else
    {

        global_root = root;

        // vector of node names giving the best node for each read
        //
        this->find_nodes_for_reads(root, &reads, mf);

        if(Settings_handle::st.is("placement-only"))
            exit(0);

        set<string> unique_nodes;
        for(int i=0;i<(int)reads.size();i++)
            unique_nodes.insert(reads.at(i).node_to_align);

        if(unique_nodes.find("discarded_read") != unique_nodes.end())
            unique_nodes.erase(unique_nodes.find("discarded_read"));

        map<string,Node*> nodes_map;
        root->get_internal_nodes(&nodes_map);

        // do one tagged node at time
        //
        for(set<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
        {
            vector<Fasta_entry> reads_for_this;

            for(int i=0;i<(int)reads.size();i++)
                if(reads.at(i).node_to_align == *sit)
                    reads_for_this.push_back(reads.at(i));

            this->sort_reads_vector(&reads_for_this);

            // remove fully overlapping reads that are mapped to this node
            //
            if(Settings_handle::st.is("rank-reads-for-nodes"))
            {
                if(Settings_handle::st.is("discard-overlapping-identical-reads"))
                {
                    this->remove_target_overlapping_identical_reads(&reads_for_this,mf);

                    cout<<"After removing overlapping ones, for node "<<*sit<<" reads remaining:\n";
                    for(int i=0;i<(int)reads_for_this.size();i++)
                        cout<<" "<<reads_for_this.at(i).name<<" "<<reads_for_this.at(i).node_score<<endl;
                    cout<<endl;
                }
                else if(Settings_handle::st.is("discard-overlapping-reads"))
                {
                    this->remove_target_overlapping_reads(&reads_for_this);

                    cout<<"After removing overlapping ones, for node "<<*sit<<" reads remaining:\n";
                    for(int i=0;i<(int)reads_for_this.size();i++)
                        cout<<" "<<reads_for_this.at(i).name<<" "<<reads_for_this.at(i).node_score<<endl;
                    cout<<endl;
                }
                else if(Settings_handle::st.is("discard-pairwise-overlapping-reads"))
                {
                    this->remove_overlapping_reads(&reads_for_this,mf);

                    cout<<"After removing overlapping ones, for node "<<*sit<<" reads remaining:\n";
                    for(int i=0;i<(int)reads_for_this.size();i++)
                        cout<<" "<<reads_for_this.at(i).name<<" "<<reads_for_this.at(i).node_score<<endl;
                    cout<<endl;
                }
            }
            else
            {
                if( Settings_handle::st.is("discard-overlapping-identical-reads") ||
                         Settings_handle::st.is("discard-overlapping-reads") )
                {
                    cout<<"\nWarning: without ranking the reads for nodes, one cannot resolve overlap between reads. The flag has no effect!\n\n";
                }
            }

            string ref_node_name = *sit;

            Node *current_root = nodes_map.find(ref_node_name)->second;
            double orig_dist = current_root->get_distance_to_parent();

            bool alignment_done = false;
            int alignments_done = 0;

            // align the remaining reads to this node
            //
            for(int i=0;i<(int)reads_for_this.size();i++)
            {
                cout<<"read "<<i+1<<"/"<<reads_for_this.size()<<"; ";

                Node * node = new Node();

                stringstream ss;
                ss<<"#"<<count<<"#";
                node->set_name(ss.str());

                current_root->set_distance_to_parent(0.001);
                node->add_left_child(current_root);

                Node * reads_node = new Node();
                this->copy_node_details(reads_node,&reads_for_this.at(i), mf->get_full_char_alphabet());

                node->add_right_child(reads_node);

                node->set_nhx_tid(node->get_left_child()->get_nhx_tid());
                node->get_right_child()->set_nhx_tid(node->get_left_child()->get_nhx_tid());

                node->align_sequences_this_node(mf,true);

//                node->get_sequence()->print_sequence();

                // check if the alignment significantly overlaps with the reference alignment
                //
                bool read_overlaps = this->read_alignment_overlaps(node, reads_for_this.at(i).name, ref_node_name);

                if(read_overlaps)
                {
                    count++;
                    current_root = node;

                    if( orig_dist > current_root->get_distance_to_parent() )
                        orig_dist -= current_root->get_distance_to_parent();

                    alignment_done = true;
                    alignments_done++;
                }
                // else delete the node; do not use the read
                else
                {
                    node->has_left_child(false);
                    delete node;
                }

            }

            current_root->set_distance_to_parent(orig_dist);


            if(alignment_done)
            {
                bool parent_found = this->correct_sites_index(current_root, ref_node_name, alignments_done, &nodes_map);

                if(!parent_found)
                {
                    global_root = current_root;
                }
            }
        }


    }
}

void Reads_alignment::add_trimming_comment(vector<Fasta_entry> *reads)
{
    vector<Fasta_entry>::iterator fit1 = reads->begin();

    for(;fit1 != reads->end();fit1++)
    {
        stringstream trimming;
        trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end;
        fit1->comment += trimming.str();
//        cout<<"s: "<<trimming.str()<<endl;
    }
}


void Reads_alignment::find_paired_reads(vector<Fasta_entry> *reads)
{

    vector<Fasta_entry>::iterator fit1 = reads->begin();

    for(;fit1 != reads->end();fit1++)
    {
        string name1 = fit1->name;
        if(name1.substr(name1.length()-2) != "/1")
        {
            continue;
        }

        vector<Fasta_entry>::iterator fit2 = fit1;
        fit2++;

        for(;fit2 != reads->end();)
        {
            string name2 = fit2->name;
            if(name2.substr(name2.length()-2) != "/2")
            {
                fit2++;
                continue;
            }

            if(name1.substr(0,name1.length()-2)==name2.substr(0,name2.length()-2))
            {
                fit1->name = name1.substr(0,name1.length()-1)+"p12";
                fit1->comment = fit2->comment;
                fit1->first_read_length = fit1->sequence.length();
                fit1->sequence += "0"+fit2->sequence;
                fit1->quality += "0"+fit2->quality;

                stringstream trimming;
                trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end<<":P2ST"<<fit2->trim_start<<":P2ET"<<fit2->trim_end;
                fit1->comment += trimming.str();

                reads->erase(fit2);

                continue;
            }

            fit2++;
        }
    }


    for(fit1 = reads->begin();fit1 != reads->end();fit1++)
    {
        if(fit1->comment.find("P1ST")==string::npos)
        {
            stringstream trimming;
            trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end;
            fit1->comment += trimming.str();
        }
    }
}

void Reads_alignment::copy_node_details(Node *reads_node,Fasta_entry *read, string full_alpha)
{
    double r_dist = Settings_handle::st.get("reads-distance").as<float>();

    reads_node->set_distance_to_parent(r_dist);
    reads_node->set_name(read->name);
    reads_node->add_name_comment(read->comment);
    reads_node->add_sequence( *read, full_alpha);

    cout<<"aligning read "<<read->name<<": "<<read->comment<<endl;
}

bool Reads_alignment::read_alignment_overlaps(Node * node, string read_name, string ref_node_name)
{
    float min_overlap = Settings_handle::st.get("min-reads-overlap").as<float>();

    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;

    for( int j=0; j < node_sequence->sites_length(); j++ )
    {
        bool read_has_site = node->has_site_at_alignment_column(j,read_name);
        bool ref_root_has_site = node->has_site_at_alignment_column(j,ref_node_name);

        if(read_has_site)
            read_length++;

        if(read_has_site && ref_root_has_site)
            aligned++;
    }

    cout<<"  alignment score "<<(float)aligned/(float)read_length<<" (simple identity score; needs improving!)"<<endl;

    if((float)aligned/(float)read_length >= min_overlap)
    {
        return true;
    }
    else
    {

        cout<<"Warning: read "<<read_name<<" dropped using the minimum overlap cut-off of "<<min_overlap<<"."<<endl;

        return false;
    }
}


bool Reads_alignment::correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map)
{

    // correct the sites index at the parent node; insertions corrected later
    //
    vector<int> sites_index;
    int index_delta = 0;

    for(int j=0; j<current_root->get_sequence()->sites_length(); j++)
    {
        if(current_root->has_site_at_alignment_column(j,ref_node_name))
        {
            sites_index.push_back(index_delta);
            index_delta = 0;
        }
        else
            index_delta++;
    }

    Node *current_parent = 0;
    map<string,Node*>::iterator mit = nodes_map->begin();
    bool parent_found = false;

    int is_left_child = true;
    for(;mit != nodes_map->end();mit++)
    {
        if(mit->second->get_left_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_left_child(current_root);
            parent_found = true;
        }

        if(mit->second->get_right_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_right_child(current_root);
            is_left_child = false;
            parent_found = true;
        }
    }


    if(parent_found)
    {
        cout<<" Parent of "<<ref_node_name<<" is "<<current_parent->get_name();
        cout<<"; "<<alignments_done<<" alignments done.";//<<endl;

        Sequence *parent_sequence = current_parent->get_sequence();

        index_delta = 0;

        for(int j=0; j<parent_sequence->sites_length(); j++)
        {
            Site *parent_site = parent_sequence->get_site_at(j);

            if(is_left_child && parent_site->get_children()->left_index >= 0)
            {
                index_delta += sites_index.at(parent_site->get_children()->left_index);
                parent_site->get_children()->left_index += index_delta;
            }
            else if(!is_left_child && parent_site->get_children()->right_index >= 0)
            {
                index_delta += sites_index.at(parent_site->get_children()->right_index);
                parent_site->get_children()->right_index += index_delta;
            }
        }

        if(index_delta>0)
            cout<<" Site index needs correcting.\n";
        else
            cout<<" Site index not changed.\n";

        if(index_delta>0)
        {
            if(is_left_child)
                current_parent->left_needs_correcting_sequence_site_index(true);
            else
                current_parent->right_needs_correcting_sequence_site_index(true);
        }

        return true;

    } // if(parent_found)
    else
    {
        cout<<" No parent for "<<ref_node_name<<" found. Assuming that this is root.\n";
        cout<<" "<<alignments_done<<" alignments done.\n";//<<endl;

        return false;
    }

}

void Reads_alignment::find_nodes_for_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{

    cout<<"Find nodes for reads.\n";
    multimap<string,string> tid_nodes;
    root->get_internal_node_names_with_tid_tag(&tid_nodes);

    ofstream pl_output;

    if(Settings_handle::st.is("placement-file"))
    {
        string fname = Settings_handle::st.get("placement-file").as<string>();
       pl_output.open(fname.append(".tsv").c_str());
    }
    for(int i=0;i<(int)reads->size();i++)
    {
        reads->at(i).node_score = -1.0;

        string tid = reads->at(i).tid;
        if(tid != "")
        {
            int matches = tid_nodes.count(tid);
            if(matches == 0)
            {
                cout<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
                reads->at(i).node_to_align = root->get_name();
            }
            else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
            {
                multimap<string,string>::iterator tit = tid_nodes.find(tid);
                if(tit != tid_nodes.end())
                {
                    cout<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                    reads->at(i).node_to_align = tit->second;
                }
            }
            else
            {
                double best_score = -HUGE_VAL;
                string best_node = root->get_name();

                map<string,Node*> nodes;
                root->get_internal_nodes(&nodes);

                multimap<string,string>::iterator tit = tid_nodes.find(tid);
                if(tit != tid_nodes.end())
                {
                    cout<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" matches the nodes:\n";
                    while(tit != tid_nodes.end())
                    {
                        map<string,Node*>::iterator nit = nodes.find(tit->second);

                        double score = this->read_match_score( nit->second, &reads->at(i), mf, best_score);

                        cout<<"   "<<tit->second<<" with score "<<score<<" (simple p-distance; needs improving!)\n";
                        if(score>best_score)
                        {
                            best_score = score;
                            best_node = tit->second;
                        }
                        tit++;
                    }
                }
                if(best_score<0.05)
                {
                    if(Settings_handle::st.is("align-bad-reads-at-root"))
                    {
                        cout<<"Best node aligns with less than 5% of identical sites. Aligning to root instead.\n";
                        reads->at(i).node_to_align = root->get_name();
                    }
                    else
                    {
                        cout<<"Best node aligns with less than 5% of identical sites. Read is discarded.\n";
                        reads->at(i).node_to_align = "discarded_read";
                    }
                }
                else
                {
                    cout<<"Best node "<<best_node<<" (score "<<best_score<<").\n";
                    reads->at(i).node_score = best_score;
                    reads->at(i).node_to_align = best_node;

                    if(Settings_handle::st.is("placement-file"))
                    {
                        pl_output<<reads->at(i).name<<" "<<best_node<<endl;
                    }

                }
            }
        }
        else
        {
            cout<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") has no tid. Aligned to root.\n";
            reads->at(i).node_to_align = root->get_name();
        }
    }

    if(Settings_handle::st.is("placement-file"))
    {
        pl_output.close();
    }
}

double Reads_alignment::read_match_score(Node *node, Fasta_entry *read, Model_factory *mf, float best_score)
{

    double r_dist = Settings_handle::st.get("reads-distance").as<float>();

    double org_dist = node->get_distance_to_parent();
    node->set_distance_to_parent(0.001);

    Node * reads_node1 = new Node();
    reads_node1->set_distance_to_parent(r_dist);
    reads_node1->set_name(read->name);
    reads_node1->add_name_comment(read->comment);
    reads_node1->add_sequence( *read, mf->get_full_char_alphabet());

    Node * tmpnode = new Node();
    tmpnode->set_name("");

    tmpnode->add_left_child(node);
    tmpnode->add_right_child(reads_node1);

    tmpnode->align_sequences_this_node(mf,true);

    node->set_distance_to_parent(org_dist);

    int matching = 0;
    int aligned = 0;

    int node_start_pos1 = -1;
    int node_end_pos1 = -1;
    int node_start_pos2 = -1;
    int node_end_pos2 = -1;
    for( int k=1; k < tmpnode->get_sequence()->sites_length()-1; k++ )
    {
        Site *site = tmpnode->get_sequence()->get_site_at(k);

        if(site->get_children()->right_index == 1)
        {
            node_start_pos1 = k;
        }
        else if(site->get_children()->right_index == read->first_read_length)
        {
            node_end_pos1 = k;
        }
        else if(site->get_children()->right_index == read->first_read_length+1)
        {
            node_start_pos2 = k;
        }
        else if(site->get_children()->right_index > 0)
        {
            node_end_pos2 = k;
        }

        if(site->get_children()->right_index>=0 && site->get_children()->left_index>=0)
        {

            Site *site1 = tmpnode->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);
            Site *site2 = tmpnode->get_left_child()->get_sequence()->get_site_at(site->get_children()->left_index);

            if(site1->get_state() == site2->get_state())
                matching++;

            aligned++;
        }
    }

    double score = (double) matching/ (double) reads_node1->get_sequence()->sites_length();
//    cout<<"  "<<score<<" : "<<matching<<" "<<aligned<<" "<<node->get_sequence()->sites_length()<<" "<<reads_node1->get_sequence()->sites_length()<<endl;

    tmpnode->has_left_child(false);
    delete tmpnode;

    if(score>best_score)
    {
        read->node_start_pos1 = node_start_pos1;
        read->node_end_pos1 =   node_end_pos1;
        read->node_start_pos2 = node_start_pos2;
        read->node_end_pos2 =   node_end_pos2;
    }

    return score;
}

void Reads_alignment::remove_target_overlapping_identical_reads(vector<Fasta_entry> *reads, Model_factory *mf)
{
    cout<<"Removing identical reads mapped at overlapping positions.\n";

    vector<Fasta_entry>::iterator ri1 = reads->begin();

    for(;ri1 != reads->end();)
    {
        vector<Fasta_entry>::iterator ri2 = ri1;
        ri2++;
        for(;ri2 != reads->end();)
        {
            bool embedded = false;
            if( Settings_handle::st.is("pair-end")
                && (  (ri2->node_start_pos1 >= ri1->node_start_pos1 && ri2->node_end_pos1 <= ri1->node_end_pos1
                    && ri2->node_start_pos2 >= ri1->node_start_pos2 && ri2->node_end_pos2 <= ri1->node_end_pos2 )
                   || (ri2->node_start_pos1 <= ri1->node_start_pos1 && ri2->node_end_pos1 >= ri1->node_end_pos1
                    && ri2->node_start_pos2 <= ri1->node_start_pos2 && ri2->node_end_pos2 >= ri1->node_end_pos2 )
                   )
              )
            {
                embedded = true;
            }
            else if(!Settings_handle::st.is("pair-end")
                    && (  (ri2->node_start_pos1 >= ri1->node_start_pos1 && ri2->node_end_pos2 <= ri1->node_end_pos2 )
                       || (ri2->node_start_pos1 <= ri1->node_start_pos1 && ri2->node_end_pos2 >= ri1->node_end_pos2 )
                       )
                   )
            {
                embedded = true;
            }


            if( embedded )
            {
                Node * node = new Node();
                node->set_name("read pair");

                this->align_two_reads(node,&(*ri1),&(*ri2),mf);

                int matching = this->reads_pairwise_matching_sites(node);

                delete node;

                int r1_length = (int)ri1->sequence.length();
                int r2_length = (int)ri2->sequence.length();

                cout<<"idnetical? "<<matching<<" "<<r1_length<<" "<<r2_length<<endl;
                if(Settings_handle::st.is("pair-end")) // remove the midpoint marker
                {
                    r1_length--;
                    r2_length--;
                }

                if( matching == r2_length )
                {
                    cout<<"Read "<<ri2->name<<" is fully embedded in read "<<ri1->name<<" and overlapping sites are identical.  Read "<<ri2->name<<" is deleted.\n";
                    reads->erase(ri2);
                }
                else if( matching == r1_length )
                {
                    cout<<"Read "<<ri1->name<<" is fully embedded in read "<<ri2->name<<" and overlapping sites are identical.  Read "<<ri1->name<<" is deleted.\n";
                    reads->erase(ri1);
                    ri1--;
                    ri2 = reads->end();
                }
                else
                {
                    ri2++;
                }

            }
            else
            {
                ri2++;
            }

        }
        ri1++;
    }
}

void Reads_alignment::remove_target_overlapping_reads(vector<Fasta_entry> *reads)
{
    cout<<"Removing reads mapped at overlapping positions.\n";

    vector<Fasta_entry>::iterator ri1 = reads->begin();

    for(;ri1 != reads->end();)
    {
        vector<Fasta_entry>::iterator ri2 = ri1;
        ri2++;
        for(;ri2 != reads->end();)
        {
            if( Settings_handle::st.is("pair-end") )
            {
                if( ri2->node_start_pos1 >= ri1->node_start_pos1 && ri2->node_end_pos1 <= ri1->node_end_pos1 &&
                    ri2->node_start_pos2 >= ri1->node_start_pos2 && ri2->node_end_pos2 <= ri1->node_end_pos2 )
                {
                    cout<<"Read "<<ri2->name<<" is fully embedded in read "<<ri1->name<<".  Read "<<ri2->name<<" is deleted.\n";
                    reads->erase(ri2);
                }
                else
                {
                    ri2++;
                }

            }
            else
            {
                if( ri2->node_start_pos1 >= ri1->node_start_pos1 && ri2->node_end_pos2 <= ri1->node_end_pos2 )
                {
                    cout<<"Read "<<ri2->name<<" is fully embedded in read "<<ri1->name<<".  Read "<<ri2->name<<" is deleted.\n";
                    reads->erase(ri2);
                }
                else
                {
                    ri2++;
                }
            }

        }
        ri1++;
    }
}

void Reads_alignment::align_two_reads(Node *node, Fasta_entry *ri1, Fasta_entry *ri2, Model_factory *mf)
{
    double r_dist = Settings_handle::st.get("reads-distance").as<float>();

    node->set_name("read pair");

    Node * reads_node1 = new Node();
    reads_node1->set_distance_to_parent(r_dist);
    reads_node1->set_name(ri1->name);
    reads_node1->add_name_comment(ri1->comment);
    reads_node1->add_sequence( *ri1, mf->get_full_char_alphabet());
    node->add_right_child(reads_node1);

    Node * reads_node2 = new Node();
    reads_node2->set_distance_to_parent(r_dist);
    reads_node2->set_name(ri2->name);
    reads_node2->add_name_comment(ri2->comment);
    reads_node2->add_sequence( *ri2, mf->get_full_char_alphabet());
    node->add_left_child(reads_node2);

    node->align_sequences_this_node(mf,true);

}

int Reads_alignment::reads_pairwise_matching_sites(Node *node)
{
    int matching = 0;

    for( int k=1; k < node->get_sequence()->sites_length()-1; k++ )
    {
        Site *site = node->get_sequence()->get_site_at(k);
        if(site->get_children()->right_index>=0 && site->get_children()->left_index>=0)
        {
            Site *site1 = node->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);
            Site *site2 = node->get_left_child()->get_sequence()->get_site_at(site->get_children()->left_index);

            if(site1->get_state() == site2->get_state())
                matching++;
        }
    }

    return matching;
}

void Reads_alignment::remove_overlapping_reads(vector<Fasta_entry> *reads, Model_factory *mf)
{

    cout<<"Removing pairwise overlapping reads.\n";

    vector<Fasta_entry>::iterator ri1 = reads->begin();

    for(;ri1 != reads->end();)
    {
        vector<Fasta_entry>::iterator ri2 = ri1;
        ri2++;
        for(;ri2 != reads->end();)
        {

            Node * node = new Node();
            node->set_name("read pair");

            this->align_two_reads(node,&(*ri1),&(*ri2),mf);

            int matching = this->reads_pairwise_matching_sites(node);

            int r1_length = (int)ri1->sequence.length();
            int r2_length = (int)ri2->sequence.length();

            if(Settings_handle::st.is("pair-end")) // remove the midpoint marker
            {
                r1_length--;
                r2_length--;
            }

            if( matching == r1_length
                && matching == r2_length )
            {
                cout<<"Reads "<<ri1->name<<" and "<<ri2->name<<" are identical. Read "<<ri2->name<<" is deleted.\n";
                reads->erase(ri2);
            }
            else if(matching == r1_length)
            {
                cout<<"Read "<<ri1->name<<" is fully embedded in read "<<ri2->name<<". Read "<<ri1->name<<" is deleted.\n";
                reads->erase(ri1);
                ri1--;
                ri2 = reads->end();
            }
            else if(matching == r2_length)
            {
                cout<<"Read "<<ri2->name<<" is fully embedded in read "<<ri1->name<<".  Read "<<ri2->name<<" is deleted.\n";
                reads->erase(ri2);
            }
            else
            {
                ri2++;
            }

            delete node;
        }

        ri1++;
    }

}
