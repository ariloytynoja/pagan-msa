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

#include "main/reads_aligner.h"
#include "utils/exonerate_queries.h"
#include "utils/text_utils.h"
#include "utils/log_output.h"
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace ppa;

Reads_aligner::Reads_aligner(){}

void Reads_aligner::align(Node *root, Model_factory *mf, int count)
{

    // Handle also the first one correctly
    //
    if(!Settings_handle::st.is("ref-seqfile"))
    {
        root->get_sequence()->is_read_sequence(true);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    string file = Settings_handle::st.get("queryfile").as<string>();

    Log_output::write_header("Aligning reads ",0);

    Fasta_reader fr;
    vector<Fasta_entry> reads;
    Log_output::write_out("Reads data file: "+file+"\n",1);

    try
    {
        fr.read(file, reads, true);
        fr.remove_gaps(&reads);
    }
    catch (ppa::IOException& e) {
        Log_output::write_out("Error reading the reads file '"+file+"'.\nExiting.\n\n",0);
        exit(1);
    }

    int data_type = fr.check_sequence_data_type(&reads);

    if(!fr.check_alphabet(&reads,data_type))
        Log_output::write_out(" Warning: Illegal characters in input reads sequences removed!\n",2);


    // Merge overlapping and trim reads if selected, otherwise jsut update the sequence comment
    //
    this->merge_and_trim_and_sort( mf, &fr,  &reads );



    //////////////////////////////////////////////////////////////////
    //                                                              //
    //     Different read placement options  (until now)            //
    //                                                              //
    //////////////////////////////////////////////////////////////////


    // Alignment of translated RNA-seq clusters; can be reverse (currently just pileup)
    if( Settings_handle::st.is("pileup-alignment") && Settings_handle::st.is("find-best-orf"))
    {
        Log_output::write_header("Aligning reads: placement with ORF search",0);
        this->loop_translated_pileup_alignment(root,&reads,mf,count);
    }

    // Pileup: no search for optimal node or TID tags in the tree; can compare reverse strand
    else if( Settings_handle::st.is("pileup-alignment") || Settings_handle::st.is("align-reads-at-root") )
    {
        Log_output::write_header("Aligning reads: simple placement",0);
        this->loop_pileup_alignment(root,&reads,mf,count);
    }

    // Query placement: search for optimal node or TID tags in the tree; can compare reverse strand
    else if( Settings_handle::st.is("new-placement") )
    {
        Log_output::write_header("Aligning reads: new placement",0);
        this->loop_query_placement(root,&reads,mf,count);
    }

//    // NEW!
//    else if(Settings_handle::st.is("tagged-search"))
//    {
//        Log_output::write_header("Aligning reads: tagged search",0);
//        this->loop_tagged_placement(root,&reads,mf,count);
//    }

//    // NEW!
//    else if(Settings_handle::st.is("upwards-search"))
//    {
//        Log_output::write_header("Aligning reads: upwards search",0);
//        this->do_upwards_search(root,&reads,mf,count);
//    }

    // Default placement
    else
    {
        Log_output::write_header("Aligning reads: default placement",0);
        this->loop_default_placement(root,&reads,mf,count);
    }
}

void Reads_aligner::loop_pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    string ref_root_name = root->get_name();
    global_root = root;

    int start_i = 0;
    if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
        start_i = 1;

    int max_attempts = Settings_handle::st.get("query-cluster-attempts").as<int>();
    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    bool compare_reverse = Settings_handle::st.is("compare-reverse") && mf->get_sequence_data_type()==Model_factory::dna;

    for(int j=0; j < max_attempts; j++)
    {

        for(int i=start_i;i<(int)reads->size();i++)
        {

            if(reads->at(i).cluster_attempts >= max_attempts)
                continue;

            stringstream ss;
            ss<<"#"<<count<<"#";

            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            this->create_temp_node(node,ss.str(), global_root, &reads->at(i),false);

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: "+reads->at(i).name+" "+reads->at(i).comment+".",0);

            node->align_sequences_this_node(mf,true,false);
            this->compute_read_overlap(node,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap,&read_identity);


            Node * node_rc;
            float read_overlap_rc = -1;
            float read_identity_rc = -1;

            if(compare_reverse)
            {
                node_rc = new Node();
                this->create_temp_node(node_rc,ss.str(), global_root, &reads->at(i),true);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read (rc): "+reads->at(i).name+" "+reads->at(i).comment+".",0);

                node_rc->align_sequences_this_node(mf,true,false);
                this->compute_read_overlap(node_rc,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap_rc,&read_identity_rc);

                Log_output::write_out("forward overlap: "+Log_output::ftos(read_overlap)+"; backward overlap: "+Log_output::ftos(read_overlap_rc)+"\n",1);
            }
            else
            {
                Log_output::write_out("forward overlap: "+Log_output::ftos(read_overlap)+"\n",1);
            }

            reads->at(i).cluster_attempts++;


            if(read_overlap > read_overlap_rc && read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                global_root = node;

                if(compare_reverse)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }

                reads->at(i).cluster_attempts = max_attempts;
            }

            else if( read_overlap_rc > min_overlap && read_identity_rc > min_identity )
            {
                count++;
                global_root = node_rc;

                node->has_left_child(false);
                delete node;

                reads->at(i).cluster_attempts = max_attempts;
            }

            else
            {
                reads->at(i).cluster_attempts++;

                node->has_left_child(false);
                delete node;

                if(compare_reverse)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }
            }
        }
    }
}

void Reads_aligner::loop_translated_pileup_alignment(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    this->define_translation_tables();

    string ref_root_name = root->get_name();
    global_root = root;

    int start_i = 0;
    if(Settings_handle::st.is("queryfile") && !Settings_handle::st.is("ref-seqfile"))
        start_i = 1;

    int max_attempts = Settings_handle::st.get("query-cluster-attempts").as<int>();
    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    for(int j=0; j < max_attempts; j++)
    {

        for(int i=start_i;i<(int)reads->size();i++)
        {

            if(reads->at(i).cluster_attempts >= max_attempts)
                continue;

            Node * best_node = 0;
            Orf *best_orf = 0;

            vector<Orf> open_frames;
            this->find_orfs(&reads->at(i),&open_frames);

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: '"+reads->at(i).name+"'",0);

            if(open_frames.size() > 0)
            {
                float best_overlap = -1;
                float best_identity = -1;

                for(int h=0;h<(int)open_frames.size();h++)
                {
                    Node * node = new Node();
                    float read_overlap = -1;
                    float read_identity = -1;

                    this->create_temp_orf_node(node,global_root, &reads->at(i), &open_frames.at(h));

                    node->align_sequences_this_node(mf,true,false);
                    this->compute_read_overlap(node,reads->at(i).name,ref_root_name,global_root->get_name(),&read_overlap,&read_identity);

                    if(read_overlap > best_overlap || (read_overlap == best_overlap && read_overlap > read_identity) )
                    {
                        best_node = node;
                        best_orf = &open_frames.at(h);
                        best_overlap = read_overlap;
                        best_identity = read_identity;
                    }
                    else
                    {
                        node->has_left_child(false);
                        delete node;
                    }
                }

                if(best_overlap > min_overlap && best_identity > min_identity)
                {
                    reads->at(i).cluster_attempts = max_attempts;
                    stringstream cs;
                    cs<<"["<<best_orf->frame<<" "<<best_orf->start+1<<"-"<<best_orf->end+1<<"]";
                    best_node->get_right_child()->add_name_comment(best_node->get_right_child()->get_name_comment()+cs.str());

                    best_node->set_nhx_tid(best_node->get_left_child()->get_nhx_tid());
                    best_node->get_right_child()->set_nhx_tid(best_node->get_left_child()->get_nhx_tid());

                    best_node->get_right_child()->set_Orf(best_orf);

                    stringstream ss;
                    ss<<"#"<<count<<"#";
                    best_node->set_name(ss.str());

                    count++;
                    global_root = best_node;

                    ss.str(string());
                    ss.precision(2);
                    ss<<"orf "<<best_orf->frame<<" ("<<best_orf->start<<"-"<<best_orf->end<<"): overlap: "<<best_overlap<<", identity: "<<best_identity<<endl;
                    Log_output::write_out(ss.str(),2);
                }
                else
                {
                    reads->at(i).cluster_attempts++;

                    best_node->has_left_child(false);
                    delete best_node;
                }
            }
        }
    }
}

void Reads_aligner::find_orfs(Fasta_entry *read,vector<Orf> *open_frames)
{

    int min_orf_length = int (Settings_handle::st.get("min-orf-coverage").as<float>() * read->sequence.length() / 3 );

    string dna = read->dna_sequence;
    int length = dna.length()-1;

    for(int i=0;i<3;i++)
    {
        string prot;
        string sequence = dna.substr(i);
        int start_site = i;
        int end_site = i+2;

        for (unsigned int j=0; j<sequence.length(); j+=3)
        {
            string codon = sequence.substr(j,3);
            if (codon_to_aa.find(codon) == codon_to_aa.end()
                    || codon_to_aa.find(codon)->second == "X" || codon_to_aa.find(codon)->second == "-" )
            {
                if((int)prot.length() >= min_orf_length)
                {
                    Orf o;
                    o.translation = prot;
                    o.frame = i+1;
                    o.start = start_site;
                    o.end = end_site;
                    o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

                    open_frames->push_back(o);                    
                }

                prot = "";
                start_site = j+i+3;
            }
            else
            {
                prot += codon_to_aa.find(codon)->second;
            }
            end_site = j+i+2;
        }


        if((int)prot.length() >= min_orf_length)
        {
            Orf o;
            o.translation = prot;
            o.frame = i+1;
            o.start = start_site;
            o.end = end_site;
            o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

            open_frames->push_back(o);
        }
    }

    dna = this->reverse_complement(read->dna_sequence);
    for(int i=0;i<3;i++)
    {
        string prot;
        string sequence = dna.substr(i);
        int start_site = i;
        int end_site = i+2;

        for (unsigned int j=0; j<sequence.length(); j+=3)
        {
            string codon = sequence.substr(j,3);
            if (codon_to_aa.find(codon) == codon_to_aa.end()
                    || codon_to_aa.find(codon)->second == "X" || codon_to_aa.find(codon)->second == "-" )
            {
                if((int)prot.length() >= min_orf_length)
                {
                    Orf o;
                    o.translation = prot;
                    o.frame = -1*(i+1);
                    o.start = length - end_site;
                    o.end = length - start_site;
                    o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

                    open_frames->push_back(o);
                }

                prot = "";
                start_site = j+i+3;
            }
            else
            {
                prot += codon_to_aa.find(codon)->second;
            }
            end_site = j+i+2;
        }

        if((int)prot.length() >= min_orf_length)
        {
            Orf o;
            o.translation = prot;
            o.frame = -1*(i+1);
            o.start = length - end_site;
            o.end = length - start_site;
            o.dna_sequence = dna.substr(start_site,end_site-start_site+1);

            open_frames->push_back(o);
        }
    }
}

string Reads_aligner::reverse_complement(string dna)
{
    string rev = dna;
    reverse(rev.begin(), rev.end());

    Text_utils tu;

    tu.replace_all(rev,'A','Z');
    tu.replace_all(rev,'T','A');
    tu.replace_all(rev,'Z','T');
    tu.replace_all(rev,'C','Z');
    tu.replace_all(rev,'G','C');
    tu.replace_all(rev,'Z','G');
    tu.replace_all(rev,'R','Z');
    tu.replace_all(rev,'Y','R');
    tu.replace_all(rev,'Z','Y');
    tu.replace_all(rev,'K','Z');
    tu.replace_all(rev,'M','K');
    tu.replace_all(rev,'Z','M');
    tu.replace_all(rev,'B','Z');
    tu.replace_all(rev,'V','B');
    tu.replace_all(rev,'Z','V');
    tu.replace_all(rev,'D','Z');
    tu.replace_all(rev,'H','D');
    tu.replace_all(rev,'Z','H');

    return rev;
}

void Reads_aligner::define_translation_tables()
{
    string codon[66] = {"TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                        "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                        "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                        "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                        "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
                        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
                        "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
                        "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
                        "NNN", "---"
                       };
    string unaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "X", "W", "R", "R", "R", "R", "S", "S", "R", "R", "G", "G", "G", "G",
                       "X", "-"
                      };
    string mtaa[66] = {"F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "M", "M", "V", "V", "V", "V",
                       "S", "S", "S", "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A",
                       "Y", "Y", "X", "X", "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                       "C", "C", "W", "W", "R", "R", "R", "R", "S", "S", "X", "X", "G", "G", "G", "G",
                       "X", "-"
                      };

    if(Settings_handle::st.is("mt-translate"))
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],mtaa[i]));
            aa_to_codon.insert(make_pair(mtaa[i],codon[i]));
        }
    }
    else
    {
        for (int i=0; i<66; i++)
        {
            codon_to_aa.insert(make_pair(codon[i],unaa[i]));
            aa_to_codon.insert(make_pair(unaa[i],codon[i]));
        }
    }
}

void Reads_aligner::loop_tagged_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

}

void Reads_aligner::do_upwards_search(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    double r_dist = Settings_handle::st.get("query-distance").as<float>();
    Evol_model model = mf->alignment_model(r_dist+0.001,false);


    for(int i=0;i<(int)reads->size();i++)
//    for(int i=0;i<1;i++)
    {
        Node * current_root = root;
        string previous_hit = "";

        while(true)
        {
            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") aligning read: '"+reads->at(i).name+"'",0);

            Node * temp_node = new Node();
            temp_node->set_name("temp");

            Node * left_node = current_root;
            float org_distance = left_node->get_distance_to_parent();
            left_node->set_distance_to_parent(0.001);
            temp_node->add_left_child(left_node);

            Node * right_node = new Node();
            this->copy_node_details(right_node,&reads->at(i));
            temp_node->add_right_child(right_node);

            temp_node->align_sequences_this_node(mf,true,false);


            vector<int> all_scores;
            for(int j=0;j<current_root->get_number_of_nodes();j++)
            {
                all_scores.push_back(0);
            }

            vector<int> one_score;
            int read_length = 0;

            //

            vector<float> all_dists;
            for(int j=0;j<current_root->get_number_of_nodes();j++)
            {
                all_dists.push_back(0);
            }

            vector<float> one_dist;
            float read_max_dist = 0;


            for(int j=1;j<temp_node->get_sequence()->sites_length()-1;j++)
            {
                Site_children *offspring = temp_node->get_sequence()->get_site_at(j)->get_children();

                int lj = offspring->left_index;
                int rj = offspring->right_index;
                if(rj>=0)
                {
                    int qs = temp_node->get_right_child()->get_sequence()->get_site_at(rj)->get_state();

                    if(lj>=0)
                    {
                        one_score.clear();
                        current_root->get_state_identity(lj,qs,&one_score,true);

                        for(int l=0;l<(int)all_scores.size();l++)
                            all_scores.at(l)+=one_score.at(l);

                        //

                        one_dist.clear();
                        current_root->get_subst_distance(lj,qs,&one_dist,&model,true);

                        for(int l=0;l<(int)all_scores.size();l++)
                            all_dists.at(l)+=one_dist.at(l);
                    }
                    read_length++;

                    read_max_dist += model.score(qs,qs);
                }
            }


            double max_score = -1;
            string max_node = current_root->get_name();
            int max_node_index = 0;

            double max_dist = -1;
            string max_dist_node = current_root->get_name();
            int max_dist_node_index = 0;

            vector<Node*> score_nodes;
            current_root->get_all_nodes(&score_nodes,true);
            for(int l=0;l<(int)all_scores.size();l++)
            {
                float score = (float)all_scores.at(l)/(float)read_length;
                if(score>=max_score)
                {
                    max_score = score;
                    max_node = score_nodes.at(l)->get_name();
                    max_node_index = l;
                }
                float dist = (float)all_dists.at(l)/(float)read_max_dist;
                if(dist>=max_dist)
                {
                    max_dist = dist;
                    max_dist_node = score_nodes.at(l)->get_name();
                    max_dist_node_index = l;
                }
            }

//            cout<<endl;
//            cout<<"s1 "<<reads->at(i).name<<": "<<max_node<<" ("<<max_score<<")\n";
//            cout<<"s2 "<<reads->at(i).name<<": "<<max_dist_node<<" ("<<max_dist<<")\n";


//            if(previous_hit=="")
////            if(true)
//            {
//                vector<Fasta_entry> aligned_seqs;
//                temp_node->get_alignment(&aligned_seqs,true);
//                ofstream out(reads->at(i).name.c_str());
//                for(int j=0;j<aligned_seqs.size();j++)
//                    out<<">"<<aligned_seqs.at(j).name<<"\n"<<aligned_seqs.at(j).sequence<<endl;
//            }

            left_node->set_distance_to_parent(org_distance);
            temp_node->has_left_child(false);
            delete temp_node;


            if(max_dist_node == previous_hit)
            {
                reads->at(i).node_score = max_dist;
                reads->at(i).node_to_align = max_dist_node;

                break;
            }
            else
            {
                previous_hit = max_dist_node;
                Node *temp = current_root->get_parent_node(max_dist_node);
                if(temp != 0)
                {
                    current_root = temp;
//                    cout<<"new root "<<current_root->get_name()<<endl;
                }
                else
                {
                    reads->at(i).node_score = max_dist;
                    reads->at(i).node_to_align = max_dist_node;
                    break;
                }
            }

        }


//        cout<<"placement "<<reads->at(i).name<<": "<<reads->at(i).node_to_align<<" ("<<reads->at(i).node_score <<")\n";
    }

}

void Reads_aligner::loop_default_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;

    // vector of node names giving the best node for each read
    //

    if(Settings_handle::st.is("very-fast-placement"))
    {
        Log_output::write_header("Placing query sequences: very fast Exonerate placement",0);
        this->find_nodes_for_all_reads_together(root, reads, mf);
    }
    else if(Settings_handle::st.is("fast-placement"))
    {
        Log_output::write_header("Placing query sequences: fast Exonerate placement",0);
        this->find_nodes_for_all_reads(root, reads, mf);
    }
    // NEW!
    else if(Settings_handle::st.is("tagged-search"))
    {
        Log_output::write_header("Placing query sequences: tagged search",0);
//        this->loop_tagged_placement(root,&reads,mf,count);
    }
    // NEW!
    else if(Settings_handle::st.is("upwards-search"))
    {
        Log_output::write_header("Placing query sequences: upwards search",0);
        this->do_upwards_search(root,reads,mf);
    }
    else
        this->find_nodes_for_reads(root, reads, mf);

    if(Settings_handle::st.is("placement-only"))
        exit(0);

    Log_output::write_header("Aligning query sequences",0);


    set<string> unique_nodeset;
    for(int i=0;i<(int)reads->size();i++)
    {
        stringstream nodestream;
        nodestream << reads->at(i).node_to_align;
        string val;
        while(nodestream >> val)
        {
            unique_nodeset.insert(val);
        }
    }

    if(unique_nodeset.find("discarded_read") != unique_nodeset.end())
        unique_nodeset.erase(unique_nodeset.find("discarded_read"));

    vector<string> unique_nodes;
    for(set<string>::iterator sit = unique_nodeset.begin(); sit != unique_nodeset.end(); sit++)
    {
        unique_nodes.push_back(*sit);
    }
    sort(unique_nodes.begin(),unique_nodes.end(),Reads_aligner::nodeIsSmaller);

    map<string,Node*> nodes_map;
    root->get_all_nodes(&nodes_map);

    // do one tagged node at time
    //
    for(vector<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
    {
        vector<Fasta_entry> reads_for_this;

        for(int i=0;i<(int)reads->size();i++)
        {
            stringstream nodestream;
            nodestream << reads->at(i).node_to_align;
            string val;
            while(nodestream >> val)
            {
                if(val == *sit)
                    reads_for_this.push_back(reads->at(i));
            }
        }

        this->sort_reads_vector(&reads_for_this);

        // remove fully overlapping reads that are mapped to this node
        //
        if(Settings_handle::st.is("rank-reads-for-nodes"))
        {
            bool print_log = false;

            if(Settings_handle::st.is("discard-overlapping-identical-reads"))
            {
                this->remove_target_overlapping_identical_reads(&reads_for_this,mf);
                print_log = true;
            }
            else if(Settings_handle::st.is("discard-overlapping-reads"))
            {
                this->remove_target_overlapping_reads(&reads_for_this);
                print_log = true;
            }
            else if(Settings_handle::st.is("discard-pairwise-overlapping-reads"))
            {
                this->remove_overlapping_reads(&reads_for_this,mf);
                print_log = true;
            }

            if(print_log && Settings::noise>1)
            {
                stringstream ss;
                ss<<"After removing overlapping ones, for node "<<*sit<<" reads remaining:\n";
                for(int i=0;i<(int)reads_for_this.size();i++)
                    ss<<" "<<reads_for_this.at(i).name<<" "<<reads_for_this.at(i).node_score<<endl;
                ss<<endl;
                Log_output::write_out(ss.str(),2);
            }

        }
        else
        {
            if( Settings_handle::st.is("discard-overlapping-identical-reads") ||
                     Settings_handle::st.is("discard-overlapping-reads") )
            {
                Log_output::write_out("\nWarning: without ranking the reads for nodes, one cannot resolve overlap between reads. The flag has no effect!\n\n",2);
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
            Node * node = new Node();

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            current_root->set_distance_to_parent(0.001);
            node->add_left_child(current_root);

            Node * reads_node = new Node();
            this->copy_node_details(reads_node,&reads_for_this.at(i));


            node->add_right_child(reads_node);

            node->set_nhx_tid(node->get_left_child()->get_nhx_tid());
            node->get_right_child()->set_nhx_tid(node->get_left_child()->get_nhx_tid());

            Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read: '"+reads_for_this.at(i).name+"'",0);

            node->align_sequences_this_node(mf,true,false);


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
            if(single_ref_sequence)
            {
                global_root = current_root;
            }
            else
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

void Reads_aligner::loop_query_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count)
{
    bool single_ref_sequence = false;
    if(root->get_number_of_leaves()==1)
        single_ref_sequence = true;

    global_root = root;
    bool compare_reverse = Settings_handle::st.is("compare-reverse") && mf->get_sequence_data_type()==Model_factory::dna;

    this->find_nodes_for_queries(root, reads, mf);

    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    if(min_overlap<0)
        min_overlap = 0;
    if(min_identity<0)
        min_identity = 0;

    Log_output::write_header("Aligning query sequences",0);


    set<string> unique_nodeset;
    for(int i=0;i<(int)reads->size();i++)
    {
        stringstream nodestream;
        nodestream << reads->at(i).node_to_align;
        string val;
        while(nodestream >> val)
        {
            unique_nodeset.insert(val);
        }
    }

    if(unique_nodeset.find("discarded_read") != unique_nodeset.end())
        unique_nodeset.erase(unique_nodeset.find("discarded_read"));

    vector<string> unique_nodes;
    for(set<string>::iterator sit = unique_nodeset.begin(); sit != unique_nodeset.end(); sit++)
    {
        unique_nodes.push_back(*sit);
    }
    sort(unique_nodes.begin(),unique_nodes.end(),Reads_aligner::nodeIsSmaller);

    map<string,Node*> nodes_map;
    root->get_all_nodes(&nodes_map);

    // do one tagged node at time
    //
    for(vector<string>::iterator sit = unique_nodes.begin(); sit != unique_nodes.end(); sit++)
    {
        vector<Fasta_entry> reads_for_this;

        for(int i=0;i<(int)reads->size();i++)
        {
            stringstream nodestream;
            nodestream << reads->at(i).node_to_align;
            string val;
            while(nodestream >> val)
            {
                if(val == *sit)
                    reads_for_this.push_back(reads->at(i));
            }
        }

        this->sort_reads_vector(&reads_for_this);


        string ref_node_name = *sit;

        Node *current_root = nodes_map.find(ref_node_name)->second;
        double orig_dist = current_root->get_distance_to_parent();

        bool alignment_done = false;
        int alignments_done = 0;

        // align the remaining reads to this node
        //
        for(int i=0;i<(int)reads_for_this.size();i++)
        {
            Node * node = new Node();
            float read_overlap = -1;
            float read_identity = -1;

            stringstream ss;
            ss<<"#"<<count<<"#";
            node->set_name(ss.str());

            if(reads_for_this.at(i).query_strand != Fasta_entry::reverse_strand)
            {
                this->create_temp_node(node,ss.str(), current_root, &reads_for_this.at(i),false);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read: '"+reads_for_this.at(i).name+"'",0);

                node->align_sequences_this_node(mf,true,false);
                this->compute_read_overlap(node,reads_for_this.at(i).name,ref_node_name,current_root->get_name(),&read_overlap,&read_identity);
            }

            Node * node_rc;
            float read_overlap_rc = -1;
            float read_identity_rc = -1;

            bool reverse_computed = false;

            if(compare_reverse && reads_for_this.at(i).query_strand != Fasta_entry::forward_strand)
            {
                node_rc = new Node();
                this->create_temp_node(node_rc,ss.str(), current_root, &reads_for_this.at(i),true);

                Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads_for_this.size())+") aligning read (rc): "+reads_for_this.at(i).name+".",0);

                node_rc->align_sequences_this_node(mf,true,false);
                this->compute_read_overlap(node_rc,reads_for_this.at(i).name,ref_node_name,current_root->get_name(),&read_overlap_rc,&read_identity_rc);

                reverse_computed = true;
            }

            Log_output::write_out("forward overlap/identity: "+Log_output::ftos(read_overlap)+"/"+Log_output::ftos(read_identity)+"; backward: "+Log_output::ftos(read_overlap_rc)+"/"+Log_output::ftos(read_identity_rc)+"\n",1);


            if(read_overlap > read_overlap_rc && read_overlap > min_overlap && read_identity > min_identity)
            {
                count++;
                current_root = node;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }
            }

            else if( read_overlap_rc > min_overlap && read_identity_rc > min_identity )
            {
                count++;
                current_root = node_rc;

                if( orig_dist > current_root->get_distance_to_parent() )
                    orig_dist -= current_root->get_distance_to_parent();

                alignment_done = true;
                alignments_done++;

                node->has_left_child(false);
                delete node;
            }

            else
            {
                node->has_left_child(false);
                delete node;

                if(reverse_computed)
                {
                    node_rc->has_left_child(false);
                    delete node_rc;
                }
            }
        }

        current_root->set_distance_to_parent(orig_dist);

        if(alignment_done)
        {
            if(single_ref_sequence)
            {
                global_root = current_root;
            }
            else
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


void Reads_aligner::merge_reads_only()
{
    string file = Settings_handle::st.get("queryfile").as<string>();

    Fasta_reader fr;
    vector<Fasta_entry> reads;
    Log_output::write_out("Reads data file: "+file+"\n",1);

    fr.read(file, reads, true);

    int data_type = fr.check_sequence_data_type(&reads);
    Model_factory mf(data_type);

    if(!fr.check_alphabet(&reads,Model_factory::dna))
        Log_output::write_out("\n Warning: Illegal characters in input sequences removed!\n",2);

    float *dna_pi = fr.base_frequencies();

    mf.dna_model(dna_pi,&Settings_handle::st);

    this->merge_paired_reads( &reads, &mf );

    string path = "outfile";
    if(Settings_handle::st.is("overlap-merge-file"))
    {
        path = Settings_handle::st.get("overlap-merge-file").as<string>();
    }

    Log_output::write_out("Reads output file: "+path+".fastq\n",1);

    fr.write_fastq(path,reads);
}

void Reads_aligner::merge_and_trim_and_sort(Model_factory *mf, Fasta_reader *fr, vector<Fasta_entry> *reads)
{

    // Merge overlapping reads
    //
    if( Settings_handle::st.is("overlap-pair-end") )
    {
        if(Settings_handle::st.is("trim-before-merge"))
        {
            fr->trim_fastq_reads(reads);
         }

        this->merge_paired_reads( reads, mf );

        if(Settings_handle::st.is("overlap-merge-file"))
        {

            string path = Settings_handle::st.get("overlap-merge-file").as<string>();

            Log_output::write_out("Reads output file: "+path+".fastq\n",1);

            fr->write_fastq(path,*reads);
        }
        if( Settings_handle::st.is("pair-end") )
        {
            Log_output::write_out(" Warning: both '--overlap-pair-end' and '--pair-end' options defined.\n"
                    " Pairing of overlapping reads may cause duplicated sequence regions.\n\n",2);
        }
    }

    // Trim reads ends
    //
    if(Settings_handle::st.is("trim-read-ends"))
    {
        Log_output::write_header("Aligning reads: trim read ends",0);
        fr->trim_fastq_reads( reads );
    }


    // Couple paired reads
    //
    if( Settings_handle::st.is("pair-end") )
    {
        Log_output::write_header("Aligning reads: find read pairs",0);
        this->find_paired_reads( reads );
    }

    // Add comment about trimming -- or just blank
    //
    vector<Fasta_entry>::iterator fit1 = reads->begin();

    for(;fit1 != reads->end();fit1++)
    {
        if(Settings_handle::st.is("trim-read-ends"))
        {
            stringstream trimming;
            trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end;
            fit1->comment += trimming.str();
        }
    }

    // Sort reads by copy number
    //
    if(Settings_handle::st.is("use-duplicate-weigths") && not Settings_handle::st.is("no-read-ordering"))
    {
        this->sort_reads_vector_by_duplicate_number( reads );
    }
}


void Reads_aligner::merge_paired_reads(vector<Fasta_entry> *reads, Model_factory *mf)
{

    vector<Fasta_entry>::iterator fit1 = reads->begin();

    Log_output::write_header("Aligning reads: finding read pairs",0);

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

                Node * node_left = new Node();
                this->copy_node_details(node_left,&(*fit1));

                Node * node_right = new Node();
                this->copy_node_details(node_right,&(*fit2));

                Node * node = new Node();
                node->set_name("merge");

                node->add_left_child(node_left);
                node->add_right_child(node_right);

                node->set_nhx_tid(node->get_left_child()->get_nhx_tid());

                Log_output::write_msg("aligning paired reads: "+fit1->name+" and "+fit2->name+"\n",0);

                node->align_sequences_this_node(mf,true,true);

                Sequence *anc_seq = node->get_sequence();
                int x = 0; int y = 0;
                stringstream l_seq("");
                stringstream r_seq("");

                int overlap = 0;
                int identical = 0;

                for(int i=1;i<anc_seq->sites_length()-1;i++)
                {
                    int path_state = anc_seq->get_site_at(i)->get_path_state();

                    if(path_state==Site::xgapped)
                    {
                        l_seq<<fit1->sequence.at(x); r_seq<<"-";
                        x++;
                    }
                    else if(path_state==Site::ygapped)
                    {
                        l_seq<<"-"; r_seq<<fit2->sequence.at(y);
                        y++;
                    }
                    else if(path_state==Site::matched)
                    {
                        l_seq<<fit1->sequence.at(x); r_seq<<fit2->sequence.at(y);

                        overlap++;
                        if(fit1->sequence.at(x) == fit2->sequence.at(y))
                            identical++;

                        x++;y++;
                    }
                    else
                    {
                        Log_output::write_out("Reads_aligner::merge_paired_reads: Error in pair-end merge alignment\n",2);
                    }
                }


                if(Settings::noise>2)
                {
                    stringstream ss;
                    ss<<"Alignment read pair "<<fit1->name<<" and "<<fit2->name<<".\n";
                    ss<<l_seq.str()<<endl<<r_seq.str()<<endl<<endl;
                    ss<<"overlap "<<overlap<<", identical "<<identical<<endl;
                    Log_output::write_out(ss.str(),3);
                }

                if( ( overlap >= Settings_handle::st.get("overlap-minimum").as<int>() &&
                      float(identical)/overlap >= Settings_handle::st.get("overlap-identity").as<float>() )
                 || ( overlap == identical && overlap >= Settings_handle::st.get("overlap-identical-minimum").as<int>() )
                    )
                {

                    Log_output::write_out("Merging "+fit1->name+" and "+fit2->name+": new name ",2);

                    fit1->name = name1.substr(0,name1.length()-1)+"m12";
                    fit1->comment = fit2->comment;

                    Log_output::write_out(fit1->name+".\n",2);

                    string seq = "";
                    string qsc = "";

                    x=0;y=0;
                    for(int i=1;i<anc_seq->sites_length()-1;i++)
                    {
                        int path_state = anc_seq->get_site_at(i)->get_path_state();

                        if(path_state==Site::xgapped)
                        {
                            seq += fit1->sequence.at(x);
                            qsc += fit1->quality.at(x);
                            x++;
                        }
                        else if(path_state==Site::ygapped)
                        {
                            seq += fit2->sequence.at(y);
                            qsc += fit2->quality.at(y);
                            y++;
                        }
                        else if(path_state==Site::matched)
                        {
                            if(static_cast<int>( fit1->quality.at(x) ) > static_cast<int>( fit2->quality.at(y) ))
                            {
                                seq += fit1->sequence.at(x);
                                qsc += fit1->quality.at(x);
                            }
                            else
                            {
                                seq += fit2->sequence.at(y);
                                qsc += fit2->quality.at(y);
                            }
                            x++;y++;
                        }
                    }

                    fit1->sequence = seq;
                    fit1->quality = qsc;

                    reads->erase(fit2);
                    delete node;

                    continue;
                }
                else
                {
                    fit2++;
                }


                delete node;
                continue;
            }

            fit2++;
        }
    }
}

void Reads_aligner::find_paired_reads(vector<Fasta_entry> *reads)
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
                Log_output::write_out("Pairing "+fit1->name+" and "+fit2->name+": new name ",2);

                fit1->name = name1.substr(0,name1.length()-1)+"p12";

                Log_output::write_out(fit1->name+".\n",2);

                fit1->comment = fit2->comment;
                fit1->first_read_length = fit1->sequence.length();
                fit1->sequence += "0"+fit2->sequence;
                if(fit1->quality!="")
                    fit1->quality += "0"+fit2->quality;

                if(Settings_handle::st.is("trim-read-ends"))
                {
                    stringstream trimming;
                    trimming << "P1ST"<<fit1->trim_start<<":P1ET"<<fit1->trim_end<<":P2ST"<<fit2->trim_start<<":P2ET"<<fit2->trim_end;
                    fit1->comment += trimming.str();
                }

                reads->erase(fit2);

                continue;
            }

            fit2++;
        }
    }


    if(Settings_handle::st.is("trim-read-ends"))
    {
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
}

void Reads_aligner::copy_node_details(Node *reads_node,Fasta_entry *read,bool turn_revcomp)
{
    double r_dist = Settings_handle::st.get("query-distance").as<float>();

    reads_node->set_distance_to_parent(r_dist);
    reads_node->set_name(read->name);
    reads_node->add_name_comment(read->comment);
    reads_node->add_sequence( *read, read->data_type, false, true, turn_revcomp);
    reads_node->get_sequence()->is_read_sequence(true);

}


bool Reads_aligner::read_alignment_overlaps(Node * node, string read_name, string ref_node_name)
{
    float min_overlap = Settings_handle::st.get("min-query-overlap").as<float>();
    float min_identity = Settings_handle::st.get("min-query-identity").as<float>();

    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;
    int matched = 0;

    if(Settings_handle::st.is("pileup-alignment") && !Settings_handle::st.is("overlap-with-reference"))
    {
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
    }
    else
    {
        for( int j=0; j < node_sequence->sites_length(); j++ )
        {
            bool read_has_site = node->has_site_at_alignment_column(j,read_name);
            bool ref_root_has_site = node->has_site_at_alignment_column(j,ref_node_name);

            if(read_has_site)
                read_length++;

            if(read_has_site && ref_root_has_site)
            {
                aligned++;

                int state_read = node->get_state_at_alignment_column(j,read_name);
                int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);
                if(state_read == state_ref)
                    matched++;
            }
        }
    }

    if( !Settings_handle::st.is("pileup-alignment") )
    {
        stringstream ss;
        ss<<"aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"];"<<
            " identical positions "<<(float)matched/(float)aligned<<" ["<<matched<<"/"<<aligned<<"]"<<endl;
        Log_output::write_out(ss.str(),2);
    }

    if(Settings_handle::st.is("pileup-alignment") && !Settings_handle::st.is("overlap-with-reference"))
    {
        if( (float)aligned/(float)read_length >= min_overlap )
            return true;
        else
            return false;
    }

    if( (float)aligned/(float)read_length >= min_overlap && (float)matched/(float)aligned >= min_identity)
    {
        return true;
    }
    else if( (float)aligned/(float)read_length < min_overlap && (float)matched/(float)aligned < min_identity )
    {
        stringstream ss;
        ss<<"Warning: read "<<read_name<<" dropped using the minimum overlap cut-off of "<<min_overlap<<
                " and the minimum identity cut-off of "<<min_identity<<"."<<endl;
        Log_output::write_out(ss.str(),2);

        return false;
    }
    else if( (float)aligned/(float)read_length < min_overlap)
    {

        stringstream ss;
        ss<<"Warning: read "<<read_name<<" dropped using the minimum overlap cut-off of "<<min_overlap<<"."<<endl;
        Log_output::write_out(ss.str(),2);

        return false;
    }
    else if( (float)matched/(float)aligned < min_identity)
    {

        stringstream ss;
        ss<<" Warning: read "<<read_name<<" dropped using the minimum identity cut-off of "<<min_identity<<"."<<endl;
        Log_output::write_out(ss.str(),2);

        return false;
    }

    return false;
}

float Reads_aligner::read_alignment_overlap(Node * node, string read_name, string ref_node_name)
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

    stringstream ss;
    ss<<"aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"]"<<endl;
    Log_output::write_out(ss.str(),2);

    return (float)aligned/(float)read_length;
}

void Reads_aligner::read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity)
{
    Sequence *node_sequence = node->get_sequence();

    int aligned = 0;
    int read_length = 0;
    int matched = 0;

    if(Settings_handle::st.is("overlap-with-reference"))
    {
        for( int j=0; j < node_sequence->sites_length(); j++ )
        {
            bool read_has_site = node->has_site_at_alignment_column(j,read_name);
            bool ref_root_has_site = node->has_site_at_alignment_column(j,ref_node_name);

            if(read_has_site)
                read_length++;

            if(read_has_site && ref_root_has_site)
            {
                aligned++;

                int state_read = node->get_state_at_alignment_column(j,read_name);
                int state_ref  = node->get_state_at_alignment_column(j,ref_node_name);
                if(state_read == state_ref)
                    matched++;
            }
        }
    }
    else
    {
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
     }

    stringstream ss;
    ss<<"aligned positions "<<(float)aligned/(float)read_length<<" ["<<aligned<<"/"<<read_length<<"]"<<endl;
    Log_output::write_out(ss.str(),2);
    ss.str(string());
    ss<<"matched positions "<<(float)matched/(float)aligned<<" ["<<matched<<"/"<<aligned<<"]"<<endl;
    Log_output::write_out(ss.str(),2);

    *overlap  = (float)aligned/(float)read_length;
    *identity = (float)matched/(float)aligned;

}

bool Reads_aligner::correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map)
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
        if(!mit->second->is_leaf() && mit->second->get_left_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_left_child(current_root);
            parent_found = true;
        }

        if(!mit->second->is_leaf() && mit->second->get_right_child()->get_name() == ref_node_name)
        {
            current_parent = mit->second;
            current_parent->add_right_child(current_root);
            is_left_child = false;
            parent_found = true;
        }
    }


    if(parent_found)
    {
        stringstream ss;
        ss<<"Parent of "<<ref_node_name<<" is "<<current_parent->get_name()<<"; "
          <<alignments_done<<" alignments done.";
        Log_output::write_out(ss.str(),3);

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
            Log_output::write_out("Site index needs correcting.\n",3);
        else
            Log_output::write_out("Site index not changed.\n",3);


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
        Log_output::write_out("No parent for "+ref_node_name+" found. Assuming that this is root.\n",2);

        return false;
    }

}

void Reads_aligner::get_target_node_names(Node *root,multimap<string,string> *tid_nodes, bool *ignore_tid_tags)
{

    if(Settings_handle::st.is("test-every-internal-node"))
    {
        root->get_internal_node_names(tid_nodes);
    }
    else if(Settings_handle::st.is("test-every-terminal-node"))
    {
        root->get_terminal_node_names(tid_nodes);
    }
    else if(Settings_handle::st.is("test-every-node"))
    {
        root->get_node_names(tid_nodes);
    }
    else
    {
        root->get_node_names_with_tid_tag(tid_nodes);
        if((int)tid_nodes->size()>0)
        {
            *ignore_tid_tags = false;
        }
        else
        {
            Log_output::write_out("No tagged nodes found. Considering all nodes!\n\n",0);
            tid_nodes->clear();
            root->get_node_names(tid_nodes);
        }
    }
}

void Reads_aligner::find_nodes_for_queries(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{
    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }


    // Get the original target nodes
    //
    multimap<string,string> tid_nodes;
    bool ignore_tid_tags = true;
    bool compare_reverse = Settings_handle::st.is("compare-reverse") && mf->get_sequence_data_type()==Model_factory::dna;

    this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

    // Handle one read at time
    //
    for(int i=0;i<(int)reads->size();i++)
    {
        reads->at(i).node_score = -1.0;
        reads->at(i).query_strand = Fasta_entry::unknown_strand;

        string tid = reads->at(i).tid;
        if( ignore_tid_tags )
            tid = "<empty>";

        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping read: '"+reads->at(i).name+" "+reads->at(i).comment+"'",0);

        // Call Exonerate to reduce the search space
        //
        map<string,hit> exonerate_hits;

        if(has_exonerate && Settings_handle::st.is("use-exonerate-local") )
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits, true,ignore_tid_tags);

            if(tid_nodes.size()>1 && Settings_handle::st.is("use-exonerate-gapped"))
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,ignore_tid_tags);
        }
        else if(has_exonerate && Settings_handle::st.is("use-exonerate-gapped") )
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,ignore_tid_tags);
        }


        // Handle (or reject) reads discarded by Exonerate
        //
        if(reads->at(i).node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("exhaustive-placement"))
            {
                reads->at(i).node_to_align == "";
                tid_nodes.clear();

                this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
            }
            else
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
                Log_output::write_out(ss.str(),2);

                continue;
            }
        }


        // Has TID or exhaustive search: now find the one with the best match
        //
        if(tid != "")
        {
            int matches = tid_nodes.count(tid);

            if( ignore_tid_tags )
                matches = tid_nodes.size();


            // Has TID but no matching node
            //
            if(matches == 0)
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_to_align = root->get_name();
            }
            // All done; continue


            // Has only one matching node and no need for ranking
            //
            else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
            {
                multimap<string,string>::iterator tit = tid_nodes.find(tid);

                if( ignore_tid_tags )
                    tit = tid_nodes.begin();

                if(tit != tid_nodes.end())
                {
                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_to_align = tit->second;

                    multimap<string,hit>::iterator ith = exonerate_hits.find(tit->second);
                    if(ith != exonerate_hits.end())
                    {
                        hit h = ith->second;

                        reads->at(i).local_qstart = h.q_start;
                        reads->at(i).local_qend = h.q_end;
                        reads->at(i).local_tstart = h.t_start;
                        reads->at(i).local_tend = h.t_end;
                    }
                }
            }
            // All done; continue


            // Has TID and matching nodes, or exhaustive search
            //
            else
            {
                double best_score = -HUGE_VAL;
                string best_node = root->get_name();
                int query_strand = Fasta_entry::forward_strand;

                map<string,Node*> nodes;
                root->get_all_nodes(&nodes);

                multimap<string,string>::iterator tit;
                int matching_nodes = 0;

                if( ignore_tid_tags )
                {
                    tit = tid_nodes.begin();
                    matching_nodes = tid_nodes.size();
                }
                else
                {
                    tit = tid_nodes.find(tid);
                    matching_nodes = tid_nodes.count(tid);
                }

                if(tit != tid_nodes.end())
                {

                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                    Log_output::write_out(ss.str(),2);

                    while(tit != tid_nodes.end() && matching_nodes>0)
                    {
                        map<string,Node*>::iterator nit = nodes.find(tit->second);
                        double score = this->read_match_score( nit->second, &reads->at(i), mf, best_score);

                        stringstream ss;
                        ss<<tit->second<<" with score "<<score<<" (simple p-distance)\n";
                        Log_output::write_out(ss.str(),2);

                        if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                        {
                            best_score = score;
                            best_node.append(" "+tit->second);
                        }
                        else if(score>best_score)
                        {
                            best_score = score;
                            best_node = tit->second;
                        }

                        if(compare_reverse)
                        {
                            Fasta_entry rev_seq = reads->at(i);
                            rev_seq.sequence = this->reverse_complement(rev_seq.sequence);

                            double score = this->read_match_score( nit->second, &rev_seq, mf, best_score);

                            stringstream ss;
                            ss<<tit->second<<"(rc) with score "<<score<<" (simple p-distance)\n";
                            Log_output::write_out(ss.str(),2);

                            if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                            {
                                best_score = score;
                                best_node.append(" "+tit->second);
                                query_strand = Fasta_entry::reverse_strand;
                            }
                            else if(score>best_score)
                            {
                                best_score = score;
                                best_node = tit->second;
                                query_strand = Fasta_entry::reverse_strand;
                            }
                        }

                        tit++;
                        matching_nodes--;
                    }
                }

                if(best_score<0.05)
                {
                    if(Settings_handle::st.is("align-bad-reads-at-root"))
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                        reads->at(i).node_to_align = root->get_name();

                    }
                    else
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                        reads->at(i).node_to_align = "discarded_read";
                    }
                }
                else
                {
                    stringstream ss;
                    ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_score = best_score;
                    reads->at(i).node_to_align = best_node;
                    reads->at(i).query_strand = query_strand;
                }
            }
            // All done; continue

        }


        // no TID, aligning at root

        else
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") has no tid. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = root->get_name();

        }
    }
}

void Reads_aligner::find_nodes_for_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{

    multimap<string,string> tid_nodes;
    bool ignore_tid_tags = true;

    this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

    ofstream pl_output;

    if(Settings_handle::st.is("placement-file"))
    {
        string fname = Settings_handle::st.get("placement-file").as<string>();
        pl_output.open(fname.append(".tsv").c_str());
    }

    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }

    for(int i=0;i<(int)reads->size();i++)
    {
        reads->at(i).node_score = -1.0;

        string tid = reads->at(i).tid;
        if( ignore_tid_tags )
            tid = "<empty>";

        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping read: '"+reads->at(i).name+"'",0);

        // Call Exonerate to reduce the search space

        map<string,hit> exonerate_hits;

        if(has_exonerate && Settings_handle::st.is("use-exonerate-local") )
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits, true,ignore_tid_tags);

            if(tid_nodes.size()>1 && Settings_handle::st.is("use-exonerate-gapped"))
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,ignore_tid_tags);
        }

        else if(has_exonerate && Settings_handle::st.is("use-exonerate-gapped"))
        {
            tid_nodes.clear();
            this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);

            if(tid_nodes.size()>1)
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_hits,false,ignore_tid_tags);
        }

        // Discarded by Exonerate

        if(reads->at(i).node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("exhaustive-placement"))
            {
                reads->at(i).node_to_align == "";
                tid_nodes.clear();

                this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
            }
            else
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
                Log_output::write_out(ss.str(),2);

                continue;
            }
        }

        // Has TID or exhaustive search
        if(tid != "")
        {
            int matches = tid_nodes.count(tid);

            if( ignore_tid_tags )
                matches = tid_nodes.size();


            // has TID but no matching node

            if(matches == 0)
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_to_align = root->get_name();

                if(Settings_handle::st.is("placement-file"))
                {
                    pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (no match)"<<endl;
                }

            }


            // has only one matching node and no need for ranking

            else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
            {
                multimap<string,string>::iterator tit = tid_nodes.find(tid);

                if( ignore_tid_tags )
                    tit = tid_nodes.begin();

                if(tit != tid_nodes.end())
                {
                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_to_align = tit->second;

                    multimap<string,hit>::iterator ith = exonerate_hits.find(tit->second);
                    if(ith != exonerate_hits.end())
                    {
                        hit h = ith->second;

                        reads->at(i).local_qstart = h.q_start;
                        reads->at(i).local_qend = h.q_end;
                        reads->at(i).local_tstart = h.t_start;
                        reads->at(i).local_tend = h.t_end;
                    }

                    if(Settings_handle::st.is("placement-file"))
                    {
                        pl_output<<reads->at(i).name<<" "<<tit->second<<" TID="<<tid<<endl;
                    }

                }
            }


            // has TID and matching nodes, or exhaustive search

            else
            {
                double best_score = -HUGE_VAL;
                string best_node = root->get_name();

                map<string,Node*> nodes;
                root->get_all_nodes(&nodes);

                multimap<string,string>::iterator tit;
                int matching_nodes = 0;

                if( ignore_tid_tags )
                {
                    tit = tid_nodes.begin();
                    matching_nodes = tid_nodes.size();
                }
                else
                {
                    tit = tid_nodes.find(tid);
                    matching_nodes = tid_nodes.count(tid);
                }

                if(tit != tid_nodes.end())
                {

                    stringstream ss;
                    ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                    Log_output::write_out(ss.str(),2);

                    while(tit != tid_nodes.end() && matching_nodes>0)
                    {
                        map<string,Node*>::iterator nit = nodes.find(tit->second);
                        double score = this->read_match_score( nit->second, &reads->at(i), mf, best_score);

                        stringstream ss;
                        ss<<tit->second<<" with score "<<score<<" (simple p-distance)\n";
                        Log_output::write_out(ss.str(),2);

                        if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                        {
                            best_score = score;
                            best_node.append(" "+tit->second);
                        }
                        else if(score>best_score)
                        {
                            best_score = score;
                            best_node = tit->second;
                        }
                        tit++;
                        matching_nodes--;
                    }
                }
                if(best_score<0.05)
                {
                    if(Settings_handle::st.is("align-bad-reads-at-root"))
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                        reads->at(i).node_to_align = root->get_name();

                        if(Settings_handle::st.is("placement-file"))
                        {
                            pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (bad) "<<endl;
                        }
                    }
                    else
                    {
                        Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                        reads->at(i).node_to_align = "discarded_read";
                    }
                }
                else
                {
                    stringstream ss;
                    ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                    Log_output::write_out(ss.str(),2);

                    reads->at(i).node_score = best_score;
                    reads->at(i).node_to_align = best_node;

                    if(Settings_handle::st.is("placement-file"))
                    {
                        pl_output<<reads->at(i).name<<" "<<best_node<<" TID="<<tid<<endl;
                    }

                }
            }
        }


        // no TID, aligning at root

        else
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") has no tid. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = root->get_name();

            if(Settings_handle::st.is("placement-file"))
            {
                pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID=NULL"<<endl;
            }
        }
    }

    if(Settings_handle::st.is("placement-file"))
    {
        pl_output.close();
    }
}


void Reads_aligner::find_nodes_for_all_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{

    ofstream pl_output;

    if(Settings_handle::st.is("placement-file"))
    {
        string fname = Settings_handle::st.get("placement-file").as<string>();
        pl_output.open(fname.append(".tsv").c_str());
    }


    multimap<string,string> all_tid_nodes;
    bool ignore_tid_tags = true;

    this->get_target_node_names(root,&all_tid_nodes,&ignore_tid_tags);

    //
    // Call Exonerate to reduce the search space
    //
    map<string, multimap<string,hit> > exonerate_hits;

    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }

    else if(all_tid_nodes.size()>0)
        er.all_local_alignments(root,reads,&all_tid_nodes,&exonerate_hits, true, ignore_tid_tags);


    for(int i=0;i<(int)reads->size();i++)
    {

        reads->at(i).node_score = -1.0;

        string tid = reads->at(i).tid;
        if( ignore_tid_tags )
            tid = "<empty>";

        multimap<string,string> tid_nodes;

        // Discarded by Exonerate
        //
        if(reads->at(i).node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("exhaustive-placement"))
            {
                this->get_target_node_names(root,&all_tid_nodes,&ignore_tid_tags);
            }

            else
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
                Log_output::write_out(ss.str(),2);

                continue;
            }
        }
        else
        {
            map<string,multimap<string,hit> >::iterator iter = exonerate_hits.find(reads->at(i).name);

            if( iter != exonerate_hits.end() )
            {
                multimap<string,hit>::iterator iter2 = iter->second.begin();

                for( ;iter2 != iter->second.end(); iter2++ )
                    tid_nodes.insert( make_pair<string,string>(iter2->first,iter2->second.node) );
            }

            if(has_exonerate && tid_nodes.size()>1)
            {
                map<string, hit> exonerate_gapped_hits;
                er.local_alignment(root,&reads->at(i),&tid_nodes,&exonerate_gapped_hits,false,ignore_tid_tags);

                tid_nodes.clear();
                map<string, hit>::iterator iter2 = exonerate_gapped_hits.begin();
                if( iter2 != exonerate_gapped_hits.end() )
                {
                    tid_nodes.insert( make_pair<string,string>(iter2->first,iter2->second.node) );
                }
            }
        }



        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping read: '"+reads->at(i).name+"'",0);


        int matches = tid_nodes.size();

        // has TID but no matching node

        if(matches == 0)
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = root->get_name();

            if(Settings_handle::st.is("placement-file"))
            {
                pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (no match)"<<endl;
            }

        }


        // has only one matching node and no need for ranking

        else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
        {
            multimap<string,string>::iterator tit = tid_nodes.begin();

            if(tit != tid_nodes.end())
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_to_align = tit->second;

                if(Settings_handle::st.is("placement-file"))
                {
                    pl_output<<reads->at(i).name<<" "<<tit->second<<" TID="<<tid<<endl;
                }

            }
        }


        // has TID and matching nodes, or exhaustive search
        //
        else
        {
            double best_score = -HUGE_VAL;
            string best_node = root->get_name();

            map<string,Node*> nodes;
            root->get_all_nodes(&nodes);

            multimap<string,string>::iterator tit = tid_nodes.begin();
            int matching_nodes = tid_nodes.size();


            if(tit != tid_nodes.end())
            {

                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                Log_output::write_out(ss.str(),2);

                while(tit != tid_nodes.end() && matching_nodes>0)
                {
                    map<string,Node*>::iterator nit = nodes.find(tit->second);
                    double score = this->read_match_score( nit->second, &reads->at(i), mf, best_score);

                    stringstream ss;
                    ss<<tit->second<<" with score "<<score<<" (simple p-distance)\n";
                    Log_output::write_out(ss.str(),2);

                    if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                    {
                        best_score = score;
                        best_node.append(" "+tit->second);
                    }
                    else if(score>best_score)
                    {
                        best_score = score;
                        best_node = tit->second;
                    }
                    tit++;
                    matching_nodes--;
                }
            }

            if(best_score<0.05)
            {
                if(Settings_handle::st.is("align-bad-reads-at-root"))
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                    reads->at(i).node_to_align = root->get_name();

                    if(Settings_handle::st.is("placement-file"))
                    {
                        pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (bad) "<<endl;
                    }
                }
                else
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                    reads->at(i).node_to_align = "discarded_read";
                }
            }
            else
            {
                stringstream ss;
                ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_score = best_score;
                reads->at(i).node_to_align = best_node;

                if(Settings_handle::st.is("placement-file"))
                {
                    pl_output<<reads->at(i).name<<" "<<best_node<<" TID="<<tid<<endl;
                }

            }
        }
    }


    if(Settings_handle::st.is("placement-file"))
    {
        pl_output.close();
    }
}

//

void Reads_aligner::find_nodes_for_all_reads_together(Node *root, vector<Fasta_entry> *reads, Model_factory *mf)
{
    ofstream pl_output;

    if(Settings_handle::st.is("placement-file"))
    {
        string fname = Settings_handle::st.get("placement-file").as<string>();
        pl_output.open(fname.append(".tsv").c_str());
    }

    multimap<string,string> all_tid_nodes;
    bool ignore_tid_tags = true;

    this->get_target_node_names(root,&all_tid_nodes,&ignore_tid_tags);

    //
    // Call Exonerate to reduce the search space
    //
    Exonerate_queries er;
    bool has_exonerate = true;
    if(!er.test_executable())
    {
        has_exonerate = false;
        Log_output::write_out("The executable for Exonerate not found! The fast placement search not used!",0);
    }

    map<string, multimap<string,hit> > exonerate_hits;

    if( has_exonerate && all_tid_nodes.size()>1)
    {
        er.all_local_alignments(root,reads,&all_tid_nodes,&exonerate_hits, true, ignore_tid_tags);
    }

    //
    // Map one read at time
    //
    for(int i=0;i<(int)reads->size();i++)
    {

        reads->at(i).node_score = -1.0;

        string tid = reads->at(i).tid;
        if( ignore_tid_tags )
            tid = "<empty>";

        multimap<string,string> tid_nodes;

        // Discarded by Exonerate
        //
        if(reads->at(i).node_to_align == "discarded_read")
        {
            if(Settings_handle::st.is("exhaustive-placement"))
            {
                this->get_target_node_names(root,&tid_nodes,&ignore_tid_tags);
            }

            else
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" was discarded by Exonerate.\n";
                Log_output::write_out(ss.str(),2);

                continue;
            }
        }
        else
        {
            map<string,multimap<string,hit> >::iterator iter = exonerate_hits.find(reads->at(i).name);

            if( iter != exonerate_hits.end() )
            {
                multimap<string,hit>::iterator iter2 = iter->second.begin();

                for( ;iter2 != iter->second.end(); iter2++ )
                    tid_nodes.insert( make_pair<string,string>(iter2->first,iter2->second.node) );
            }
        }



        Log_output::write_msg("("+Log_output::itos(i+1)+"/"+Log_output::itos(reads->size())+") mapping read: '"+reads->at(i).name+"'",0);


        int matches = tid_nodes.size();

        // has TID but no matching node

        if(matches == 0)
        {
            stringstream ss;
            ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" has no matching node. Aligned to root.\n";
            Log_output::write_out(ss.str(),2);

            reads->at(i).node_to_align = root->get_name();

            if(Settings_handle::st.is("placement-file"))
            {
                pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (no match)"<<endl;
            }

        }


        // has only one matching node and no need for ranking

        else if(matches == 1 && !Settings_handle::st.is("rank-reads-for-nodes") )
        {
            multimap<string,string>::iterator tit = tid_nodes.begin();

            if(tit != tid_nodes.end())
            {
                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" ("<<i+1<<"/"<<reads->size()<<") with the tid "<<tid<<" only matches the node "<<tit->second<<"."<<endl;
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_to_align = tit->second;

                if(Settings_handle::st.is("placement-file"))
                {
                    pl_output<<reads->at(i).name<<" "<<tit->second<<" TID="<<tid<<endl;
                }

            }
        }


        // has TID and matching nodes, or exhaustive search
        //
        else
        {
            double best_score = -HUGE_VAL;
            string best_node = root->get_name();

            map<string,Node*> nodes;
            root->get_all_nodes(&nodes);

            multimap<string,string>::iterator tit = tid_nodes.begin();
            int matching_nodes = tid_nodes.size();


            if(tit != tid_nodes.end())
            {

                stringstream ss;
                ss<<"Read "<<reads->at(i).name<<" with TID "<<tid<<" matches "<<matching_nodes<<" nodes.\n";
                Log_output::write_out(ss.str(),2);

                while(tit != tid_nodes.end() && matching_nodes>0)
                {
                    map<string,Node*>::iterator nit = nodes.find(tit->second);
                    double score = this->read_match_score( nit->second, &reads->at(i), mf, best_score);

                    stringstream ss;
                    ss<<tit->second<<" with score "<<score<<" (simple p-distance)\n";
                    Log_output::write_out(ss.str(),2);

                    if(score==best_score && !Settings_handle::st.is("one-placement-only") && !Settings_handle::st.is("exhaustive-placement"))
                    {
                        best_score = score;
                        best_node.append(" "+tit->second);
                    }
                    else if(score>best_score)
                    {
                        best_score = score;
                        best_node = tit->second;
                    }
                    tit++;
                    matching_nodes--;
                }
            }

            if(best_score<0.05)
            {
                if(Settings_handle::st.is("align-bad-reads-at-root"))
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Aligning to root instead.\n",2);

                    reads->at(i).node_to_align = root->get_name();

                    if(Settings_handle::st.is("placement-file"))
                    {
                        pl_output<<reads->at(i).name<<" "<<root->get_name()<<" TID="<<tid<<" (bad) "<<endl;
                    }
                }
                else
                {
                    Log_output::write_out("Best node aligns with less than 5% of identical sites. Read is discarded.\n",2);

                    reads->at(i).node_to_align = "discarded_read";
                }
            }
            else
            {
                stringstream ss;
                ss<<"best node "<<best_node<<" (score "<<best_score<<").\n";
                Log_output::write_out(ss.str(),2);

                reads->at(i).node_score = best_score;
                reads->at(i).node_to_align = best_node;

                if(Settings_handle::st.is("placement-file"))
                {
                    pl_output<<reads->at(i).name<<" "<<best_node<<" TID="<<tid<<endl;
                }

            }
        }
    }


    if(Settings_handle::st.is("placement-file"))
    {
        pl_output.close();
    }
}

double Reads_aligner::read_match_score(Node *node, Fasta_entry *read, Model_factory *mf, float best_score)
{

    double r_dist = Settings_handle::st.get("query-distance").as<float>();

    double org_dist = node->get_distance_to_parent();
    node->set_distance_to_parent(0.001);

    Node * reads_node1 = new Node();
    reads_node1->set_distance_to_parent(r_dist);
    reads_node1->set_name(read->name);
    reads_node1->add_name_comment(read->comment);
    reads_node1->add_sequence( *read, read->data_type, false, true);

    Node * tmpnode = new Node();
    tmpnode->set_name("(tmp)");

    tmpnode->add_left_child(node);
    tmpnode->add_right_child(reads_node1);

    tmpnode->align_sequences_this_node(mf,true);

    node->set_distance_to_parent(org_dist);

    double score = 0;

    if(tmpnode->node_has_sequence_object)
    {

        // For scoring (below)
        Evol_model model = mf->alignment_model(r_dist+0.001,false);

        int matching = 0;
        int aligned = 0;

        float subst_score = 0;
        float max_subst_score_l = 0;
        float max_subst_score_r = 0;

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

                subst_score += model.score(site1->get_state(),site2->get_state());
                max_subst_score_l += model.score(site2->get_state(),site2->get_state());

                aligned++;
            }

            if(site->get_children()->right_index>=0)
            {
                Site *site1 = tmpnode->get_right_child()->get_sequence()->get_site_at(site->get_children()->right_index);
                max_subst_score_r += model.score(site1->get_state(),site1->get_state());
            }

        }

        double score_s = (double) matching/ (double) reads_node1->get_sequence()->sites_length();
        double score_l = (double) subst_score/ (double) max_subst_score_l;
        double score_r = (double) subst_score/ (double) max_subst_score_r;

        score = score_r;

        if(Settings_handle::st.is("use-identity-score"))
            score = score_s;
        else if(Settings_handle::st.is("use-target-normalised-score"))
            score = score_l;

        if(score>best_score)
        {
            read->node_start_pos1 = node_start_pos1;
            read->node_end_pos1 =   node_end_pos1;
            read->node_start_pos2 = node_start_pos2;
            read->node_end_pos2 =   node_end_pos2;
        }
    }

    tmpnode->has_left_child(false);
    delete tmpnode;

    return score;
}

void Reads_aligner::remove_target_overlapping_identical_reads(vector<Fasta_entry> *reads, Model_factory *mf)
{
    Log_output::write_header("Aligning reads: remove reads mapped at overlapping positions",0);

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

                if(Settings_handle::st.is("pair-end")) // remove the midpoint marker
                {
                    r1_length--;
                    r2_length--;
                }

                if( matching == r2_length )
                {
                    Log_output::write_out("Read "+ri2->name+" is fully embedded in read "+ri1->name+" and overlapping sites are identical.  Read "+ri2->name+" is deleted.\n",2);
                    reads->erase(ri2);
                }
                else if( matching == r1_length )
                {
                    Log_output::write_out("Read "+ri1->name+" is fully embedded in read "+ri2->name+" and overlapping sites are identical.  Read "+ri1->name+" is deleted.\n",2);
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

void Reads_aligner::remove_target_overlapping_reads(vector<Fasta_entry> *reads)
{
    Log_output::write_header("Aligning reads: remove reads mapped at overlapping positions",0);

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
                    Log_output::write_out("Read "+ri2->name+" is fully embedded in read "+ri1->name+".  Read "+ri2->name+" is deleted.\n",2);
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
                    Log_output::write_out("Read "+ri2->name+" is fully embedded in read "+ri1->name+".  Read "+ri2->name+" is deleted.\n",2);
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

void Reads_aligner::align_two_reads(Node *node, Fasta_entry *ri1, Fasta_entry *ri2, Model_factory *mf)
{
    double r_dist = Settings_handle::st.get("query-distance").as<float>();

    node->set_name("read pair");

    Node * reads_node1 = new Node();
    reads_node1->set_distance_to_parent(r_dist);
    reads_node1->set_name(ri1->name);
    reads_node1->add_name_comment(ri1->comment);
    reads_node1->add_sequence( *ri1, Model_factory::dna);
    node->add_right_child(reads_node1);

    Node * reads_node2 = new Node();
    reads_node2->set_distance_to_parent(r_dist);
    reads_node2->set_name(ri2->name);
    reads_node2->add_name_comment(ri2->comment);
    reads_node2->add_sequence( *ri2, Model_factory::dna);
    node->add_left_child(reads_node2);

    node->align_sequences_this_node(mf,true);

}

int Reads_aligner::reads_pairwise_matching_sites(Node *node)
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

void Reads_aligner::remove_overlapping_reads(vector<Fasta_entry> *reads, Model_factory *mf)
{

    Log_output::write_header("Aligning reads: remove overlapping reads",0);

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
                Log_output::write_out("Reads "+ri1->name+" and "+ri2->name+" are identical.  Read "+ri2->name+" is deleted.\n",2);
                reads->erase(ri2);
            }
            else if(matching == r1_length)
            {
                Log_output::write_out("Read "+ri1->name+" is fully embedded in read "+ri2->name+".  Read "+ri1->name+" is deleted.\n",2);
                reads->erase(ri1);
                ri1--;
                ri2 = reads->end();
            }
            else if(matching == r2_length)
            {
                Log_output::write_out("Read "+ri2->name+" is fully embedded in read "+ri1->name+".  Read "+ri2->name+" is deleted.\n",2);
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
