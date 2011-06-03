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

#ifndef MATRIX_POINTER_H
#define MATRIX_POINTER_H

#include <cmath>
#include "utils/model_factory.h"

namespace ppa {

using namespace std;
using namespace ppa;

struct Matrix_pointer
{
    double score;
    double full_score;
    double fwd_score;
    double bwd_score;
    int x_ind;
    int y_ind;
    int x_edge_ind;
    int y_edge_ind;
    int matrix;
    int path_index;

    Matrix_pointer() : score(-HUGE_VAL), full_score(0), fwd_score(0), bwd_score(0), x_ind(-1), y_ind(-1),
                        x_edge_ind(-1), y_edge_ind(-1), matrix(-1), path_index(-1) {}
    Matrix_pointer(double s,int x, int y, int m) : score(s), full_score(0), fwd_score(0), bwd_score(0),
                        x_ind(x), y_ind(y), x_edge_ind(-1), y_edge_ind(-1), matrix(m), path_index(-1) {}
};

struct Path_pointer
{

    Matrix_pointer mp;
    bool real_site;

    float branch_length_increase;
    int branch_count_increase;

    Path_pointer(Matrix_pointer m, bool real, float l, int i) : mp(m), real_site(real), branch_length_increase(l), branch_count_increase(i) {}
    Path_pointer(Matrix_pointer m, bool real) : mp(m), real_site(real), branch_length_increase(0), branch_count_increase(0) {}
    Path_pointer(bool real) : real_site(real), branch_length_increase(0), branch_count_increase(0) {}

};

struct Edge_history
{
    int left_real_site_index;
    int right_real_site_index;
    int left_skip_site_index;
    int right_skip_site_index;

    int real_site_index;
    int match_site_index;
    int path_state;

    Edge_history(int lr,int rr,int ls=-1,int rs=-1) : left_real_site_index(lr), right_real_site_index(rr),
                                                      left_skip_site_index(ls), right_skip_site_index(rs),
                                                      real_site_index(-1), match_site_index(0), path_state(-1) {}
};


}
#endif // MATRIX_POINTER_H

#ifndef BASIC_ALIGNMENT_H
#define BASIC_ALIGNMENT_H

#include "utils/evol_model.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
//#include "main/graph_reconstruction.h"
#include "boost/multi_array.hpp"
#include <string>

namespace ppa {

using namespace std;
using namespace ppa;


class Basic_alignment
{
protected:
    enum Matrix_pt {x_mat,y_mat,m_mat};
    enum Gap_type {normal_gap,end_gap,pair_break_gap};

    Sequence *left;
    Sequence *right;
    Sequence *ancestral_sequence;

    float left_branch_length;
    float right_branch_length;

    Evol_model *model;
//    string full_char_alphabet;

    vector<Path_pointer> path;

    // parameters that may be set by program arguments
    //
    float ins_del_ratio;   // (not used so far)
    float del_ins_ratio;   // (not used so far)

    float max_allowed_skip_distance;
    int   max_allowed_skip_branches;
    int   max_allowed_match_skip_branches;
    float branch_skip_weight;
    float branch_skip_probability;
    bool  weighted_branch_skip_penalty;

    bool no_terminal_edges;
    bool pair_end_reads;
    bool edges_for_skipped_flanked_by_gaps;
    int x_length;
    int y_length;
    int x_read1_length;
    int y_read1_length;

    bool compute_full_score;
    bool weight_edges;
    double bwd_full_probability; // for control

    void build_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path,bool is_reads_sequence=false);

    void create_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path,bool is_reads_sequence);
    void create_ancestral_edges(Sequence *sequence);
    void check_skipped_boundaries(Sequence *sequence);

    void delete_edge_range(Sequence *sequence,int edge_ind,int skip_start_site);

    void transfer_child_edge(Sequence *sequence, Edge *child, vector<int> *child_index, float branch_length,
                             bool connects_neighbour_site = false, bool adjust_posterior_weight = true, float branch_weight = 1.0);
    void transfer_child_edge(Sequence *sequence, Edge edge, Edge *child, float branch_length,
                             bool connects_neighbour_site = false, bool adjust_posterior_weight = true, float branch_weight = 1.0);

    /*********************************/

    void compute_site_consensus(Site *site,Sequence *left,int l_pos,Sequence *right,int r_pos,bool is_dna)
    {
        if(!is_dna)
        {
            if(l_pos>=0)
            {
                if(!left->is_terminal_sequence())
                    site->set_sumAmino(left->get_site_at(l_pos)->get_sumAmino());
                else if(left->is_read_sequence())
                    site->set_sumAmino(1);
            }
            if(r_pos>=0)
            {
                if(!right->is_terminal_sequence())
                    site->set_sumAmino(right->get_site_at(r_pos)->get_sumAmino());
                else if(right->is_read_sequence())
                    site->set_sumAmino(1);
            }

            return;
        }

        int lsA = 0; int lsC = 0; int lsG = 0; int lsT = 0;
        int rsA = 0; int rsC = 0; int rsG = 0; int rsT = 0;

        if(l_pos>=0)
        {
            if(!left->is_terminal_sequence())
            {
                lsA += left->get_site_at(l_pos)->get_sumA();
                lsC += left->get_site_at(l_pos)->get_sumC();
                lsG += left->get_site_at(l_pos)->get_sumG();
                lsT += left->get_site_at(l_pos)->get_sumT();
            }
            else if(left->is_read_sequence())
            {
                int s = left->get_site_at(l_pos)->get_state();
                if(s == 0)
                    lsA += 1;
                else if(s == 1)
                    lsC += 1;
                else if(s == 2)
                    lsG += 1;
                else if(s == 3)
                    lsT += 1;
                else if(s == 4) {
                    lsA += 1;
                    lsG += 1;    }
                else if(s == 5) {
                    lsC += 1;
                    lsT += 1;    }
                else if(s == 6) {
                    lsA += 1;
                    lsC += 1;    }
                else if(s == 7) {
                    lsG += 1;
                    lsT += 1;    }
                else if(s == 8) {
                    lsA += 1;
                    lsT += 1;    }
                else if(s == 9) {
                    lsC += 1;
                    lsG += 1;    }
                else if(s == 10) {
                    lsC += 1;
                    lsG += 1;
                    lsT += 1;    }
                else if(s == 11) {
                    lsA += 1;
                    lsG += 1;
                    lsT += 1;    }
                else if(s == 12) {
                    lsA += 1;
                    lsC += 1;
                    lsT += 1;    }
                else if(s == 13) {
                    lsA += 1;
                    lsC += 1;
                    lsG += 1;    }
                else if(s == 14) {
                    lsA += 1;
                    lsC += 1;
                    lsG += 1;
                    lsT += 1;    }
                else
                    cout<<"compute_site_consensus: no such option (l)\n";
            }
        }

        if(r_pos>=0)
        {
            if(!right->is_terminal_sequence())
            {
                rsA += right->get_site_at(r_pos)->get_sumA();
                rsC += right->get_site_at(r_pos)->get_sumC();
                rsG += right->get_site_at(r_pos)->get_sumG();
                rsT += right->get_site_at(r_pos)->get_sumT();
            }
            else if(right->is_read_sequence())
            {
                int s = right->get_site_at(r_pos)->get_state();
                if(s == 0)
                    rsA += 1;
                else if(s == 1)
                    rsC += 1;
                else if(s == 2)
                    rsG += 1;
                else if(s == 3)
                    rsT += 1;
                else if(s == 4) {
                    rsA += 1;
                    rsG += 1;    }
                else if(s == 5) {
                    rsC += 1;
                    rsT += 1;    }
                else if(s == 6) {
                    rsA += 1;
                    rsC += 1;    }
                else if(s == 7) {
                    rsG += 1;
                    rsT += 1;    }
                else if(s == 8) {
                    rsA += 1;
                    rsT += 1;    }
                else if(s == 9) {
                    rsC += 1;
                    rsG += 1;    }
                else if(s == 10) {
                    rsC += 1;
                    rsG += 1;
                    rsT += 1;    }
                else if(s == 11) {
                    rsA += 1;
                    rsG += 1;
                    rsT += 1;    }
                else if(s == 12) {
                    rsA += 1;
                    rsC += 1;
                    rsT += 1;    }
                else if(s == 13) {
                    rsA += 1;
                    rsC += 1;
                    rsG += 1;    }
                else if(s == 14) {
                    rsA += 1;
                    rsC += 1;
                    rsG += 1;
                    rsT += 1;    }
                else
                    cout<<"compute_site_consensus: no such option (r)\n";
            }
        }

        if(lsA+lsC+lsG+lsT+rsA+rsC+rsG+rsT>0)
        {
            int sA = lsA+rsA;
            int sC = lsC+rsC;
            int sG = lsG+rsG;
            int sT = lsT+rsT;

            site->set_sumA(sA);
            site->set_sumC(sC);
            site->set_sumG(sG);
            site->set_sumT(sT);

            if(Settings_handle::st.is("use-consensus"))
            {
                if(sA>sC && sA>sG && sA>sT)
                    site->set_state(0);
                else if(sC>sA && sC>sG && sC>sT)
                    site->set_state(1);
                else if(sG>sA && sG>sC && sG>sT)
                    site->set_state(2);
                else if(sT>sA && sT>sC && sT>sG)
                    site->set_state(3);
                else if(sA>sC && sA==sG && sA>sT)
                    site->set_state(4);
                else if(sC>sA && sC>sG && sC==sT)
                    site->set_state(5);
                else if(sA==sC && sA>sG && sA>sT)
                    site->set_state(6);
                else if(sG>sA && sG>sC && sG==sT)
                    site->set_state(7);
                else if(sA>sC && sA>sG && sA==sT)
                    site->set_state(8);
                else if(sC>sA && sC==sG && sC>sT)
                    site->set_state(9);
                else if(sC>sA && sC==sG && sC==sT)
                    site->set_state(10);
                else if(sA>sC && sA==sG && sA==sT)
                    site->set_state(11);
                else if(sA==sC && sA>sG && sA==sT)
                    site->set_state(12);
                else if(sA==sC && sA==sG && sA>sT)
                    site->set_state(13);
                else if(sA==sC && sA==sG && sA==sT)
                    site->set_state(14);
                else
                    cout<<"compute_site_consensus: no such option (s)"<<sA<<" "<<sC<<" "<<sG<<" "<<sT<<"\n";
            }
        }
    }

    /*********************************/

    void debug_msg(std::string msg,int noise_level)
    {
        if(Settings::noise>noise_level)
            cout<<msg<<endl;
    }

    std::string itos(int i) // convert int to string
    {
        std::stringstream s;
        s << i;
        return s.str();
    }

    std::string ftos(float f) // convert float to string
    {
        std::stringstream s;
        s << f;
        return s.str();
    }

    /*********************************/

    bool first_is_bigger(double a,double b)
    {
        if (a==-HUGE_VAL && b==-HUGE_VAL)
            return false;
        else if (a>b)
            return true;
        else if (a<b)
            return false;
        else
            return false;
//            if ((double)rand()/(double)RAND_MAX>0.5)
//                return true;
        return false;
    }


    /********************************************/

    double get_log_edge_weight(Edge *edge) { return (*this.*log_edge_weight)(edge); }
    double (ppa::Basic_alignment::*log_edge_weight)(Edge *edge);

    double get_edge_weight(Edge *edge) { return (*this.*edge_weight)(edge); }
    double (ppa::Basic_alignment::*edge_weight)(Edge *edge);

    double edge_posterior_weight(Edge *edge) { return edge->get_posterior_weight(); }
    double edge_log_posterior_weight(Edge *edge) { return edge->get_log_posterior_weight(); }

    double edge_equal_weight(Edge *) const { return 1.0; }
    double edge_log_equal_weight(Edge *) const { return 0.0; }

    /********************************************/

    double get_transformed_edge_weight(double w) { return (*this.*transform_edge_weight)(w); }
    double (ppa::Basic_alignment::*transform_edge_weight)(double w);

    double square_root_edge_weight(double w) { return sqrt( w );}
    double cube_root_edge_weight(double w) { return exp((1.0/3.0)*(log(w)));}
    double plain_edge_weight(double w) { return w;}

    /********************************************/

    float get_log_gap_open_penalty(int prev_site, bool is_x_matrix)
    {
        if(no_terminal_edges)
        {
            if(prev_site==0)
            {
                return 0;
            }

            if(pair_end_reads)
            {
                if(is_x_matrix && prev_site == x_read1_length)
                {
                    return 0;
                }
                else if(!is_x_matrix && prev_site == y_read1_length)
                {
                    return 0;
                }
            }
        }

        return model->log_gap_open();
    }

    float get_log_gap_close_penalty(int this_site, bool is_x_matrix)
    {
        if(no_terminal_edges)
        {
            if(is_x_matrix && this_site==x_length)
            {
                return 0;
            }
            else if(!is_x_matrix && this_site==y_length)
            {
                return 0;
            }

            if(pair_end_reads)
            {
                if(is_x_matrix && this_site == x_read1_length+1)
                {
                    return 0;
                }
                else if(!is_x_matrix && this_site == y_read1_length+1)
                {
                    return 0;
                }
            }
        }

        return model->log_gap_close();
    }

    /********************************************/

    void set_basic_settings()
    {
        del_ins_ratio = 1;//model->del_prob/model->ins_prob;
        ins_del_ratio = 1;//model->ins_prob/model->del_prob;

        edges_for_skipped_flanked_by_gaps = false;

        weighted_branch_skip_penalty = false; // by default, use penalty *per node*

        max_allowed_skip_distance = 0.5;
        if(Settings_handle::st.is("branch-length-confirm-insertion"))
            max_allowed_skip_distance = Settings_handle::st.get("branch-length-confirm-insertion").as<float>();

        max_allowed_skip_branches = 10;
        if(Settings_handle::st.is("any-skips-confirm-insertion"))
            max_allowed_skip_branches = Settings_handle::st.get("any-skips-confirm-insertion").as<int>();

        max_allowed_match_skip_branches = 5;
        if(Settings_handle::st.is("match-skips-confirm-insertion"))
            max_allowed_match_skip_branches = Settings_handle::st.get("match-skips-confirm-insertion").as<int>();

        branch_skip_weight = 1;
        if(Settings_handle::st.is("branch-skip-weight-per-distance"))
        {
            branch_skip_weight = Settings_handle::st.get("branch-skip-weight-per-distance").as<float>();
            weighted_branch_skip_penalty = true;
        }

        branch_skip_probability = 0.2;
        if(Settings_handle::st.is("branch-skip-penalty-per-branch"))
        {
            branch_skip_probability = Settings_handle::st.get("branch-skip-penalty-per-branch").as<float>();
            weighted_branch_skip_penalty = false;
        }

        weight_edges = false;
        if( Settings_handle::st.is("weight-sampled-edges") && Settings_handle::st.get("sample-additional-paths").as<int>() > 0)
            weight_edges = true;

        compute_full_score = false;
        if( Settings_handle::st.is("full-probability") ||
            Settings_handle::st.is("mpost-posterior-plot-file") ||
            Settings_handle::st.is("sample-path") ||
            Settings_handle::st.get("sample-additional-paths").as<int>() > 0 )
            compute_full_score = true;

        no_terminal_edges = false;
        if( Settings_handle::st.is("no-terminal-edges") )
            no_terminal_edges = true;

        pair_end_reads = false;

        cout << noshowpos;
    }

    void set_reads_alignment_settings()
    {

        no_terminal_edges = true;

        max_allowed_skip_distance = 5;
        max_allowed_skip_branches = 50000;
        max_allowed_match_skip_branches = 50000;

        branch_skip_weight = 1;
        branch_skip_probability = 1;

        if(Settings_handle::st.is("pair-end"))
            pair_end_reads = true;
    }

    void set_reference_alignment_settings()
    {

        max_allowed_skip_distance = 5;
        max_allowed_skip_branches = 50000;
        max_allowed_match_skip_branches = 50000;

    }

    /********************************************/

    void mark_no_gap_penalty_sites(Sequence *left, Sequence *right)
    {

        x_length = left->sites_length();
        y_length = right->sites_length();
        x_read1_length = -1;
        y_read1_length = -1;

        if(pair_end_reads)
        {
            // this is a clumsy way to transfer the information of the break point
            // in a combined pair-end read but has to for now
            for(int i=0;i<left->sites_length();i++)
            {
                if( left->get_site_at(i)->get_site_type() == Site::break_start_site )
                {
                    x_read1_length = i;
                    left->get_site_at(i)->set_site_type( Site::real_site );
                }

                if( left->get_site_at(i)->get_site_type() == Site::break_stop_site )
                {
                    left->get_site_at(i)->set_site_type( Site::real_site );
                    break;
                }
            }

            for(int i=0;i<right->sites_length();i++)
            {
                if( right->get_site_at(i)->get_site_type() == Site::break_start_site )
                {
                    y_read1_length = i;
                    right->get_site_at(i)->set_site_type( Site::real_site );
                }

                if( right->get_site_at(i)->get_site_type() == Site::break_stop_site )
                {
                    right->get_site_at(i)->set_site_type( Site::real_site );
                    break;
                }
            }
        }

    }

    /********************************************/

    string print_pairwise_alignment(vector<Site> *sites);

    void debug_print_input_sequences(int noise_level)
    {
        if(Settings::noise>noise_level)
            print_input_sequences();
    }

    void print_input_sequences()
    {

        vector<string> *full_char_alphabet = Model_factory::get_dna_full_character_alphabet();

        if(model->get_data_type() == Model_factory::protein)
            full_char_alphabet = Model_factory::get_protein_full_character_alphabet();
        else if(model->get_data_type() == Model_factory::dna && Settings_handle::st.is("codons"))
            full_char_alphabet = Model_factory::get_codon_full_character_alphabet();

        cout<<"sequences:"<<endl<<" ";
        for(int i=1;i<left->sites_length()-1;i++)
            cout<<Model_factory::get_ancestral_character_alphabet_at(left->get_site_at(i)->get_state());
        cout<<endl<<" ";
        for(int i=1;i<right->sites_length()-1;i++)
            cout<<Model_factory::get_ancestral_character_alphabet_at(right->get_site_at(i)->get_state());
        cout<<endl;

        if(Settings::noise>4)
        {
            cout<<"\nLEFT";left->print_sequence(left->get_sites());
            cout<<"\nRIGHT";right->print_sequence(right->get_sites());
            cout<<endl;
        }
    }

    void print_path(vector<Path_pointer> *path)
    {
        cout<<endl;
        for(unsigned int i=0;i<path->size();i++)
        {
            Path_pointer *tsite =  &path->at(i);
            if(tsite->real_site)
                cout<<i<<" r: ";
            else
                cout<<i<<" s: ";

            cout<<tsite->mp.x_ind<<" "<<tsite->mp.y_ind<<" "<<tsite->mp.matrix<<": "<<tsite->mp.score<<" "<<
                   fixed<<log(tsite->mp.full_score)<<" "<<fixed<<log(tsite->mp.bwd_score)<<endl;
        }
    }

    /********************************************/

public:
    Basic_alignment();

    Sequence* get_simple_sequence() { return ancestral_sequence; }

};

}
#endif // BASIC_ALIGNMENT_H
