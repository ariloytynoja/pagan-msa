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

#ifndef SIMPLE_ALIGNMENT_H
#define SIMPLE_ALIGNMENT_H

#include "utils/evol_model.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
#include "boost/multi_array.hpp"
#include <string>

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

    Matrix_pointer() : score(-HUGE_VAL), full_score(0), fwd_score(0), bwd_score(0), x_ind(-1), y_ind(-1),
                        x_edge_ind(-1), y_edge_ind(-1), matrix(-1) {}
    Matrix_pointer(double s,int x, int y, int m) : score(s), full_score(0), fwd_score(0), bwd_score(0),
                        x_ind(x), y_ind(y), x_edge_ind(-1), y_edge_ind(-1), matrix(m) {}
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

class Simple_alignment
{
    enum Matrix_pt {x_mat,y_mat,m_mat};
    enum Gap_type {normal_gap,end_gap,pair_break_gap};

    typedef boost::multi_array<Matrix_pointer, 2> align_array;
    typedef align_array::index index;

    align_array *match;
    align_array *xgap;
    align_array *ygap;

    typedef boost::multi_array_types::index_range range;
    align_array::index_gen indices;
    typedef align_array::array_view<1>::type align_slice;

    vector<Matrix_pointer> *mvectp;
    vector<Matrix_pointer> *xvectp;
    vector<Matrix_pointer> *yvectp;

    vector<int> *left_child_site_to_path_index_p;
    vector<int> *right_child_site_to_path_index_p;

    vector<int> *left_child_site_to_last_path_index_p;
    vector<int> *right_child_site_to_last_path_index_p;

    vector<int> *path_to_left_child_site_index_p;
    vector<int> *path_to_right_child_site_index_p;

    vector<int> *opposite_index_left_child_p;
    vector<int> *opposite_index_right_child_p;

    Sequence *left;
    Sequence *right;
    Sequence *ancestral_sequence;

    float left_branch_length;
    float right_branch_length;

    Evol_model *model;
    string full_char_alphabet;

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

    void initialise_array_corner();
    void initialise_array_corner_bwd();
    void initialise_full_arrays(int ll,int rl);

    void compute_fwd_scores(int i,int j);
    void compute_known_fwd_scores(int i,int j,int mat);
    void compute_bwd_full_score(int i,int j);
    void compute_posterior_score(int i,int j,double full_score);

    void backtrack_new_path(vector<Path_pointer> *path,Path_pointer pp);
    void backtrack_new_vector_path(vector<Path_pointer> *path,Path_pointer fp,vector<Matrix_pointer> *simple_path);
    void sample_new_path(vector<Path_pointer> *path,Path_pointer pp);
    void build_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path);

    void create_ancestral_sequence(Sequence *sequence,vector<Path_pointer> *path);
    void create_ancestral_edges(Sequence *sequence);
    void check_skipped_boundaries(Sequence *sequence);

    void delete_edge_range(Sequence *sequence,int edge_ind,int skip_start_site);

    void transfer_child_edge(Sequence *sequence, Edge *child, vector<int> *child_index, float branch_length,
                             bool connects_neighbour_site = false, bool adjust_posterior_weight = true, float branch_weight = 1.0);
    void transfer_child_edge(Sequence *sequence, Edge edge, Edge *child, float branch_length,
                             bool connects_neighbour_site = false, bool adjust_posterior_weight = true, float branch_weight = 1.0);

    /*********************************/

//    void iterate_bwd_edges_for_gap(Site * site,align_slice *x_slice,align_slice *y_slice,align_slice *m_slice,
//                                   Matrix_pointer *max,bool is_x_matrix, bool is_edge_cell = false);
    void iterate_bwd_edges_for_gap(Site * site,align_slice *x_slice,align_slice *y_slice,align_slice *m_slice,
                                   Matrix_pointer *max,bool is_x_matrix, int gap_type = Simple_alignment::normal_gap);
    void iterate_bwd_edges_for_match(Site * left_site,Site * right_site,Matrix_pointer *max);
    void iterate_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max);
    void iterate_known_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max, int matrix);

    void iterate_fwd_edges_for_gap(Site * site,align_slice *g_slice,
                                   Matrix_pointer *max_s,Matrix_pointer *max_d,Matrix_pointer *max_m);
    void iterate_fwd_edges_for_match(Site * left_site,Site * right_site,
                                     Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);

    /*********************************/

    void iterate_bwd_edges_for_known_gap(Site * left_site,Site * right_site,vector<Matrix_pointer> *z_slice,vector<Matrix_pointer> *w_slice,
                                                     vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type,bool alignment_end=false);

    void iterate_bwd_edges_for_known_double_gap(Site * site,vector<Matrix_pointer> *z_slice,vector<Matrix_pointer> *w_slice,
                                                            vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type);

    void iterate_bwd_edges_for_known_match(Site * left_site,Site * right_site,Matrix_pointer *max,int prev_mat);

    void iterate_bwd_edges_for_vector_end(Site * left_site,Site * right_site,Matrix_pointer *max,int last_matrix);

    /*********************************/

    void iterate_bwd_edges_for_sampled_gap(int site_index1,int site_index2,Matrix_pointer *bwd_p,bool is_x_matrix);
    void iterate_bwd_edges_for_sampled_match(int left_index,int right_index,Matrix_pointer *bwd_p);

    void iterate_bwd_edges_for_sampled_end_corner(Site * left_site,Site * right_site,Matrix_pointer *sample);

    void add_sample_m_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double m_match = 0);
    void add_sample_x_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double x_match = 0);
    void add_sample_y_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double y_match = 0);

    void add_sample_gap_ext(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_double(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_open(Edge * edge,align_slice *m_slice,Matrix_pointer *bwd_p,bool is_x_matrix);
    void add_sample_gap_close(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix);

    /*********************************/

    void score_m_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_x_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_y_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);

//    void score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix, bool is_edge_cell = false);
    void score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix, int gap_type = Simple_alignment::normal_gap);
    void score_gap_double(Edge *edge,align_slice *w_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_open(Edge *edge,align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_close(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix);


    /*********************************/
    void score_m_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max, bool allow_any=false);
    void score_x_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max);
    void score_y_match_v(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max);

    void score_gap_ext_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix,int gap_type,bool alignment_end = false);
    void score_gap_double_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *w_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_open_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *m_slice,Matrix_pointer *max,bool is_x_matrix,bool alignment_end = false);
    void score_gap_close_v(Edge *left_edge,Edge *right_edge,vector<Matrix_pointer> *z_slice,Matrix_pointer *max,bool is_x_matrix);

    /*********************************/

    void score_gap_ext_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max);
    void score_gap_double_bwd(Edge *edge,align_slice *w_slice,Matrix_pointer *max);
    void score_gap_open_bwd(Edge *edge,align_slice *m_slice,Matrix_pointer *max);

    void score_match_bwd(Edge * left_edge,Edge * right_edge,
                                       Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);

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

    int get_expected_path_length(vector<Path_pointer> *path)
    {
        int m = 0;

        for(unsigned int i=1;i<path->size()-2;i++)
            if(path->at(i).mp.matrix == Simple_alignment::m_mat)
                ++m;

        return (left->sites_length()+right->sites_length()-m);
    }

    void insert_gap_path_pointer(vector<Path_pointer> *path, int i, int j, int matrix,float branch_length)
    {
        Matrix_pointer mp(-1,i,j,matrix);
        if(matrix == Simple_alignment::x_mat)
        {
            mp.fwd_score = (*xgap)[i][j].fwd_score;
            mp.bwd_score = (*xgap)[i][j].bwd_score;
            mp.full_score = (*xgap)[i][j].full_score;
        }
        else
        {
            mp.fwd_score = (*ygap)[i][j].fwd_score;
            mp.bwd_score = (*ygap)[i][j].bwd_score;
            mp.full_score = (*ygap)[i][j].full_score;
        }
        Path_pointer pp( mp, false, branch_length,1 );
        path->insert(path->begin(),pp);
    }

    void insert_gap_vector_path_pointer(vector<Path_pointer> *path, int i, int j, int matrix,float branch_length,int k)
    {
//        cout<<"si "<<i<<" "<<j<<" "<<matrix<<endl;
        Matrix_pointer mp(-1,i,j,matrix);
        if(matrix == Simple_alignment::x_mat)
        {
            mp.fwd_score = (*xvectp)[k].fwd_score;
            mp.bwd_score = (*xvectp)[k].bwd_score;
            mp.full_score = (*xvectp)[k].full_score;
        }
        else
        {
            mp.fwd_score = (*yvectp)[k].fwd_score;
            mp.bwd_score = (*yvectp)[k].bwd_score;
            mp.full_score = (*yvectp)[k].full_score;
        }
        Path_pointer pp( mp, false, branch_length,1 );
        path->insert(path->begin(),pp);
    }

    void insert_preexisting_gap(vector<Path_pointer> *path,int *i, int *j, int x_ind, int y_ind)
    {
        bool add_one_more = false;
        while(x_ind<*i)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0)
                left->get_edges()->at(ind).is_used(true);

            if(Settings::noise>6) /*DEBUG*/
                cout<<"x skip ("<<*i<<","<<*j<<") ["<<x_ind<<","<<y_ind<<"] "<<endl;

            this->insert_gap_path_pointer(path,*i-1,*j,Simple_alignment::x_mat,left_branch_length);
            --*i;

            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0)
                left->get_edges()->at(ind).is_used(true);
        }

        add_one_more = false;

        while(y_ind<*j)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0)
                right->get_edges()->at(ind).is_used(true);

            if(Settings::noise>6) /*DEBUG*/
                cout<<"y skip ("<<*i<<","<<*j<<") ["<<x_ind<<","<<y_ind<<"] "<<endl;

            this->insert_gap_path_pointer(path,*i,*j-1,Simple_alignment::y_mat,right_branch_length);
            --*j;
            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0)
                right->get_edges()->at(ind).is_used(true);
        }
    }

    void insert_preexisting_vector_gap(vector<Path_pointer> *path,int *i, int *j, int x_ind, int y_ind, int *k)
    {
        bool add_one_more = false;
        while(x_ind<*i)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0)
                left->get_edges()->at(ind).is_used(true);

            if(Settings::noise>6) /*DEBUG*/
                cout<<"x skip ("<<*i<<","<<*j<<") ["<<x_ind<<","<<y_ind<<"] "<<endl;

            this->insert_gap_vector_path_pointer(path,*i-1,*j,Simple_alignment::x_mat,left_branch_length,*k);
            --*i;
            --*k;

            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*i,*i+1);
            int ind = left->get_bwd_edge_index_at_site(*i+1,&edge);
            if(ind>=0)
                left->get_edges()->at(ind).is_used(true);
        }

        add_one_more = false;

        while(y_ind<*j)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0)
                right->get_edges()->at(ind).is_used(true);

            if(Settings::noise>6) /*DEBUG*/
                cout<<"y skip ("<<*i<<","<<*j<<") ["<<x_ind<<","<<y_ind<<"] "<<endl;

            this->insert_gap_vector_path_pointer(path,*i,*j-1,Simple_alignment::y_mat,right_branch_length,*k);
            --*j;
            --*k;

            add_one_more = true;
        }

        if(add_one_more)
        {
            Edge edge(*j,*j+1);
            int ind = right->get_bwd_edge_index_at_site(*j+1,&edge);
            if(ind>=0)
                right->get_edges()->at(ind).is_used(true);
        }
    }

    void insert_new_path_pointer(vector<Path_pointer> *path,int *i, int *j,Path_pointer pp)
    {
        if(Settings::noise>6) /*DEBUG*/
            cout<<pp.mp.matrix<<" ("<<*i<<","<<*j<<") ["<<pp.mp.x_ind<<","<<pp.mp.y_ind<<"] "<<pp.mp.score<<endl;;

        if( *i>0 || *j>0 )
            path->insert(path->begin(),pp);
    }


    /********************************************/

    void merge_sampled_sequence(Sequence *ancestral_sequence, Sequence *sampled_sequence);

    /********************************************/

    double get_log_edge_weight(Edge *edge) { return (*this.*log_edge_weight)(edge); }
    double (ppa::Simple_alignment::*log_edge_weight)(Edge *edge);

    double get_edge_weight(Edge *edge) { return (*this.*edge_weight)(edge); }
    double (ppa::Simple_alignment::*edge_weight)(Edge *edge);

    double edge_posterior_weight(Edge *edge) { return edge->get_posterior_weight(); }
    double edge_log_posterior_weight(Edge *edge) { return edge->get_log_posterior_weight(); }

    double edge_equal_weight(Edge *) const { return 1.0; }
    double edge_log_equal_weight(Edge *) const { return 0.0; }

    /********************************************/

    double get_transformed_edge_weight(double w) { return (*this.*transform_edge_weight)(w); }
    double (ppa::Simple_alignment::*transform_edge_weight)(double w);

    double square_root_edge_weight(double w) { return sqrt( w );}
    double cube_root_edge_weight(double w) { return exp((1.0/3.0)*(log(w)));}
    double plain_edge_weight(double w) { return w;}

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

//        cout<<"x_length "<<x_length<<"; y_length "<<y_length<<"; x_read1_length "<<x_read1_length<<"; y_read1_length "<<y_read1_length<<endl;
    }

    static int plot_number;
    void plot_posterior_probabilities_up();
    void plot_posterior_probabilities_down();

    void debug_print_matrices(int noise_level)
    {
        if(Settings::noise>noise_level)
            print_matrices();
    }

    void print_matrices();
    void print_sequences(vector<Site> *sites);
    string print_pairwise_alignment(vector<Site> *sites);

    void debug_print_input_sequences(int noise_level)
    {
        if(Settings::noise>noise_level)
            print_input_sequences();
    }

    void print_input_sequences()
    {

        cout<<"sequences:"<<endl<<" ";
        for(int i=1;i<left->sites_length()-1;i++)
            cout<<full_char_alphabet.at(left->get_site_at(i)->get_state());
        cout<<endl<<" ";
        for(int i=1;i<right->sites_length()-1;i++)
            cout<<full_char_alphabet.at(right->get_site_at(i)->get_state());
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

            cout<<i<<" ";
            if(tsite->real_site)
                cout<<" r: ";
            else
                cout<<" s: ";

            cout<<tsite->mp.x_ind<<" "<<tsite->mp.y_ind<<" "<<tsite->mp.matrix<<": "<<tsite->mp.score<<" "<<fixed<<log(tsite->mp.full_score)<<" "<<fixed<<log(tsite->mp.bwd_score)<<endl;
        }
    }

public:
    Simple_alignment();

    void align(Sequence *left_sequence,Sequence *right_sequence,Evol_model *model,
               float left_branch_length=0,float right_branch_length=0,bool is_reads_sequence=false);

    void read_alignment(Sequence *left_sequence,Sequence *right_sequence,Evol_model *model,
                        float left_branch_length=0,float right_branch_length=0);


    void make_alignment_path(vector<Matrix_pointer> *simple_path);
    void make_alignment_path2(vector<Matrix_pointer> *simple_path);

    Sequence* get_simple_sequence() { return ancestral_sequence; }

    vector<Path_pointer> get_alignment_path() { return path; }
};

}
#endif // SIMPLE_ALIGNMENT_H
