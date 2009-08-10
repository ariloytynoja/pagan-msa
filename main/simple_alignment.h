#ifndef SIMPLE_ALIGNMENT_H
#define SIMPLE_ALIGNMENT_H

#include "utils/dna_model.h"
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
    double bwd_score;
    int x_ind;
    int y_ind;
    int x_edge_ind;
    int y_edge_ind;
    int matrix;

    Matrix_pointer() : score(-HUGE_VAL), full_score(0), bwd_score(0), x_ind(-1), y_ind(-1),
                        x_edge_ind(-1), y_edge_ind(-1), matrix(-1) {}
    Matrix_pointer(double s,int x, int y, int m) : score(s), full_score(0), bwd_score(0),
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

    typedef boost::multi_array<Matrix_pointer, 2> align_array;
    typedef align_array::index index;

    align_array *match;
    align_array *xgap;
    align_array *ygap;

    typedef boost::multi_array_types::index_range range;
    align_array::index_gen indices;
    typedef align_array::array_view<1>::type align_slice;

    Sequence *left;
    Sequence *right;
    Sequence *ancestral_sequence;

    float left_branch_length;
    float right_branch_length;

    Dna_model *model;
    string full_dna_alphabet;

    vector<Path_pointer> path;

    // parameters that may be set by program arguments
    //
    float ins_del_ratio;   // (not used by now)
    float del_ins_ratio;   // (not used by now)

    float max_allowed_skip_distance;
    int   max_allowed_skip_branches;
    int   max_allowed_match_skip_branches;
    float branch_skip_weight;
    float branch_skip_probability;
    bool  weighted_branch_skip_penalty;

    bool compute_full_score;
    double bwd_full_probability; // for control

    void initialise_array_corner();
    void initialise_array_corner_bwd();

    void find_best_transition(int i,int j);
    void compute_bwd_full_score(int i,int j);

    void backtrack_new_path(vector<Path_pointer> *path,Path_pointer pp);
    void build_ancestral_sequence(vector<Path_pointer> *path);

    void create_ancestral_sequence(vector<Path_pointer> *path);
    void create_ancestral_edges();
    void check_skipped_boundaries();

    void delete_edge_range(int edge_ind,int skip_start_site);

    void transfer_child_edge(Edge *child, vector<int> *child_index, float branch_length,
                             bool adjust_posterior_weight = true, float branch_weight = 1.0);
    void transfer_child_edge(Edge edge, Edge *child, float branch_length,
                             bool adjust_posterior_weight = true, float branch_weight = 1.0);

    /*********************************/

    void iterate_bwd_edges_for_gap(Site * site,align_slice *x_slice,align_slice *y_slice,align_slice *m_slice,
                                   Matrix_pointer *max,bool is_x_matrix);
    void iterate_bwd_edges_for_match(Site * left_site,Site * right_site,Matrix_pointer *max);
    void iterate_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max);

    void iterate_fwd_edges_for_gap(Site * site,align_slice *g_slice,
                                   Matrix_pointer *max_s,Matrix_pointer *max_d,Matrix_pointer *max_m);
    void iterate_fwd_edges_for_match(Site * left_site,Site * right_site,
                                     Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);


    /*********************************/

    void score_m_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_x_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);
    void score_y_match(Edge * left_edge,Edge * right_edge,double log_match,Matrix_pointer *max, double match = 0);

    void score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_double(Edge *edge,align_slice *w_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_open(Edge *edge,align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix);
    void score_gap_close(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix);

    /*********************************/

    void score_gap_ext_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max);
    void score_gap_double_bwd(Edge *edge,align_slice *w_slice,Matrix_pointer *max);
    void score_gap_open_bwd(Edge *edge,align_slice *m_slice,Matrix_pointer *max);

    void score_match_bwd(Edge * left_edge,Edge * right_edge,
                                       Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m);

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
            if ((double)rand()/(double)RAND_MAX>0.5)
                return true;
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

    void insert_new_path_pointer(vector<Path_pointer> *path,int *i, int *j,Path_pointer pp)
    {
        if(Settings::noise>6) /*DEBUG*/
            cout<<pp.mp.matrix<<" ("<<*i<<","<<*j<<") ["<<pp.mp.x_ind<<","<<pp.mp.y_ind<<"] "<<pp.mp.score<<endl;;

        if( *i>0 || *j>0 )
            path->insert(path->begin(),pp);
    }


    /********************************************/


    double get_log_edge_weight(Edge *edge) { return (*this.*log_edge_weight)(edge); }
    double (ppa::Simple_alignment::*log_edge_weight)(Edge *edge);

    double get_edge_weight(Edge *edge) { return (*this.*edge_weight)(edge); }
    double (ppa::Simple_alignment::*edge_weight)(Edge *edge);

    double edge_posterior_weight(Edge *edge) { return edge->get_posterior_weight(); }
    double edge_log_posterior_weight(Edge *edge) { return edge->get_log_posterior_weight(); }

    double edge_equal_weight(Edge *) const { return 1.0; }
    double edge_log_equal_weight(Edge *) const { return 0.0; }

    void set_basic_settings()
    {
        del_ins_ratio = 1;//model->del_prob/model->ins_prob;
        ins_del_ratio = 1;//model->ins_prob/model->del_prob;

        weighted_branch_skip_penalty = false; // by default, use penalty *per node*

        max_allowed_skip_distance = 0.5;
        if(Settings_handle::st.is("branch-length-confirm-insertion"))
            max_allowed_skip_distance = Settings_handle::st.get("branch-length-confirm-insertion").as<float>();

        max_allowed_skip_branches = 5;
        if(Settings_handle::st.is("any-skips-confirm-insertion"))
            max_allowed_skip_branches = Settings_handle::st.get("any-skips-confirm-insertion").as<int>();

        max_allowed_match_skip_branches = 2;
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
            branch_skip_weight = Settings_handle::st.get("branch-skip-penalty-per-branch").as<float>();
            weighted_branch_skip_penalty = false;
        }

        compute_full_score = false;
        if(Settings_handle::st.is("full-probability"))
        {
            compute_full_score = true;
        }

        cout << noshowpos;
    }

    static int plot_number;
    void plot_posterior_probabilities(Matrix_pointer max_end);

    void print_matrices();
    void print_sequences(vector<Site> *sites);

    void print_input_sequences()
    {

        cout<<"sequences:"<<endl<<" ";
        for(int i=1;i<left->sites_length()-1;i++)
            cout<<full_dna_alphabet.at(left->get_site_at(i)->get_state());
        cout<<endl<<" ";
        for(int i=1;i<right->sites_length()-1;i++)
            cout<<full_dna_alphabet.at(right->get_site_at(i)->get_state());
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

    void align(Sequence *left_sequence,Sequence *right_sequence,Dna_model *model,float left_branch_length=0,float right_branch_length=0);


    Sequence* get_simple_sequence() { return ancestral_sequence; }

    vector<Path_pointer> get_alignment_path() { return path; }
};

}
#endif // SIMPLE_ALIGNMENT_H
