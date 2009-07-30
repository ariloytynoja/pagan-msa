#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <string>
#include <cmath>

using namespace std;

namespace ppa{

struct Edge
{
private:
    int index;

    int start_site_index;
    int end_site_index;

    float posterior_weight;
    float log_posterior_weight;

    int next_fwd_edge_index;
    int next_bwd_edge_index;

    bool used_in_alignment;

    int branch_count_since_last_used;
    float branch_distance_since_last_used;
    int branch_count_as_skipped_edge;


public:
    Edge(int s, int e): index(-1), start_site_index(s), end_site_index(e),
                               posterior_weight(1.0), log_posterior_weight(0), next_fwd_edge_index(-1), next_bwd_edge_index(-1)
                               , used_in_alignment(false)
                               ,branch_count_since_last_used(0),branch_distance_since_last_used(0)
                               ,branch_count_as_skipped_edge(0){}

    Edge(int s, int e, float w): index(-1), start_site_index(s), end_site_index(e),
                                        posterior_weight(w), log_posterior_weight(log(w)), next_fwd_edge_index(-1), next_bwd_edge_index(-1)
                                        , used_in_alignment(false)
                                        ,branch_count_since_last_used(0),branch_distance_since_last_used(0)
                                        ,branch_count_as_skipped_edge(0){}

    bool is_used() { return used_in_alignment; }
    void is_used(bool set) { used_in_alignment = set; }

    int get_branch_count_since_last_used() { return branch_count_since_last_used; }
    float get_branch_distance_since_last_used() { return branch_distance_since_last_used; }
    int get_branch_count_as_skipped_edge() { return branch_count_as_skipped_edge; }

    void set_branch_count_since_last_used(int s) { branch_count_since_last_used = s; }
    void set_branch_distance_since_last_used(float s) { branch_distance_since_last_used = s; }
    void set_branch_count_as_skipped_edge(int s) { branch_count_as_skipped_edge = s; }

    void increase_branch_count_as_skipped_edge() { ++branch_count_as_skipped_edge; }

    int get_start_site_index() { return start_site_index; }
    int get_end_site_index() { return end_site_index; }

    void set_start_site_index(int s) { start_site_index = s; }
    void set_end_site_index(int s) { end_site_index = s; }

    int get_index() { return index; }
    void set_index(int i) { index = i; }

    int get_next_fwd_edge_index() { return next_fwd_edge_index; }
    void set_next_fwd_edge_index(int i) { next_fwd_edge_index = i; }

    int get_next_bwd_edge_index() { return next_bwd_edge_index; }
    void set_next_bwd_edge_index(int i) { next_bwd_edge_index = i; }

    double get_posterior_weight() { return posterior_weight; }
    double get_log_posterior_weight() { return log_posterior_weight; }

    void set_weight(float w) { posterior_weight = w; log_posterior_weight = log(w); }

    bool operator==(const Edge& b)
    {
        return start_site_index==b.start_site_index && end_site_index==b.end_site_index;
    }

    bool operator!=(const Edge& b)
    {
        return !(start_site_index==b.start_site_index && end_site_index==b.end_site_index);
    }

    bool operator<(const Edge& b)
    {
        return start_site_index<b.start_site_index && end_site_index<b.end_site_index;
    }

    bool operator>(const Edge& b)
    {
        return start_site_index>b.start_site_index && end_site_index>b.end_site_index;
    }

};

/**************************************/

struct Site_children
{
    int left_index;
    int right_index;

public:

    bool operator==(const Site_children& b)
    {
        return left_index==b.left_index && right_index==b.right_index;
    }

    bool operator!=(const Site_children& b)
    {
        return !(left_index==b.left_index && right_index==b.right_index);
    }

    bool operator<(const Site_children& b)
    {
        return left_index<b.left_index && right_index<b.right_index;
    }

    bool operator>(const Site_children& b)
    {
        return left_index>b.left_index && right_index>b.right_index;
    }
};

/**************************************/

struct Site
{
    int index;
    Site_children children;

    int character_state;

    int site_type;
    enum Site_type {start_site,real_site,stop_site};

    int path_state;
    enum Path_state {ends_site,terminal,matched,xgapped,ygapped,xskipped,yskipped};

    vector<Edge> *edges;

    int first_fwd_edge_index;
    int current_fwd_edge_index;

    int first_bwd_edge_index;
    int current_bwd_edge_index;

    float posterior_support;

    int branch_count_since_last_used;
    float branch_distance_since_last_used;


public:
    Site(vector<Edge> *e,int type=Site::real_site,int p_state=Site::terminal):index(-1),character_state(-1),
            site_type(type),path_state(p_state),edges(e),first_fwd_edge_index(-1),current_fwd_edge_index(-1),
            first_bwd_edge_index(-1),current_bwd_edge_index(-1),posterior_support(1),
            branch_count_since_last_used(0),branch_distance_since_last_used(0) {}

    int get_branch_count_since_last_used() { return branch_count_since_last_used; }
    float get_branch_distance_since_last_used() { return branch_distance_since_last_used; }

    void set_branch_count_since_last_used(int s) { branch_count_since_last_used = s; }
    void set_branch_distance_since_last_used(float s) { branch_distance_since_last_used = s; }

    void set_edge_vector(vector<Edge> *e) { edges = e; }

    void set_state(int c) { character_state = c; }
    int  get_state() { return character_state; }

    void set_path_state(int c) { path_state = c; }
    int  get_path_state() { return path_state; }

    void set_index(int i) { index = i; }
    int  get_index() { return index; }

    int get_site_type() { return site_type; }

    /**************************************/

    void set_empty_children()
    {
        children.left_index = -1; children.right_index = -1;
    }

    void set_children(int li, int ri)
    {
        children.left_index = li; children.right_index = ri;
    }

    Site_children * get_children() { return &children; }

    /**************************************/

    void set_first_fwd_edge_index(int i)
    {
        first_fwd_edge_index = current_fwd_edge_index = i;
    }

    void set_first_bwd_edge_index(int i)
    {
        first_bwd_edge_index = current_bwd_edge_index = i;
    }

    /**************************************/

    void add_new_fwd_edge_index(int i)
    {
        if(first_fwd_edge_index<0)
        {
            this->set_first_fwd_edge_index(i);
            return;
        }
        int prev_fwd_edge_index = current_fwd_edge_index;
        current_fwd_edge_index = i;
        edges->at(prev_fwd_edge_index).set_next_fwd_edge_index(current_fwd_edge_index);
    }

    void add_new_bwd_edge_index(int i) {
        if(first_bwd_edge_index<0)
        {
            this->set_first_bwd_edge_index(i);
            return;
        }
        int prev_bwd_edge_index = current_bwd_edge_index;
        current_bwd_edge_index = i;
        edges->at(prev_bwd_edge_index).set_next_bwd_edge_index(current_bwd_edge_index);
    }

    /**************************************/

    bool has_fwd_edge()
    {
        return first_fwd_edge_index >= 0;
    }

    Edge * get_first_fwd_edge()
    {
        current_fwd_edge_index = first_fwd_edge_index;
        return &edges->at(current_fwd_edge_index);
    }

    bool has_next_fwd_edge()
    {
        return edges->at(current_fwd_edge_index).get_next_fwd_edge_index() >= 0;
    }

    Edge * get_next_fwd_edge()
    {
        if(edges->at(current_fwd_edge_index).get_next_fwd_edge_index() < 0)
            return 0;
        current_fwd_edge_index = edges->at(current_fwd_edge_index).get_next_fwd_edge_index();

        return &edges->at(current_fwd_edge_index);
    }

    /**************************************/

    bool has_bwd_edge()
    {
        return first_bwd_edge_index >= 0;
    }

    Edge * get_first_bwd_edge()
    {
        current_bwd_edge_index = first_bwd_edge_index;
        return &edges->at(current_bwd_edge_index);
    }

    bool has_next_bwd_edge()
    {
        return edges->at(current_bwd_edge_index).get_next_bwd_edge_index() >= 0;
    }

    Edge * get_next_bwd_edge()
    {
        if(edges->at(current_bwd_edge_index).get_next_bwd_edge_index() < 0)
            return 0;
        current_bwd_edge_index = edges->at(current_bwd_edge_index).get_next_bwd_edge_index();
        return &edges->at(current_bwd_edge_index);
    }

    bool contains_bwd_edge(Edge *copy)
    {
        if( this->has_bwd_edge() )
        {
            Edge *edge = this->get_first_bwd_edge();
            if(*copy == *edge)
                return true;
            while(this->has_next_bwd_edge())
            {
                edge = this->get_next_bwd_edge();
                if(*copy == *edge)
                    return true;
            }
        }
        return false;
    }

    void delete_bwd_edge(int edge_ind)
    {
        if(this->has_bwd_edge())
        {
            Edge *edge = this->get_first_bwd_edge();

            if(edge->get_index() == edge_ind)
            {
                if(this->has_next_bwd_edge())
                {
                    edge = this->get_next_bwd_edge();
                    int next_ind = edge->get_index();
                    current_bwd_edge_index = first_bwd_edge_index = next_ind;
                }
                else
                {
                    current_bwd_edge_index = first_bwd_edge_index = -1;
                }
                return;
            }

            int prev_ind = edge->get_index();

            while(this->has_next_bwd_edge())
            {
                edge = this->get_next_bwd_edge();
                if(edge->get_index() == edge_ind)
                {
                    if(this->has_next_bwd_edge())
                    {
                        edge = this->get_next_bwd_edge();
                        int next_ind = edge->get_index();
                        edges->at(prev_ind).set_next_bwd_edge_index(next_ind);
                    }
                    else
                    {
                        edges->at(prev_ind).set_next_bwd_edge_index(-1);
                    }
                }
            }
        }
    }

    void delete_fwd_edge(int edge_ind)
    {
        if(this->has_fwd_edge())
        {
            Edge *edge = this->get_first_fwd_edge();

            if(edge->get_index() == edge_ind)
            {
                if(this->has_next_fwd_edge())
                {
                    edge = this->get_next_fwd_edge();
                    int next_ind = edge->get_index();
                    current_fwd_edge_index = first_fwd_edge_index = next_ind;
                }
                else
                {
                    current_fwd_edge_index = first_fwd_edge_index = -1;
                }
                return;
            }

            int prev_ind = edge->get_index();

            while(this->has_next_fwd_edge())
            {
                edge = this->get_next_fwd_edge();
                if(edge->get_index() == edge_ind)
                {
                    if(this->has_next_fwd_edge())
                    {
                        edge = this->get_next_fwd_edge();
                        int next_ind = edge->get_index();
                        edges->at(prev_ind).set_next_fwd_edge_index(next_ind);
                    }
                    else
                    {
                        edges->at(prev_ind).set_next_fwd_edge_index(-1);
                    }
                }
            }
        }
    }

};

/**************************************/

class Sequence
{
    int curr_site_index;
    int prev_site_index;
    int curr_edge_index;

    vector<Site> sites;
    vector<Edge> edges;
    string full_dna_alphabet;

public:

    void initialise_indeces() {
        curr_site_index = prev_site_index = curr_edge_index = 0;
    }

    int get_new_site_index() {curr_site_index++; return curr_site_index; }
    int get_new_edge_index() {curr_edge_index++; return curr_edge_index; }

    int get_current_site_index() { return curr_site_index; }
    int get_previous_site_index() { return prev_site_index; }
    int get_current_edge_index() { return curr_edge_index; }

    void push_back_site(Site site)
    {
        site.set_index(sites.size());
        sites.push_back(site);

        prev_site_index = sites.size()-2;
        curr_site_index = sites.size()-1;
    }

    void push_back_edge(Edge edge)
    {
        edge.set_index(edges.size());
        edges.push_back(edge);

        curr_edge_index = edges.size()-1;
    }

    bool contains_edge(Edge edge)
    {
        vector<Edge>::iterator it = edges.begin();
        for(;it!=edges.end();it++)
        {
            if(*it==edge)
                return true;
        }
        return false;
    }

    Site *get_previous_site()
    {
        return &sites.at(prev_site_index);
    }

    Site *get_current_site()
    {
        return &sites.at(curr_site_index);
    }

    int get_bwd_edge_index_at_site(int site,Edge *copy)
    {
        if(!sites.at(site).has_bwd_edge())
            return -1;

        Edge *edge = sites.at(site).get_first_bwd_edge();
        if(*edge == *copy)
            return edge->get_index();

        while(sites.at(site).has_next_bwd_edge())
        {
            edge = sites.at(site).get_next_bwd_edge();
            if(*edge == *copy)
                return edge->get_index();
        }
        return -1;
    }

    bool contains_this_bwd_edge_at_site(int site, Edge *copy)
    {
        if(this->get_bwd_edge_index_at_site(site,copy)<0)
            return false;

        return true;
    }

    int get_fwd_edge_index_at_site(int site,Edge *copy)
    {
        if(!sites.at(site).has_fwd_edge())
            return -1;

        Edge *edge = sites.at(site).get_first_fwd_edge();
        if(*edge == *copy)
            return edge->get_index();

        while(sites.at(site).has_next_fwd_edge())
        {
            edge = sites.at(site).get_next_fwd_edge();
            if(*edge == *copy)
                return edge->get_index();
        }
        return -1;
    }

    bool contains_this_fwd_edge_at_site(int site, Edge *copy)
    {
        if(this->get_fwd_edge_index_at_site(site,copy)<0)
            return false;

        return true;
    }

    void print_sequence(vector<Site> *sites);
    void print_sequence() { this->print_sequence(this->get_sites()); }


    Sequence(const string &seq_string,const string &alphabet);
    Sequence(const vector<Site>* s, const vector<Edge>* e, const string& alphabet);
    Sequence(const int length,const string& alphabet);

    vector<Site> *get_sites() { return &sites; }
    vector<Edge> *get_edges() { return &edges; }

    int sites_length() { return sites.size(); }
    int edges_length() { return edges.size(); }

    Site *get_site_at(int i) { return &sites.at(i); }
    Edge *get_first_bwd_edge_at(int i) { return sites.at(i).get_first_bwd_edge(); }

    string get_full_alphabet() { return full_dna_alphabet; }

    void delete_all_bwd_edges_at_site(int index)
    {
        Site *site = this->get_site_at(index);

        if(site->has_bwd_edge())
        {
            Edge *edge = site->get_first_bwd_edge();
            this->get_site_at(edge->get_start_site_index())->delete_fwd_edge(edge->get_index());

            while(site->has_next_bwd_edge())
            {
                edge = site->get_next_bwd_edge();
                this->get_site_at(edge->get_start_site_index())->delete_fwd_edge(edge->get_index());
            }
        }
        site->set_first_bwd_edge_index(-1);
    }

    void delete_all_fwd_edges_at_site(int index)
    {
        Site *site = this->get_site_at(index);

        if(site->has_fwd_edge())
        {
            Edge *edge = site->get_first_fwd_edge();
            this->get_site_at(edge->get_end_site_index())->delete_bwd_edge(edge->get_index());

            while(site->has_next_fwd_edge())
            {
                edge = site->get_next_fwd_edge();
                this->get_site_at(edge->get_end_site_index())->delete_bwd_edge(edge->get_index());
            }
        }
        site->set_first_fwd_edge_index(-1);
    }


};
}

#endif // SEQUENCE_H
