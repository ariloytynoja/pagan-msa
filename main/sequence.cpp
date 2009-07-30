#include <iostream>
#include "utils/settings.h"
#include "main/sequence.h"
#include <iomanip>

using namespace std;
using namespace ppa;

Sequence::Sequence(const string& seq_string,const string& alphabet)
{

    this->initialise_indeces();

    sites.reserve(seq_string.size()+2);
    edges.reserve(seq_string.size()+3);

    full_dna_alphabet = alphabet;

    Site first_site( &edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_empty_children();
    this->push_back_site(first_site);

//    Site site( &edges );
//    site.set_state( full_dna_alphabet.find( seq_string.at(0) ) );
//    site.set_empty_children();
//    this->push_back_site(site);

    Edge first_edge( -1,this->get_current_site_index() );
    this->push_back_edge(first_edge);

//    this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );

    for(unsigned int i=0;i<seq_string.size();i++)
    {
//cout<<i<<endl;
        Site site( &edges );
        site.set_state( full_dna_alphabet.find( seq_string.at(i) ) );
        site.set_empty_children();
        this->push_back_site(site);

        Edge edge( this->get_previous_site_index(),this->get_current_site_index() );
        this->push_back_edge(edge);

        this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
        this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );
    }

    Site last_site( &edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_empty_children();
    this->push_back_site(last_site);

    Edge last_edge( this->get_previous_site_index(),this->get_current_site_index() );
    this->push_back_edge(last_edge);

    this->get_previous_site()->set_first_fwd_edge_index( this->get_current_edge_index() );
    this->get_current_site()->set_first_bwd_edge_index( this->get_current_edge_index() );


    if(Settings::noise>5)
    {
        this->print_sequence(&sites);
    }
}

Sequence::Sequence(const int length,const string& alphabet)
{

    this->initialise_indeces();

    sites.reserve(length+2);
    edges.reserve(length+3);

    full_dna_alphabet = alphabet;
}


Sequence::Sequence(const vector<Site>* s, const vector<Edge>* e, const string& alphabet)
{
    sites = *s;
    edges = *e;
    full_dna_alphabet = alphabet;
    for(unsigned int i=0;i<sites.size();i++)
    {
        sites.at(i).set_edge_vector(&edges);
    }

    if(Settings::noise>5)
    {
        this->print_sequence(&sites);
    }
}

void Sequence::print_sequence(vector<Site> *sites)
{
    cout<<endl;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);
//cout<<i<<" "<<tsite->get_state()<<endl;

        if(tsite->get_site_type()==Site::real_site)
            cout<<tsite->get_index()<<" "<<full_dna_alphabet.at(tsite->get_state());
        else
            cout<<tsite->get_index()<<" +";

//        if(tsite->has_fwd_edge())
//        {
//            Edge *tedge = tsite->get_first_fwd_edge();
//            cout<<" F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
//            while(tsite->has_next_fwd_edge())
//            {
//                tedge = tsite->get_next_fwd_edge();
//                cout<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
//            }
//        }
//        if(tsite->has_bwd_edge())
//        {
//            Edge *tedge = tsite->get_first_bwd_edge();
//            cout<<"; B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
//            while(tsite->has_next_bwd_edge())
//            {
//                tedge = tsite->get_next_bwd_edge();
//                cout<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
//            }
//        }

        cout<<"\t";
        cout << setprecision (2);

        if(true)
        {
            if(tsite->has_fwd_edge())
            {
                Edge *tedge = tsite->get_first_fwd_edge();
                cout<<" F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                        <<" "<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";
                while(tsite->has_next_fwd_edge())
                {
                    tedge = tsite->get_next_fwd_edge();
                    cout<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                        <<" "<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";
                }
            }
            cout<<"; \t";
        }
        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            cout<<"B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                    <<" "<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";
            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                cout<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index()<<" ["<<tedge->get_log_posterior_weight()
                    <<" "<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";
            }
        }
        cout << setprecision (4);
        
        cout<<"\n";

    }
}
