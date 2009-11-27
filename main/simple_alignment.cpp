#include "main/simple_alignment.h"
#include "main/sequence.h"
#include "utils/exceptions.h"
#include "utils/node.h"
#include <iomanip>
#include <fstream>

using namespace std;
using namespace ppa;

Simple_alignment::Simple_alignment()
{
}

void Simple_alignment::align(Sequence *left_sequence,Sequence *right_sequence,Dna_model *dna_model,float l_branch_length,float r_branch_length)
{

    left = left_sequence;
    right = right_sequence;
    model = dna_model;
    full_dna_alphabet = model->get_full_alphabet();

    left_branch_length = l_branch_length;
    right_branch_length = r_branch_length;


    this->debug_print_input_sequences(2);


    // set the basic parameters (copy from Settings)
    //
    this->set_basic_settings();

    // Set the edge weighting scheme
    //
    log_edge_weight = &ppa::Simple_alignment::edge_log_posterior_weight;
    edge_weight = &ppa::Simple_alignment::edge_posterior_weight;


    int left_length = left_sequence->sites_length();
    int right_length = right_sequence->sites_length();

    this->debug_msg("Simple_alignment: lengths: "+this->itos(left_length)+" "+this->itos(right_length),1);

    align_array M(boost::extents[left_length-1][right_length-1]);
    align_array X(boost::extents[left_length-1][right_length-1]);
    align_array Y(boost::extents[left_length-1][right_length-1]);

    this->debug_msg("Simple_alignment: matrix created",1);

    match = &M;
    xgap = &X;
    ygap = &Y;

    this->initialise_array_corner();

    int j_max = match->shape()[1];
    int i_max = match->shape()[0];

    Site * left_site = left->get_site_at(i_max);
    Site * right_site = right->get_site_at(j_max);



    // Dynamic programming loop
    // Note that fwd full probability is integrated in
    // the main dp loop but bwd comutation is done separately
    //
    for(int j=0;j<j_max;j++)
    {
        for(int i=0;i<i_max;i++)
        {
            this->compute_fwd_viterbi_path(i,j);
        }
    }
    this->debug_msg("Simple_alignment: matrix filled",1);



    // Find the incoming edge in the end corner; also, get the full_fwd probability
    //
    Matrix_pointer max_end;
    this->iterate_bwd_edges_for_end_corner(left_site,right_site,&max_end);
    this->debug_msg("Simple_alignment: corner found",1);



    // If needed, compute the bwd posterior probabilities
    //
    if(compute_full_score)
    {
        this->initialise_array_corner_bwd();

        for(int j=j_max-1;j>=0;j--)
        {
            for(int i=i_max-1;i>=0;i--)
            {
                this->compute_bwd_full_score(i,j);
            }
        }

        // Check that full probability computation worked
        //
        Matrix_pointer max_start = (*match)[0][0];

        this->debug_msg(" bwd full probability: "+this->ftos(log(max_start.bwd_score))
                        +" ["+this->ftos(max_start.bwd_score)+"]",1);

        double ratio = max_end.full_score/max_start.bwd_score;

        if(ratio<0.99 || ratio>1.01)
            this->debug_msg("Problem in computation? fwd: "+this->ftos(max_end.full_score)
                            +", bwd: "+this->ftos(max_start.bwd_score),0);

        if(Settings::noise>5)
            this->print_matrices();

    }



    // It is more convenient to build the sequence forward;
    // thus, backtrace the path into a vector ('path') and use that in the next step
    //
    // Backtrack the Viterbi path
    if( ! Settings_handle::st.is("sample-path") )
    {
        Path_pointer pp(max_end,true);
        this->backtrack_new_path(&path,pp);

        this->debug_msg("Simple_alignment: path found",1);


        // Now build the sequence forward following the path saved in a vector;
        //
        ancestral_sequence = new Sequence(path.size(),model->get_full_alphabet());
        this->build_ancestral_sequence(ancestral_sequence,&path);

        this->debug_msg("Simple_alignment: sequence built",1);
    }
    // Sample a path from the posterior probabilities
    else
    {
        Matrix_pointer sample_end;
        this->iterate_bwd_edges_for_sampled_end_corner(left_site,right_site,&sample_end);

        Path_pointer sp(sample_end,true);
        vector<Path_pointer> sample_path;

        this->sample_new_path(&sample_path,sp);

        this->debug_msg("Simple_alignment: path sampled",1);


        // Now build the sequence forward following the path saved in a vector;
        //
        ancestral_sequence = new Sequence(sample_path.size(),model->get_full_alphabet());
        this->build_ancestral_sequence(ancestral_sequence,&sample_path);

        this->debug_msg("Simple_alignment: sequence sampled and built",1);
    }

    // Find the incoming edge in the end corner
    //
    if( compute_full_score && Settings_handle::st.is("sample-additional-paths") )
    {
        int iter = Settings_handle::st.get("sample-additional-paths").as<int>();
        for(int i=0;i<iter;i++)
        {
            Matrix_pointer sample_end;
            this->iterate_bwd_edges_for_sampled_end_corner(left_site,right_site,&sample_end);

            Path_pointer sp(sample_end,true);
            vector<Path_pointer> sample_path;

            this->sample_new_path(&sample_path,sp);

            this->debug_msg("Simple_alignment: additional path sampled",1);

            // Now build the sequence forward following the path saved in a vector;
            //
            Sequence *sampled_sequence = new Sequence(sample_path.size(),model->get_full_alphabet());
            this->build_ancestral_sequence(sampled_sequence,&sample_path);

            this->debug_msg("Simple_alignment: additional sequence sampled and built",1);

//            create a unique index of sites in sequences 'ancestral' and 'sampled'
//             - gaps need to contain info of a site in other child sequence that they follow
//            check if sites in 'sampled' are included in 'ancestral'
//             - add if not, add also in unique index
//            add edges that are missing
//            finally, reorder sites for the next alignment
//             - note that (some) pointers are *index-based*

            cout<<"\nancestral\n";
            ancestral_sequence->print_path();
            cout<<"\nsampled\n";
            sampled_sequence->print_path();
            cout<<"\nancestral\n";
            ancestral_sequence->print_sequence();
            cout<<"\nsampled\n";
            sampled_sequence->print_sequence();

            this->merge_sampled_sequence(ancestral_sequence, sampled_sequence);

            cout<<"\nmerged ancestral\n";
            ancestral_sequence->print_sequence();

            delete sampled_sequence;
        }
    }


    // Make the posterior probability plots
    //
    if(Settings_handle::st.is("mpost-posterior-plot-file") &&
       Settings_handle::st.is("full-probability") )
    {
        double full_score = max_end.full_score;

        for(int j=j_max-1;j>=0;j--)
        {
            for(int i=i_max-1;i>=0;i--)
            {
                this->compute_posterior_score(i,j,full_score);
            }
        }

        if(Settings_handle::st.is("plot-slope-up") )
            this->plot_posterior_probabilities_up();
        else
            this->plot_posterior_probabilities_down();
    }

}


/********************************************/

void Simple_alignment::merge_sampled_sequence(Sequence *ancestral_sequence, Sequence *sampled_sequence)
{
    ancestral_sequence->initialise_unique_index();
    sampled_sequence->initialise_unique_index();

    vector<int> sample_index_for_added;
    vector<int> sample_index_in_ancestral;

    vector<Unique_index> *index = sampled_sequence->get_unique_index();
    for(int i=0;i<(int) index->size();i++)
    {
        int anc_index = ancestral_sequence->unique_index_of_term( &index->at(i) );

        if( anc_index>=0 )
        {
            sample_index_in_ancestral.push_back(anc_index);
        }
        else
        {
            // Make a copy of the site in sampled sequence that is missing from
            // ancestral sequence and add it there; can't make a carbon copy as
            // indeces (neighbour sites, edges) aren't correct
            //
            Site *sample_site = sampled_sequence->get_site_at( index->at(i).site_index );

            Site site( ancestral_sequence->get_edges() );
            ancestral_sequence->copy_site_details(sample_site, &site);
            ancestral_sequence->push_back_site( site );

            ancestral_sequence->add_term_in_unique_index( index->at(i) );

            // Keep track of sites added so that edges can be added below
            //
            sample_index_in_ancestral.push_back( ancestral_sequence->get_current_site_index() );
            sample_index_for_added.push_back( sample_index_in_ancestral.size()-1 );
        }

    }

    // If sites were added, add also their edges
    //
    for(int i=0;i<(int) sample_index_for_added.size();i++)
    {
        Site *sample_site = sampled_sequence->get_site_at( sample_index_for_added.at( i ) );

        if( sample_site->has_bwd_edge() )
        {
            Edge *sample_edge = sample_site->get_first_bwd_edge();

            int edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
            int edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

            Edge edge( edge_start_site, edge_end_site );

            if(!ancestral_sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
            {

                ancestral_sequence->copy_edge_details( sample_edge, &edge );
                ancestral_sequence->push_back_edge(edge);

                ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
            }

            while( sample_site->has_next_bwd_edge() )
            {
                sample_edge = sample_site->get_next_bwd_edge();

                edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
                edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

                Edge edge( edge_start_site, edge_end_site );

                if(!ancestral_sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
                {

                    ancestral_sequence->copy_edge_details( sample_edge, &edge );
                    ancestral_sequence->push_back_edge(edge);

                    ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                    ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
                }
            }
        }

        // The same for fwd edges
        if( sample_site->has_fwd_edge() )
        {
            Edge *sample_edge = sample_site->get_first_fwd_edge();

            int edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
            int edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

            Edge edge( edge_start_site, edge_end_site );

            if(!ancestral_sequence->contains_this_fwd_edge_at_site(edge.get_start_site_index(),&edge))
            {

                ancestral_sequence->copy_edge_details( sample_edge, &edge );
                ancestral_sequence->push_back_edge(edge);

                ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
            }

            while( sample_site->has_next_fwd_edge() )
            {
                sample_edge = sample_site->get_next_fwd_edge();

                edge_start_site = sample_index_in_ancestral.at ( sample_edge->get_start_site_index() );
                edge_end_site = sample_index_in_ancestral.at ( sample_edge->get_end_site_index() );

                Edge edge( edge_start_site, edge_end_site );

                if(!ancestral_sequence->contains_this_fwd_edge_at_site(edge.get_start_site_index(),&edge))
                {

                    ancestral_sequence->copy_edge_details( sample_edge, &edge );
                    ancestral_sequence->push_back_edge(edge);

                    ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
                    ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
                }
            }
        }
    }

    if( sample_index_for_added.size()>0 )
    {
        // Now reorder the ancestral sequence

        ancestral_sequence->sort_sites_vector();

        ancestral_sequence->remap_edges_vector();

    }
}

/********************************************/

void Simple_alignment::initialise_array_corner()
{
    double small = -HUGE_VAL;

    ((*match)[0][0]).score = 0.0;
    ((*match)[0][0]).full_score = 1.0;
    ((*xgap)[0][0]).score = small;
    ((*ygap)[0][0]).score = small;


    /*MISSING FEATURE: allow extending existing gap (for fragments) */
}

/********************************************/

void Simple_alignment::initialise_array_corner_bwd()
{
    int i_max = match->shape()[0];
    int j_max = match->shape()[1];

    Site * left_site = left->get_site_at(i_max);
    Site * right_site = right->get_site_at(j_max);

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {
        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        double left_edge_wght = this->get_edge_weight(left_edge);
        int left_index = left_edge->get_start_site_index();

        double right_edge_wght = this->get_edge_weight(right_edge);
        int right_index = right_edge->get_start_site_index();

        (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                right_edge_wght = this->get_edge_weight(right_edge);
                right_index = right_edge->get_start_site_index();

                (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;
            }
        }

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;

            while(left_site->has_next_bwd_edge())
            {
                left_edge = left_site->get_next_bwd_edge();

                left_edge_wght = this->get_edge_weight(left_edge);
                left_index = left_edge->get_start_site_index();

                (*match)[left_index][right_index].bwd_score = model->non_gap() * left_edge_wght * right_edge_wght;
            }
        }
    }

    if( left_site->has_bwd_edge() )
    {
        Edge * left_edge = left_site->get_first_bwd_edge();

        double left_edge_wght = this->get_edge_weight(left_edge);
        int left_index = left_edge->get_start_site_index();

        (*xgap)[left_index][j_max-1].bwd_score = model->gap_close() * left_edge_wght;

        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();

            left_edge_wght = this->get_edge_weight(left_edge);
            left_index = left_edge->get_start_site_index();

            (*xgap)[left_index][j_max-1].bwd_score = model->gap_close() * left_edge_wght;
        }
    }

    if(right_site->has_bwd_edge())
    {
        Edge * right_edge = right_site->get_first_bwd_edge();

        double right_edge_wght = this->get_edge_weight(right_edge);
        int right_index = right_edge->get_start_site_index();

        (*ygap)[i_max-1][right_index].bwd_score = model->gap_close() * right_edge_wght;

        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();

            right_edge_wght = this->get_edge_weight(right_edge);
            right_index = right_edge->get_start_site_index();

            (*ygap)[i_max-1][right_index].bwd_score = model->gap_close() * right_edge_wght;
        }
    }

}

void Simple_alignment::compute_fwd_viterbi_path(int i,int j)
{
    if(i==0 && j==0)
        return;

    int left_index = i;
    int right_index = j;

    Site * left_site;
    Site * right_site;

    Matrix_pointer *max_x = &(*xgap)[i][j];
    Matrix_pointer *max_y = &(*ygap)[i][j];
    Matrix_pointer *max_m = &(*match)[i][j];


    if(left_index>0)
    {
        left_site = left->get_site_at(left_index);

        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][j] ];
        align_slice y_slice = (*ygap)[ indices[ range( 0,ygap->shape()[0] ) ][j] ];
        align_slice m_slice = (*match)[ indices[ range( 0,match->shape()[0] ) ][j] ];

        this->iterate_bwd_edges_for_gap(left_site,&x_slice,&y_slice,&m_slice,max_x,true);
        max_x->y_ind = j;

    }
    else
    {
        max_x->x_ind = -1;
        max_x->y_ind = -1;
        max_x->matrix = -1;

        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

    if(right_index>0)
    {
        right_site = right->get_site_at(right_index);

        align_slice x_slice = (*xgap)[ indices[i][ range( 0,xgap->shape()[1] ) ] ];
        align_slice y_slice = (*ygap)[ indices[i][ range( 0,ygap->shape()[1] ) ] ];
        align_slice m_slice = (*match)[ indices[i][ range( 0,match->shape()[1] ) ] ];

        this->iterate_bwd_edges_for_gap(right_site,&y_slice,&x_slice,&m_slice,max_y,false);
        max_y->x_ind = i;

    }
    else
    {
        max_y->x_ind = -1;
        max_y->y_ind = -1;
        max_y->matrix = -1;

        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

    if(left_index>0 && right_index>0)
    {
        left_site = left->get_site_at(left_index);
        right_site = right->get_site_at(right_index);

        this->iterate_bwd_edges_for_match(left_site,right_site,max_m);

    }
    else
    {
        max_m->x_ind = -1;
        max_m->y_ind = -1;
        max_m->matrix = -1;
    }

}

/********************************************/

void Simple_alignment::compute_bwd_full_score(int i,int j)
{
    int i_end = match->shape()[0]-1;
    int j_end = match->shape()[1]-1;

    int left_index = i;
    int right_index = j;


    Site * left_site;
    Site * right_site;

    Matrix_pointer *max_x = &(*xgap)[i][j];
    Matrix_pointer *max_y = &(*ygap)[i][j];
    Matrix_pointer *max_m = &(*match)[i][j];

    if(left_index==i_end && right_index==j_end)
    {
        return;
    }

    if(left_index<i_end)
    {
        left_site = left->get_site_at(left_index);

        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][j] ];

        this->iterate_fwd_edges_for_gap(left_site,&x_slice,max_x,max_y,max_m);

    }

    if(right_index<j_end)
    {
        right_site = right->get_site_at(right_index);

        align_slice y_slice = (*ygap)[ indices[i][ range( 0,ygap->shape()[1] ) ] ];

        this->iterate_fwd_edges_for_gap(right_site,&y_slice,max_y,max_x,max_m);

    }

    if(left_index<i_end && right_index<j_end)
    {
        left_site = left->get_site_at(left_index);
        right_site = right->get_site_at(right_index);

        this->iterate_fwd_edges_for_match(left_site,right_site,max_x,max_y,max_m);

    }
}


void Simple_alignment::compute_posterior_score(int i,int j,double full_score)
{
    (*match)[i][j].full_score = (*match)[i][j].full_score * (*match)[i][j].bwd_score / full_score;
    (*xgap)[i][j].full_score = (*xgap)[i][j].full_score * (*xgap)[i][j].bwd_score / full_score;
    (*ygap)[i][j].full_score = (*ygap)[i][j].full_score * (*ygap)[i][j].bwd_score / full_score;
}

/********************************************/

void Simple_alignment::backtrack_new_path(vector<Path_pointer> *path,Path_pointer fp)
{
    vector<Edge> *left_edges = left->get_edges();
    vector<Edge> *right_edges = right->get_edges();

    int vit_mat = fp.mp.matrix;
    int x_ind = fp.mp.x_ind;
    int y_ind = fp.mp.y_ind;

    left_edges->at(fp.mp.x_edge_ind).is_used(true);
    right_edges->at(fp.mp.y_edge_ind).is_used(true);

    int j = match->shape()[1]-1;
    int i = match->shape()[0]-1;

    // Pre-existing gaps in the end skipped over
    //
    this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

    this->insert_new_path_pointer(path,&i,&j,fp);

    // Actual alignment path
    //
    while(j>=0)
    {
        while(i>=0)
        {
            if(vit_mat == Simple_alignment::m_mat)
            {

                vit_mat = (*match)[i][j].matrix;
                x_ind = (*match)[i][j].x_ind;
                y_ind = (*match)[i][j].y_ind;

                left_edges->at((*match)[i][j].x_edge_ind).is_used(true);
                right_edges->at((*match)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*match)[i][j], true );

                i--; j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::x_mat)
            {
                vit_mat = (*xgap)[i][j].matrix;
                x_ind = (*xgap)[i][j].x_ind;
                y_ind = (*xgap)[i][j].y_ind;

                left_edges->at((*xgap)[i][j].x_edge_ind).is_used(true);

                Path_pointer pp( (*xgap)[i][j], true );

                i--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::y_mat)
            {
                vit_mat = (*ygap)[i][j].matrix;
                x_ind = (*ygap)[i][j].x_ind;
                y_ind = (*ygap)[i][j].y_ind;

                right_edges->at((*ygap)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*ygap)[i][j], true );

                j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else
            {
                cout<<"incorrect backward pointer: "<<vit_mat<<endl;
                exit(-1);
            }

            if(i<1 && j<1)
                break;

        }

        if(i<1 && j<1)
            break;

    }

    /*DEBUG*/
    if(Settings::noise>4)
    {
        cout<<"\npath"<<endl;
        for(unsigned int i=0;i<path->size();i++)
            cout<<path->at(i).mp.matrix;
        cout<<endl;
    }
    /*DEBUG*/

}

/********************************************/

void Simple_alignment::sample_new_path(vector<Path_pointer> *path,Path_pointer fp)
{
    vector<Edge> *left_edges = left->get_edges();
    vector<Edge> *right_edges = right->get_edges();

    int vit_mat = fp.mp.matrix;
    int x_ind = fp.mp.x_ind;
    int y_ind = fp.mp.y_ind;

    left_edges->at(fp.mp.x_edge_ind).is_used(true);
    right_edges->at(fp.mp.y_edge_ind).is_used(true);

    int j = match->shape()[1]-1;
    int i = match->shape()[0]-1;

    // Pre-existing gaps in the end skipped over
    //
    this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

    this->insert_new_path_pointer(path,&i,&j,fp);

    Matrix_pointer bwd_p;

    // Actual alignment path
    //
    while(j>=0)
    {
        while(i>=0)
        {
//            cout<<i<<" "<<j<<": "<<x_ind<<" "<<y_ind<<" "<<vit_mat<<endl;
//
            if(vit_mat == Simple_alignment::m_mat)
            {
                this->iterate_bwd_edges_for_sampled_match(i,j,&bwd_p);

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

                left_edges->at(bwd_p.x_edge_ind).is_used(true);
                right_edges->at(bwd_p.y_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                i--; j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::x_mat)
            {
                this->iterate_bwd_edges_for_sampled_gap(i,j,&bwd_p,true);
                bwd_p.y_ind = j;

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

                left_edges->at(bwd_p.x_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                i--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::y_mat)
            {
                this->iterate_bwd_edges_for_sampled_gap(j,i,&bwd_p,false);
                bwd_p.x_ind = i;

                vit_mat = bwd_p.matrix;
                x_ind = bwd_p.x_ind;
                y_ind = bwd_p.y_ind;

                right_edges->at(bwd_p.y_edge_ind).is_used(true);

                Path_pointer pp( bwd_p, true );

                j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else
            {
                cout<<"incorrect backward pointer: "<<vit_mat<<endl;
                exit(-1);
            }

            if(i<1 && j<1)
                break;

        }

        if(i<1 && j<1)
            break;

    }

    /*DEBUG*/
    if(Settings::noise>4)
    {
        cout<<"\npath"<<endl;
        for(unsigned int i=0;i<path->size();i++)
            cout<<path->at(i).mp.matrix;
        cout<<endl;
    }
    /*DEBUG*/

}

/********************************************/

void Simple_alignment::build_ancestral_sequence(Sequence *sequence, vector<Path_pointer> *path)
{

    // The path is given as input.
    // This will create sites with correct child sites.
    this->create_ancestral_sequence(sequence,path);

    // This will add the edges connecting sites.
    this->create_ancestral_edges(sequence);

    // This will do a check-up and delete edges if needed:
    // " To mimic PRANK+F: one has to scan the sequence again to find boundaries 'Match/Skipped' and 'Skipped/Matched'
    //   and record them on the edges. If they are above a limit, the edges (and all between them) are deleted. "
    this->check_skipped_boundaries(sequence);


    if(Settings::noise>4)
    {
        cout<<"ANCESTRAL SEQUENCE:\n";
        sequence->print_sequence();
    }

    if(Settings::noise>6)
    {
        this->print_path(path);
    }

}

void Simple_alignment::create_ancestral_sequence(Sequence *sequence, vector<Path_pointer> *path)
{

    vector<Edge> *edges = sequence->get_edges();

    Site first_site( edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_children(0,0);

    sequence->push_back_site(first_site);

    int l_pos = 1;
    int r_pos = 1;

    for(unsigned int i=0;i<path->size();i++)
    {

        Site site( edges );
        site.set_empty_children();

        if(path->at(i).mp.matrix == Simple_alignment::x_mat)
        {
            int lc = left->get_site_at(l_pos)->get_state();
            site.set_state( lc );

            if(path->at(i).real_site)
                site.set_path_state( Site::xgapped );
            else
            {
                site.set_path_state( Site::xskipped );
                site.set_branch_count_since_last_used(
                        left->get_site_at(l_pos)->get_branch_count_since_last_used()+1 );
                site.set_branch_distance_since_last_used(
                        left->get_site_at(l_pos)->get_branch_distance_since_last_used()+left_branch_length );
            }
            site.get_children()->left_index = l_pos;

            l_pos++;
        }
        else if(path->at(i).mp.matrix == Simple_alignment::y_mat)
        {
            int rc = right->get_site_at(r_pos)->get_state();
            site.set_state( rc );

            if(path->at(i).real_site)
                site.set_path_state( Site::ygapped );
            else
            {
                site.set_path_state( Site::yskipped );
                site.set_branch_count_since_last_used(
                        right->get_site_at(r_pos)->get_branch_count_since_last_used()+1 );
                site.set_branch_distance_since_last_used(
                        right->get_site_at(r_pos)->get_branch_distance_since_last_used()+right_branch_length );
            }
            site.get_children()->right_index = r_pos;

            r_pos++;
        }
        else if(path->at(i).mp.matrix == Simple_alignment::m_mat)
        {
            int lc = left->get_site_at(l_pos)->get_state();
            int rc = right->get_site_at(r_pos)->get_state();
            site.set_state( model->parsimony_state(lc,rc) );

            site.set_path_state( Site::matched );

            site.set_children(l_pos,r_pos);

            l_pos++; r_pos++;
        }

        sequence->push_back_site(site);

        /*DEBUG*/
        if(Settings::noise>6)
            cout<<i<<": m "<<path->at(i).mp.matrix<<" si "<<sequence->get_current_site_index()<<": l "<<l_pos<<", r "<<r_pos<<" (st "<<site.get_state()<<")"<<endl;
        /*DEBUG*/

    }

    Site last_site( edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_children(left->sites_length()-1,right->sites_length()-1);
    sequence->push_back_site(last_site);

}

void Simple_alignment::create_ancestral_edges(Sequence *sequence)
{

    vector<Site> *sites = sequence->get_sites();

    vector<int> left_child_index;
    vector<int> right_child_index;

    // First create an index for the child's sites in the parent
    //
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        Site *lsite;
        Site *rsite;

        if(offspring->left_index>=0)
        {
            lsite = left->get_site_at(offspring->left_index);
            left_child_index.push_back(i);
        }
        if(offspring->right_index>=0)
        {
            rsite = right->get_site_at(offspring->right_index);
            right_child_index.push_back(i);
        }
    }

    // print the index
    if(Settings::noise>5)
    {
        cout<<"Child sequence site indeces:"<<endl;
        for(unsigned int i=0;i<left_child_index.size();i++)
            cout<<left_child_index.at(i)<<" ";

        cout<<endl;
        for(unsigned int i=0;i<right_child_index.size();i++)
            cout<<right_child_index.at(i)<<" ";
        cout<<endl;
    }

    // Then copy the edges of child sequences in their parent.
    // Additionally, create edges for cases where skipped gap is flanked by a new gap.
    //
    Edge_history prev(-1,-1);

    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *psite = &sites->at(i);
        int pstate = psite->get_path_state();

        Site_children *offspring = psite->get_children();

        // left sequence is matched
        if(offspring->left_index>=0)
        {
            Site *tsite = left->get_site_at(offspring->left_index);

            if( tsite->has_bwd_edge() )
            {
                Edge *child = tsite->get_first_bwd_edge();
                this->transfer_child_edge(sequence,child, &left_child_index, left_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(sequence,child, &left_child_index, left_branch_length );
                }
            }

            // these create edges to/from skipped sites flanked by gaps.
            if( (pstate == Site::matched || pstate == Site::ends_site ) && prev.left_skip_site_index >= 0)
            {
                // edge from the skipped site to the *next* site
                // as no better info is available, *this* edge is "copied" to one coming to the current site
                Edge query(prev.left_skip_site_index,prev.left_skip_site_index+1);
                int ind = left->get_fwd_edge_index_at_site(prev.left_skip_site_index,&query);

                if(ind>=0)
                {
                    Edge *child = &left->get_edges()->at(ind);
                    Edge edge( left_child_index.at(prev.left_skip_site_index), i );
                    this->transfer_child_edge(sequence,edge, child, left_branch_length );
                }

                prev.left_skip_site_index = -1;
            }
            else if(pstate == Site::xskipped && ( prev.path_state == Site::xgapped || prev.path_state == Site::ygapped ) )
            {
                // the same here: use this as a template for an extra edge
                Edge query(offspring->left_index-1, offspring->left_index);
                int ind = left->get_bwd_edge_index_at_site(offspring->left_index,&query);

                if(ind>=0)
                {
                    Edge *child = &left->get_edges()->at(ind);
                    Edge edge( prev.match_site_index, i );
                    this->transfer_child_edge(sequence,edge, child, left_branch_length );
                }

            }

            if(pstate == Site::xskipped)
                prev.left_skip_site_index = offspring->left_index;
            else
                prev.left_real_site_index = offspring->left_index;

            if(pstate == Site::matched)
                prev.match_site_index = i;
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right->get_site_at(offspring->right_index);

            if( tsite->has_bwd_edge() )
            {
                Edge *child = tsite->get_first_bwd_edge();
                this->transfer_child_edge(sequence,child, &right_child_index, right_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(sequence,child, &right_child_index, right_branch_length );
                }
            }

            if( (pstate == Site::matched || pstate == Site::ends_site ) && prev.right_skip_site_index >= 0)
            {
                Edge query(prev.right_skip_site_index,prev.right_skip_site_index+1);
                int ind = right->get_fwd_edge_index_at_site(prev.right_skip_site_index,&query);

                if(ind>=0)
                {
                    Edge *child = &right->get_edges()->at(ind);
                    Edge edge( right_child_index.at(prev.right_skip_site_index), i );
                    this->transfer_child_edge(sequence,edge, child, right_branch_length );
                }

                prev.right_skip_site_index = -1;
            }
            else if(pstate == Site::yskipped && ( prev.path_state == Site::xgapped || prev.path_state == Site::ygapped ) )
            {
                Edge query(offspring->right_index-1, offspring->right_index);
                int ind = right->get_bwd_edge_index_at_site(offspring->right_index,&query);

                if(ind>=0)
                {
                    Edge *child = &right->get_edges()->at(ind);
                    Edge edge( prev.match_site_index, i );
                    this->transfer_child_edge(sequence,edge, child, right_branch_length );
                }

            }

            if(pstate == Site::yskipped)
                prev.right_skip_site_index = offspring->right_index;
            else
                prev.right_real_site_index = offspring->right_index;

        }
        prev.path_state = pstate;
    }
}

void Simple_alignment::check_skipped_boundaries(Sequence *sequence)
{

    vector<Site> *sites = sequence->get_sites();

    // First, find 'Match/Skipped' and 'Skipped/Matched' boundaries and update the counts
    //
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if( tsite->has_bwd_edge() )
        {
            Edge *edge = tsite->get_first_bwd_edge();

            while( tsite->has_next_bwd_edge() )
            {
                Edge *another = tsite->get_next_bwd_edge();
                if( another->get_start_site_index() > edge->get_start_site_index() )
                    edge = another;
            }

            Site *psite = &sites->at( edge->get_start_site_index() );

            if( ( psite->get_path_state()==Site::matched || psite->get_path_state()==Site::start_site )
                && ( tsite->get_path_state()==Site::xskipped || tsite->get_path_state()==Site::yskipped ) )
            {
                edge->increase_branch_count_as_skipped_edge();
            }
        }

        if( tsite->has_fwd_edge() )
        {
            Edge *edge = tsite->get_first_fwd_edge();

            while( tsite->has_next_fwd_edge() )
            {
                Edge *another = tsite->get_next_fwd_edge();
                if( another->get_start_site_index() < edge->get_start_site_index() )
                    edge = another;
            }

            Site *nsite = &sites->at( edge->get_end_site_index() );

            if( ( tsite->get_path_state()==Site::xskipped || tsite->get_path_state()==Site::yskipped )
                && ( nsite->get_path_state()==Site::matched || nsite->get_path_state()==Site::ends_site ) )
            {
                edge->increase_branch_count_as_skipped_edge();
            }
        }
    }

    // Then, see if any pair of boundaries (covering a skipped gap) is above the limit. Delete the range.
    //
    bool non_skipped = true;
    int skip_start = -1;
    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);
        int tstate = tsite->get_path_state();

        if( non_skipped && ( tstate == Site::xskipped || tstate == Site::yskipped ) )
        {

            if( tsite->has_bwd_edge() )
            {
                Edge *edge = tsite->get_first_bwd_edge();
                while( tsite->has_next_bwd_edge() )
                {
                    Edge *another = tsite->get_next_bwd_edge();
                    if( another->get_start_site_index() > edge->get_start_site_index() )
                        edge = another;
                }
                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                {
                    skip_start = i;
                }
            }

            non_skipped = false;
        }

        if(!non_skipped && skip_start>=0 && tstate == Site::matched)
        {

            int edge_ind = -1;
            if( tsite->has_bwd_edge() )
            {
                Edge *edge = tsite->get_first_bwd_edge();

                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                    edge_ind = edge->get_index();

                while( tsite->has_next_bwd_edge() )
                {
                    edge = tsite->get_next_bwd_edge();

                    if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                        edge_ind = edge->get_index();

                }
            }

            if(edge_ind>=0)
            {
                if(Settings::noise>2) cout<<"delete: "<<edge_ind<<" "<<skip_start<<" "<<i<<endl;
                this->delete_edge_range(sequence,edge_ind,skip_start);
            }

            non_skipped = true;
            skip_start = -1;
        }

        if(tstate == Site::xgapped || tstate == Site::ygapped || tstate == Site::matched)
        {
            non_skipped = true;
            skip_start = -1;
        }
    }
}

void Simple_alignment::delete_edge_range(Sequence *sequence,int edge_ind,int skip_start_site)
{
    vector<Edge> *edges = sequence->get_edges();

    Edge *edge = &edges->at(edge_ind);

    int this_site_index = edge->get_start_site_index();

    // if not the start of the range, delete further
    while(this_site_index >= skip_start_site)
    {
//        cout<<"delete this: "<<this_site_index<<" ("<<skip_start_site<<")\n";
        sequence->delete_all_bwd_edges_at_site(this_site_index);
        sequence->delete_all_fwd_edges_at_site(this_site_index);
        --this_site_index;
    }

}

void Simple_alignment::transfer_child_edge(Sequence *sequence,Edge *child, vector<int> *child_index, float branch_length, bool adjust_posterior_weight, float branch_weight)
{
    Edge edge( child_index->at( child->get_start_site_index() ), child_index->at( child->get_end_site_index() ) );

    this->transfer_child_edge(sequence, edge, child, branch_length, adjust_posterior_weight, branch_weight);
}

void Simple_alignment::transfer_child_edge(Sequence *sequence, Edge edge, Edge *child, float branch_length, bool adjust_posterior_weight, float branch_weight)
{

    // No identical copies
    if(sequence->get_site_at( edge.get_end_site_index() )->contains_bwd_edge( &edge ) )
        return;

    // Limits for copying old edges:
    //  first, number of nodes since last used
    if(child->get_branch_count_since_last_used()+1 > max_allowed_skip_branches)
        return;

    //  then, total branch distance since last used
    if(child->get_branch_distance_since_last_used()+branch_length > max_allowed_skip_distance)
        return;


    // Comparison of distance and node count since last used to find boundaries of path branches.
    // Only start and end of an alternative path should be penalised; continuation on a path not.
    //
    float dist_start = sequence->get_site_at(edge.get_start_site_index())->get_branch_distance_since_last_used();
    float dist_end   = sequence->get_site_at(edge.get_end_site_index()  )->get_branch_distance_since_last_used();

    int count_start = sequence->get_site_at(edge.get_start_site_index())->get_branch_count_since_last_used();
    int count_end   = sequence->get_site_at(edge.get_end_site_index()  )->get_branch_count_since_last_used();

    // Sites on the two ends of an edge have different history: branch point that should be penalised
    if( dist_start != dist_end || count_start != count_end )
    {
        edge.set_branch_distance_since_last_used( max(dist_start,dist_end) );
        edge.set_branch_count_since_last_used( max(count_start,count_end) );

        if(adjust_posterior_weight)
            if(weighted_branch_skip_penalty)
                edge.set_weight( branch_weight * child->get_posterior_weight() * ( exp( -1.0 * this->branch_skip_weight * branch_length ) ) );
            else
                edge.set_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_probability );
        else
            edge.set_weight(child->get_posterior_weight());
    }
    // Edge is not used: just update the history
    else if(!child->is_used())
    {
        edge.set_branch_distance_since_last_used( child->get_branch_distance_since_last_used()+branch_length );
        edge.set_branch_count_since_last_used( child->get_branch_count_since_last_used()+1 );

        if(adjust_posterior_weight)
            if(weighted_branch_skip_penalty)
                edge.set_weight( branch_weight * child->get_posterior_weight() * ( exp( -1.0 * this->branch_skip_weight * branch_length ) ) );
            else
                edge.set_weight( branch_weight * child->get_posterior_weight() * this->branch_skip_probability );
        else
            edge.set_weight(child->get_posterior_weight());

    }

    if(!sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
    {

        edge.set_branch_count_as_skipped_edge( child->get_branch_count_as_skipped_edge() );
        sequence->push_back_edge(edge);

        sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( sequence->get_current_edge_index() );
        sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( sequence->get_current_edge_index() );
    }
}

/********************************************/

void Simple_alignment::iterate_bwd_edges_for_gap(Site * site,align_slice *z_slice,align_slice *w_slice,
                                                 align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix)
{
    if(site->has_bwd_edge()) {

        Edge * edge = site->get_first_bwd_edge();

        this->score_gap_ext(edge,z_slice,max,is_x_matrix);
        this->score_gap_double(edge,w_slice,max,is_x_matrix);
        this->score_gap_open(edge,m_slice,max,is_x_matrix);

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();

            this->score_gap_ext(edge,z_slice,max,is_x_matrix);
            this->score_gap_double(edge,w_slice,max,is_x_matrix);
            this->score_gap_open(edge,m_slice,max,is_x_matrix);
        }
    }
}

/********************************************/

void Simple_alignment::iterate_bwd_edges_for_match(Site * left_site,Site * right_site,Matrix_pointer *max)
{
    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        // match score & gap close scores for this match
        //
        double log_match_score = model->log_score(left_site->character_state,right_site->character_state);
        double m_log_match = 2*model->log_non_gap() + log_match_score;
        double x_log_match = model->log_gap_close() + model->log_non_gap() + log_match_score;
        double y_log_match = model->log_gap_close() + model->log_non_gap() + log_match_score;

        // start corner
//        if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//        {
////            m_log_match -= model->log_non_gap();
//        }
        // error in correction for bwd computation so let's forget it now.

        double m_match = 0;
        double x_match = 0;
        double y_match = 0;

        if(compute_full_score)
        {
            double match_score = model->score(left_site->character_state,right_site->character_state);
            m_match = model->non_gap() * model->non_gap() * match_score;

            // start corner
//            if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//            {
////                m_match /= model->non_gap();
//            }

            x_match = model->gap_close() * model->non_gap() * match_score;
            y_match = model->gap_close() * model->non_gap() * match_score;

        }

        this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
        this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
        this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
            this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
            this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

        }

        // left site extra edges first
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
            this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
            this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
                this->score_x_match(left_edge,right_edge,x_log_match,max,x_match);
                this->score_y_match(left_edge,right_edge,y_log_match,max,y_match);
            }
        }

    }
}

/********************************************/

void Simple_alignment::iterate_bwd_edges_for_end_corner(Site * left_site,Site * right_site,Matrix_pointer *max)
{

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        // match score & gap close scores for this match
        //
        double m_log_match = model->log_non_gap();
        double m_match = model->non_gap();

        this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
        double best_score = max->score;

        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
        this->score_gap_close(left_edge,&x_slice,max,true);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->y_ind = xgap->shape()[1]-1;
        }

        align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
        this->score_gap_close(right_edge,&y_slice,max,false);

        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->x_ind = ygap->shape()[0]-1;
        }

        // first right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
            this->score_gap_close(right_edge,&y_slice,max,false);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->x_ind = ygap->shape()[0]-1;
            }
        }

        // left site extra edges then
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
            }

            align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
            this->score_gap_close(left_edge,&x_slice,max,true);

            if(this->first_is_bigger(max->score,best_score) )
            {
                best_score = max->score;
                max->y_ind = xgap->shape()[1]-1;
            }


            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                }

                align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
                this->score_gap_close(right_edge,&y_slice,max,false);

                if(this->first_is_bigger(max->score,best_score) )
                {
                    best_score = max->score;
                    max->x_ind = ygap->shape()[0]-1;
                }

            }
        }

    }

    if(!compute_full_score && Settings::noise>5)
        this->print_matrices();


    if(Settings::noise>1)
    {
        cout<<"viterbi score: "<<setprecision(8)<<max->score<<" ["<<exp(max->score)<<"] ("<<max->matrix<<")";
        if(compute_full_score)
            cout<<", full probability: "<<log(max->full_score)<<" ["<<max->full_score<<"]";
        cout<<setprecision(4)<<endl; /*DEBUG*/
    }
}



/********************************************/

void Simple_alignment::iterate_fwd_edges_for_gap(Site * site,align_slice *g_slice,
                                   Matrix_pointer *max_s,Matrix_pointer *max_d,Matrix_pointer *max_m)
{
    if(site->has_fwd_edge()) {

        Edge * edge = site->get_first_fwd_edge();
        int slice_end = (int)g_slice->shape()[0];

        if(edge->get_end_site_index() < slice_end )
        {
            this->score_gap_ext_bwd(edge,g_slice,max_s);
            this->score_gap_double_bwd(edge,g_slice,max_d);
            this->score_gap_open_bwd(edge,g_slice,max_m);
        }
        while(site->has_next_fwd_edge())
        {
            edge = site->get_next_fwd_edge();

            if(edge->get_end_site_index() < slice_end )
            {
                this->score_gap_ext_bwd(edge,g_slice,max_s);
                this->score_gap_double_bwd(edge,g_slice,max_d);
                this->score_gap_open_bwd(edge,g_slice,max_m);
            }
        }
    }
}

/********************************************/

void Simple_alignment::iterate_fwd_edges_for_match(Site * left_site,Site * right_site,
                                     Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m)
{
    if(left_site->has_fwd_edge() && right_site->has_fwd_edge())
    {

        Edge * left_edge = left_site->get_first_fwd_edge();
        Edge * right_edge = right_site->get_first_fwd_edge();

        // PROBLEM: match score has to be computed for the cell where coming from, not where now!

        int left_end = (int)match->shape()[0];
        int right_end = (int)match->shape()[1];

        if(left_edge->get_end_site_index() < left_end
           && right_edge->get_end_site_index() < right_end )
        {
            this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
        }

        // then right site extra edges
        //
        while(right_site->has_next_fwd_edge())
        {
            right_edge = right_site->get_next_fwd_edge();
            left_edge = left_site->get_first_fwd_edge();

            if(left_edge->get_end_site_index() < left_end
               && right_edge->get_end_site_index() < right_end )
            {
                this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
            }
        }

        // left site extra edges first
        //
        while(left_site->has_next_fwd_edge())
        {
            left_edge = left_site->get_next_fwd_edge();
            right_edge = right_site->get_first_fwd_edge();

            if(left_edge->get_end_site_index() < left_end
               && right_edge->get_end_site_index() < right_end )
            {
                this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
            }

            while(right_site->has_next_fwd_edge())
            {

                right_edge = right_site->get_next_fwd_edge();

                if(left_edge->get_end_site_index() < left_end
                   && right_edge->get_end_site_index() < right_end )
                {
                    this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
                }
            }
        }
    }
}

/********************************************/

/********************************************/

void Simple_alignment::iterate_bwd_edges_for_sampled_gap(int site_index1,int site_index2,Matrix_pointer *sample_p,bool is_x_matrix)
{

    align_slice *z_slice;
    align_slice *w_slice;
    align_slice *m_slice;

    Site * site;

    if(is_x_matrix)
    {
        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][site_index2] ];
        align_slice y_slice = (*ygap)[ indices[ range( 0,ygap->shape()[0] ) ][site_index2] ];
        align_slice M_slice = (*match)[ indices[ range( 0,match->shape()[0] ) ][site_index2] ];

        z_slice = &x_slice;
        w_slice = &y_slice;
        m_slice = &M_slice;

        site = left->get_site_at(site_index1);
    }
    else
    {
        align_slice x_slice = (*xgap)[ indices[site_index2][ range( 0,xgap->shape()[1] ) ] ];
        align_slice y_slice = (*ygap)[ indices[site_index2][ range( 0,ygap->shape()[1] ) ] ];
        align_slice M_slice = (*match)[ indices[site_index2][ range( 0,match->shape()[1] ) ] ];

        z_slice = &y_slice;
        w_slice = &x_slice;
        m_slice = &M_slice;

        site = right->get_site_at(site_index1);
    }

    if(site->has_bwd_edge()) {

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;

        Edge * edge = site->get_first_bwd_edge();        

        Matrix_pointer bwd_p;
        this->add_sample_gap_ext(edge,z_slice,&bwd_p,is_x_matrix);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_gap_double(edge,w_slice,&bwd_p,is_x_matrix);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_gap_open(edge,m_slice,&bwd_p,is_x_matrix);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();

            this->add_sample_gap_ext(edge,z_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_gap_double(edge,w_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_gap_open(edge,m_slice,&bwd_p,is_x_matrix);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*sample_p) = bwd_pointers.at(i);
        if(is_x_matrix)
        {
            sample_p->y_ind = xgap->shape()[1]-1;
        }
        else
        {
            sample_p->x_ind = ygap->shape()[0]-1;
        }

    }
}

/********************************************/

void Simple_alignment::iterate_bwd_edges_for_sampled_match(int left_index,int right_index,Matrix_pointer *sample_p)
{
    Site * left_site = left->get_site_at(left_index);
    Site * right_site = right->get_site_at(right_index);

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;


        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        double match_score = model->score(left_site->character_state,right_site->character_state);
        double m_match = model->non_gap() * model->non_gap() * match_score;

            // start corner
//            if(left_edge->get_start_site_index() == 1 || right_edge->get_start_site_index() == 1)
//            {
////                m_match /= model->non_gap();
//            }

        double x_match = model->gap_close() * model->non_gap() * match_score;
        double y_match = model->gap_close() * model->non_gap() * match_score;


        Matrix_pointer bwd_p;

        this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );

        // then right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

        }

        // left site extra edges first
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

                this->add_sample_x_match(left_edge,right_edge,&bwd_p,x_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

                this->add_sample_y_match(left_edge,right_edge,&bwd_p,y_match);
                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );
            }
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*sample_p) = bwd_pointers.at(i);

    }
}

/********************************************/


void Simple_alignment::iterate_bwd_edges_for_sampled_end_corner(Site * left_site,Site * right_site,Matrix_pointer *end_p)
{

    if(left_site->has_bwd_edge() && right_site->has_bwd_edge())
    {

        Edge * left_edge = left_site->get_first_bwd_edge();
        Edge * right_edge = right_site->get_first_bwd_edge();

        vector<Matrix_pointer> bwd_pointers;
        double sum_score = 0;

        // match score & gap close scores for this match
        //
        double m_match = model->non_gap();

        Matrix_pointer bwd_p;

        this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
        this->add_sample_gap_close(left_edge,&x_slice,&bwd_p,true);
        bwd_p.y_ind = xgap->shape()[1]-1;

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
        this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
        bwd_p.x_ind = ygap->shape()[0]-1;

        sum_score  += bwd_p.score;
        bwd_pointers.push_back( bwd_p );


        // first right site extra edges
        //
        while(right_site->has_next_bwd_edge())
        {
            right_edge = right_site->get_next_bwd_edge();
            left_edge = left_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
            this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
            bwd_p.x_ind = ygap->shape()[0]-1;

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );

        }

        // left site extra edges then
        //
        while(left_site->has_next_bwd_edge())
        {
            left_edge = left_site->get_next_bwd_edge();
            right_edge = right_site->get_first_bwd_edge();

            this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
            this->add_sample_gap_close(left_edge,&x_slice,&bwd_p,true);
            bwd_p.y_ind = xgap->shape()[1]-1;

            sum_score  += bwd_p.score;
            bwd_pointers.push_back( bwd_p );


            while(right_site->has_next_bwd_edge())
            {

                right_edge = right_site->get_next_bwd_edge();

                this->add_sample_m_match(left_edge,right_edge,&bwd_p,m_match);

                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );


                align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
                this->add_sample_gap_close(right_edge,&y_slice,&bwd_p,false);
                bwd_p.x_ind = ygap->shape()[0]-1;

                sum_score  += bwd_p.score;
                bwd_pointers.push_back( bwd_p );

            }
        }

        double random_v = (  sum_score * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        int i = 0;
        double sum = bwd_pointers.at(i).score;
        while(sum < random_v)
        {
            ++i;
            sum += bwd_pointers.at(i).score;
        }

        (*end_p) = bwd_pointers.at(i);
    }

}


/********************************************/


/********************************************/

void Simple_alignment::score_m_match(Edge * left_edge,Edge * right_edge,double m_log_match,Matrix_pointer *max,double m_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*match)[left_prev_index][right_prev_index].score + m_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Simple_alignment::m_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*match)[left_prev_index][right_prev_index].full_score * m_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->full_score += this_full_score;
    }

}

void Simple_alignment::score_x_match(Edge * left_edge,Edge * right_edge,double x_log_match,Matrix_pointer *max, double x_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*xgap)[left_prev_index][right_prev_index].score + x_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Simple_alignment::x_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*xgap)[left_prev_index][right_prev_index].full_score * x_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->full_score += this_full_score;
    }
}

void Simple_alignment::score_y_match(Edge * left_edge,Edge * right_edge,double y_log_match,Matrix_pointer *max, double y_match)
{
    double left_edge_wght = this->get_log_edge_weight(left_edge);
    int left_prev_index = left_edge->get_start_site_index();

    double right_edge_wght = this->get_log_edge_weight(right_edge);
    int right_prev_index = right_edge->get_start_site_index();

    double this_score = (*ygap)[left_prev_index][right_prev_index].score + y_log_match + left_edge_wght + right_edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->x_ind = left_prev_index;
        max->y_ind = right_prev_index;
        max->x_edge_ind = left_edge->get_index();
        max->y_edge_ind = right_edge->get_index();
        max->matrix = Simple_alignment::y_mat;
    }

    if(compute_full_score)
    {
        double this_full_score =   (*ygap)[left_prev_index][right_prev_index].full_score * y_match
                                   * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);
        max->full_score += this_full_score;
    }
}

/************************************************/

void Simple_alignment::score_gap_ext(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix)
{
    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].score + model->log_gap_ext() + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Simple_alignment::x_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Simple_alignment::y_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*z_slice)[prev_index].full_score * model->gap_ext() * this->get_edge_weight(edge);
        max->full_score += this_full_score;
    }
}

void Simple_alignment::score_gap_double(Edge *edge,align_slice *w_slice,Matrix_pointer *max,bool is_x_matrix)
{
    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*w_slice)[prev_index].score + model->log_gap_close() + model->log_gap_open() + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Simple_alignment::y_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Simple_alignment::x_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*w_slice)[prev_index].full_score * model->gap_close() * model->gap_open() * this->get_edge_weight(edge);
        max->full_score += this_full_score;
    }

}

void Simple_alignment::score_gap_open(Edge *edge,align_slice *m_slice,Matrix_pointer *max,bool is_x_matrix)
{
    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*m_slice)[prev_index].score + model->log_non_gap() + model->log_gap_open() + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        max->matrix = Simple_alignment::m_mat;
        if(is_x_matrix)
        {
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*m_slice)[prev_index].full_score * model->non_gap() * model->gap_open() * this->get_edge_weight(edge);
        max->full_score += this_full_score;
    }

}

void Simple_alignment::score_gap_close(Edge *edge,align_slice *z_slice,Matrix_pointer *max,bool is_x_matrix)
{
    double edge_wght = this->get_log_edge_weight(edge);
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].score + model->log_gap_close() + edge_wght;

    if(this->first_is_bigger(this_score,max->score) )
    {
        max->score = this_score;
        if(is_x_matrix)
        {
            max->matrix = Simple_alignment::x_mat;
            max->x_ind = prev_index;
            max->x_edge_ind = edge->get_index();
        }
        else
        {
            max->matrix = Simple_alignment::y_mat;
            max->y_ind = prev_index;
            max->y_edge_ind = edge->get_index();
        }
    }

    if(compute_full_score)
    {
        double this_full_score =  (*z_slice)[prev_index].full_score * model->gap_close() * this->get_edge_weight(edge);
        max->full_score += this_full_score;
    }

}


/********************************************/
void Simple_alignment::score_match_bwd(Edge * left_edge,Edge * right_edge,
                                       Matrix_pointer *max_x,Matrix_pointer *max_y,Matrix_pointer *max_m)
{
    int left_prev_index = left_edge->get_end_site_index();
    int right_prev_index = right_edge->get_end_site_index();

    double m_match = model->non_gap() * model->non_gap();// * match_score;
    double x_match = model->gap_close() * model->non_gap();// * match_score;
    double y_match = model->gap_close() * model->non_gap();// * match_score;

    double match_score =
        model->score(left->get_site_at(left_prev_index)->get_state(),right->get_site_at(right_prev_index)->get_state());

    double this_full_score =   (*match)[left_prev_index][right_prev_index].bwd_score * match_score
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    max_x->bwd_score += this_full_score  * x_match;
    max_y->bwd_score += this_full_score  * y_match;
    max_m->bwd_score += this_full_score  * m_match;
}


/********************************************/

void Simple_alignment::score_gap_ext_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*z_slice)[prev_index].bwd_score * model->gap_ext() * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

void Simple_alignment::score_gap_double_bwd(Edge *edge,align_slice *w_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*w_slice)[prev_index].bwd_score * model->gap_close() * model->gap_open() * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

void Simple_alignment::score_gap_open_bwd(Edge *edge,align_slice *m_slice,Matrix_pointer *max)
{
    int prev_index = edge->get_end_site_index();

    double this_full_score =  (*m_slice)[prev_index].bwd_score * model->non_gap() * model->gap_open() * this->get_edge_weight(edge);
    max->bwd_score += this_full_score;
}

/********************************************/



/********************************************/

void Simple_alignment::add_sample_m_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double m_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*match)[left_prev_index][right_prev_index].full_score * m_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Simple_alignment::m_mat;
}

void Simple_alignment::add_sample_x_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double x_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*xgap)[left_prev_index][right_prev_index].full_score * x_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Simple_alignment::x_mat;
}

void Simple_alignment::add_sample_y_match(Edge * left_edge,Edge * right_edge,Matrix_pointer *bwd_p,double y_match)
{
    int left_prev_index = left_edge->get_start_site_index();
    int right_prev_index = right_edge->get_start_site_index();

    double this_score =   (*ygap)[left_prev_index][right_prev_index].full_score * y_match
                               * this->get_edge_weight(left_edge) * this->get_edge_weight(right_edge);

    bwd_p->score = this_score;
    bwd_p->x_ind = left_prev_index;
    bwd_p->y_ind = right_prev_index;
    bwd_p->x_edge_ind = left_edge->get_index();
    bwd_p->y_edge_ind = right_edge->get_index();
    bwd_p->matrix = Simple_alignment::y_mat;
}

/************************************************/

void Simple_alignment::add_sample_gap_ext(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].full_score * model->gap_ext() * this->get_edge_weight(edge);
    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Simple_alignment::x_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Simple_alignment::y_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }
}

void Simple_alignment::add_sample_gap_double(Edge * edge,align_slice *w_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*w_slice)[prev_index].full_score * model->gap_close() * model->gap_open() * this->get_edge_weight(edge);

    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Simple_alignment::y_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Simple_alignment::x_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}

void Simple_alignment::add_sample_gap_open(Edge * edge,align_slice *m_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*m_slice)[prev_index].full_score * model->non_gap() * model->gap_open() * this->get_edge_weight(edge);

    bwd_p->score = this_score;
    bwd_p->matrix = Simple_alignment::m_mat;

    if(is_x_matrix)
    {
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}

void Simple_alignment::add_sample_gap_close(Edge * edge,align_slice *z_slice,Matrix_pointer *bwd_p,bool is_x_matrix)
{
    int prev_index = edge->get_start_site_index();

    double this_score =  (*z_slice)[prev_index].full_score * model->gap_close() * this->get_edge_weight(edge);

    bwd_p->score = this_score;

    if(is_x_matrix)
    {
        bwd_p->matrix = Simple_alignment::x_mat;
        bwd_p->x_ind = prev_index;
        bwd_p->x_edge_ind = edge->get_index();
    }
    else
    {
        bwd_p->matrix = Simple_alignment::y_mat;
        bwd_p->y_ind = prev_index;
        bwd_p->y_edge_ind = edge->get_index();
    }

}


/********************************************/




/********************************************/
int Simple_alignment::plot_number = 1;

void Simple_alignment::plot_posterior_probabilities_down()
{
    string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

    string path = file;
    path.append(".mp");

    string path2 = file;
    path2.append(".tex");

    ofstream output(path.c_str(), (ios::out|ios::app));
    if (! output) { throw IOException ("Simple_alignment::plot_posterior_probabilities. Failed to open file"); }

    ofstream output2(path2.c_str(), (ios::out|ios::app));
    if (! output2) { throw IOException ("Simple_alignment::plot_posterior_probabilities. Failed to open file"); }

    output<<"beginfig("<<plot_number<<");\n"
            "u := 1mm;\ndefaultscale := 1.5pt/fontsize(defaultfont);\n"
            "path l[]; path r[]; path sqr; sqr := unitsquare scaled u;\n"
            "pickup pencircle scaled 0.1; labeloffset := 1bp; \n";

    for(unsigned int j=1;j<match->shape()[1];j++)
    {
        for(unsigned int i=1;i<match->shape()[0];i++)
        {
            if((*match)[i][j].full_score>0){
                int score = int( abs( log( (*match)[i][j].full_score ) ) );
                float red = 1;

                float green = score*7;
                if(green>255) green = 255;
                green /= 255;

                float blue = (score-39)*7;
                if(blue<0) blue = 0;
                if(blue>255) blue = 255;
                blue /= 255;

                output<<"fill sqr shifted ("<<i<<"*u,-"<<j<<"*u)\n";
                output<<"withcolor ("<<red<<","<<green<<","<<blue<<");\n";
            }
            output<<"draw sqr shifted ("<<i<<"*u,-"<<j<<"*u);\n";

        }
    }

    output<<"label.lft(\"#"<<plot_number<<"#\" infont defaultfont scaled 0.2,(0mm,1.5mm));\n";

    vector<Site> *left_sites = left->get_sites();
    vector<Site> *right_sites = right->get_sites();
    string full_alphabet = left->get_full_alphabet();

    for(unsigned int i=0;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"l"<<i<<" = circle"<<"(("<<0.5+i<<"mm,1.5mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 0;
            string place = "edgetop";
            if(start+1==stop)
                place = "edgebot";
            else if(start+2==stop)
                angle = 45;
            else if(start+3==stop)
                angle = 35;
            else if(start+4<=stop)
                angle = 25;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 0;
                place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 45;
                else if(start+3==stop)
                    angle = 35;
                else if(start+4<=stop)
                    angle = 25;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    for(unsigned int i=0;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"r"<<i<<" = circle"<<"((-0.5mm,"<<0.5-i<<"mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 270;
            string place = "edgelft";
            if(start+1==stop)
                place = "edgergt";
            else if(start+2==stop)
                angle = 225;
            else if(start+3==stop)
                angle = 235;
            else if(start+4<=stop)
                angle = 245;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 270;
                place = "edgelft";
                if(start+1==stop)
                    place = "edgergt";
                else if(start+2==stop)
                    angle = 225;
                else if(start+3==stop)
                    angle = 235;
                else if(start+4<=stop)
                    angle = 245;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    vector<Site> *sites = ancestral_sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if(tsite->get_site_type()==Site::real_site)
        {

            Site_children *offspring = sites->at(i).get_children();
            int lc = offspring->left_index;
            int rc = offspring->right_index;

            if(lc>=0 && rc>=0)
            {
                output<<"draw sqr shifted ("<<lc<<"*u,-"<<rc<<"*u) withcolor (0,0,1) withpen pencircle scaled 0.5;\n";
            }
        }
    }

    output<<"endfig;\n";
    output.close();

    output2<<"\\includegraphics[width=0.9\\columnwidth]{"<<file<<"."<<plot_number<<"}\n\n\\bigskip\n";
    output2<<"~\n\n\\bigskip\n";

    output2.close();

    Simple_alignment::plot_number += 1;
}

void Simple_alignment::plot_posterior_probabilities_up()
{
    string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

    string path = file;
    path.append(".mp");

    string path2 = file;
    path2.append(".tex");

    ofstream output(path.c_str(), (ios::out|ios::app));
    if (! output) { throw IOException ("Simple_alignment::plot_posterior_probabilities. Failed to open file"); }

    ofstream output2(path2.c_str(), (ios::out|ios::app));
    if (! output2) { throw IOException ("Simple_alignment::plot_posterior_probabilities. Failed to open file"); }

    output<<"beginfig("<<plot_number<<");\n"
            "u := 1mm;\ndefaultscale := 1.5pt/fontsize(defaultfont);\n"
            "path l[]; path r[]; path sqr; sqr := unitsquare scaled u;\n"
            "pickup pencircle scaled 0.1; labeloffset := 1bp; \n";

    for(unsigned int j=1;j<match->shape()[1];j++)
    {
        for(unsigned int i=1;i<match->shape()[0];i++)
        {
            if((*match)[i][j].full_score>0){
                int score = int( abs( log( (*match)[i][j].full_score ) ) );
                float red = 1;

                float green = score*7;
                if(green>255) green = 255;
                green /= 255;

                float blue = (score-39)*7;
                if(blue<0) blue = 0;
                if(blue>255) blue = 255;
                blue /= 255;

                output<<"fill sqr shifted ("<<i<<"*u,"<<j<<"*u)\n";
                output<<"withcolor ("<<red<<","<<green<<","<<blue<<");\n";
            }
            output<<"draw sqr shifted ("<<i<<"*u,"<<j<<"*u);\n";

        }
    }

    output<<"label.lft(\"#"<<plot_number<<"#\" infont defaultfont scaled 0.2,(0mm,-0.5mm));\n";

    vector<Site> *left_sites = left->get_sites();
    vector<Site> *right_sites = right->get_sites();
    string full_alphabet = left->get_full_alphabet();

    for(unsigned int i=0;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"l"<<i<<" = circle"<<"(("<<0.5+i<<"mm,-0.5mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<left_sites->size();i++)
    {
        Site *tsite =  &left_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 0;
            string place = "edgebot";
            if(start+1==stop)
                place = "edgetop";
            else if(start+2==stop)
                angle = 45;
            else if(start+3==stop)
                angle = 35;
            else if(start+4<=stop)
                angle = 25;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 0;
                place = "edgebot";
                if(start+1==stop)
                    place = "edgetop";
                else if(start+2==stop)
                    angle = 45;
                else if(start+3==stop)
                    angle = 35;
                else if(start+4<=stop)
                    angle = 25;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    for(unsigned int i=0;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = Node::get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        output<<"r"<<i<<" = circle"<<"((-0.5mm,"<<0.5+i<<"mm),\""<<c<<"\","<<color<<");\n";
    }

    for(unsigned int i=1;i<right_sites->size();i++)
    {
        Site *tsite =  &right_sites->at(i);

        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            int start = tedge->get_start_site_index();
            int stop  = tedge->get_end_site_index();

            int angle = 90;
            string place = "edgelft";
            if(start+1==stop)
                place = "edgergt";
            else if(start+2==stop)
                angle = 135;
            else if(start+3==stop)
                angle = 125;
            else if(start+4<=stop)
                angle = 115;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                start = tedge->get_start_site_index();
                stop  = tedge->get_end_site_index();

                angle = 90;
                place = "edgelft";
                if(start+1==stop)
                    place = "edgergt";
                else if(start+2==stop)
                    angle = 135;
                else if(start+3==stop)
                    angle = 125;
                else if(start+4<=stop)
                    angle = 115;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    vector<Site> *sites = ancestral_sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);

        if(tsite->get_site_type()==Site::real_site)
        {

            Site_children *offspring = sites->at(i).get_children();
            int lc = offspring->left_index;
            int rc = offspring->right_index;

            if(lc>=0 && rc>=0)
            {
                output<<"draw sqr shifted ("<<lc<<"*u,"<<rc<<"*u) withcolor (0,0,1) withpen pencircle scaled 0.5;\n";
            }
        }
    }

    output<<"endfig;\n";
    output.close();

    output2<<"\\includegraphics[width=0.9\\columnwidth]{"<<file<<"."<<plot_number<<"}\n\n\\bigskip\n";
    output2<<"~\n\n\\bigskip\n";

    output2.close();

    Simple_alignment::plot_number += 1;
}

/********************************************/

void Simple_alignment::print_matrices()
{
    cout << fixed << noshowpos << setprecision (4);

    cout<<"m"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<(*match)[i][j].matrix<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"m"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<setw(8)<<fixed<<(*match)[i][j].score<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    if(compute_full_score)
    {
        cout<<"m"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*match)[i][j].full_score<<" ";
                cout<<setw(8)<<log((*match)[i][j].full_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;

        cout<<"m"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*match)[i][j].bwd_score<<" ";
                cout<<setw(8)<<log((*match)[i][j].bwd_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }

    //    cout<<"m"<<endl;
//    for(unsigned int j=0;j<match->shape()[1];j++)
//    {
//        for(unsigned int i=0;i<match->shape()[0];i++)
//        {
//            cout<<setw(6)<<"("<<(*match)[i][j].x_ind<<","<<(*match)[i][j].y_ind<<") ";
//        }
//        cout<<endl;
//    }
//    cout<<endl;

    cout<<"x"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<(*xgap)[i][j].matrix<<" ";
        }
        cout<<endl;
    }
    cout<<endl;


    cout<<"x"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<setw(8)<<fixed<<(*xgap)[i][j].score<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    if(compute_full_score)
    {
        cout<<"x"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*xgap)[i][j].full_score<<" ";
                cout<<setw(8)<<log((*xgap)[i][j].full_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;

        cout<<"x"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*xgap)[i][j].bwd_score<<" ";
                cout<<setw(8)<<log((*xgap)[i][j].bwd_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }

//    cout<<"x"<<endl;
//    for(unsigned int j=0;j<match->shape()[1];j++)
//    {
//        for(unsigned int i=0;i<match->shape()[0];i++)
//        {
//            cout<<setw(6)<<"("<<(*xgap)[i][j].x_ind<<","<<(*xgap)[i][j].y_ind<<") ";
//        }
//        cout<<endl;
//    }
//    cout<<endl;

        cout<<"y"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<(*ygap)[i][j].matrix<<" ";
        }
        cout<<endl;
    }
    cout<<endl;



    cout<<"y"<<endl;
    for(unsigned int j=0;j<match->shape()[1];j++)
    {
        for(unsigned int i=0;i<match->shape()[0];i++)
        {
            cout<<setw(8)<<fixed<<(*ygap)[i][j].score<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    if(compute_full_score)
    {
        cout<<"y"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*ygap)[i][j].full_score<<" ";
                cout<<setw(8)<<log((*ygap)[i][j].full_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;

        cout<<"y"<<endl;
        for(unsigned int j=0;j<match->shape()[1];j++)
        {
            for(unsigned int i=0;i<match->shape()[0];i++)
            {
//                cout<<setw(8)<<scientific<<(*ygap)[i][j].bwd_score<<" ";
                cout<<setw(8)<<log((*ygap)[i][j].bwd_score)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }

//    cout<<"y"<<endl;
//    for(unsigned int j=0;j<match->shape()[1];j++)
//    {
//        for(unsigned int i=0;i<match->shape()[0];i++)
//        {
//            cout<<setw(6)<<"("<<(*ygap)[i][j].x_ind<<","<<(*ygap)[i][j].y_ind<<") ";
//        }
//        cout<<endl;
//    }
//    cout<<endl;
}

void Simple_alignment::print_sequences(vector<Site> *sites)
{
    /*DEBUG -->*/
    cout<<endl;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite = &sites->at(i);
        cout<<i<<" "<<tsite->get_state()<<endl;

        if(tsite->get_site_type()==Site::real_site)
            cout<<full_dna_alphabet.at(tsite->get_state())<<": ";
        else
            cout<<"+: ";

        Site_children *offspring = sites->at(i).get_children();
        int lc = -1;
        int rc = -1;
        Site *lsite = 0;
        Site *rsite = 0;


        if(offspring->left_index>=0)
        {
            lsite = left->get_site_at(offspring->left_index);
            lc = lsite->get_state();
        }
        if(offspring->right_index>=0)
        {
            rsite = right->get_site_at(offspring->right_index);
            rc = rsite->get_state();
        }

        if(lc>=0 && rc>=0)
            cout<<full_dna_alphabet.at(lc)<<full_dna_alphabet.at(rc)<<" ";
        else if(lc>=0)
            cout<<full_dna_alphabet.at(lc)<<"- ";
        else if(rc>=0)
            cout<<"-"<<full_dna_alphabet.at(rc)<<" ";


        cout<<" (P) ";
        if(tsite->has_fwd_edge())
        {
            Edge *tedge = tsite->get_first_fwd_edge();
            cout<<" F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
            while(tsite->has_next_fwd_edge())
            {
                tedge = tsite->get_next_fwd_edge();
                cout<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
            }
        }
        if(tsite->has_bwd_edge())
        {
            Edge *tedge = tsite->get_first_bwd_edge();
            cout<<"; B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
            while(tsite->has_next_bwd_edge())
            {
                tedge = tsite->get_next_bwd_edge();
                cout<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
            }
        }

        if(lc>=0)
        {
            tsite = lsite;
            cout<<"; (L) ";
            if(tsite->has_fwd_edge())
            {
                Edge *tedge = tsite->get_first_fwd_edge();
                cout<<"; F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                while(tsite->has_next_fwd_edge())
                {
                    tedge = tsite->get_next_fwd_edge();
                    cout<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                }
            }
            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                cout<<"; B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    cout<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                }
            }
        }

        if(rc>=0)
        {
            tsite = rsite;
            cout<<"; (R) ";
            if(tsite->has_fwd_edge())
            {
                Edge *tedge = tsite->get_first_fwd_edge();
                cout<<"; F "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                while(tsite->has_next_fwd_edge())
                {
                    tedge = tsite->get_next_fwd_edge();
                    cout<<"; f "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                }
            }
            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                cout<<"; B "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    cout<<"; b "<<tedge->get_start_site_index()<<" "<<tedge->get_end_site_index();
                }
            }
        }
        cout<<endl;
    }
    /*<-- DEBUG*/

}
