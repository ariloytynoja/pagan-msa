#include "main/simple_alignment.h"
#include "main/sequence.h"
#include <iomanip>

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


    if(Settings::noise>2)  /*DEBUG*/
        this->print_input_sequences();


    // set the basic parameters (copy from Settings)
    //
    this->set_basic_settings();

    // Set the edge weighting scheme
    //
    log_edge_weight = &ppa::Simple_alignment::edge_log_posterior_weight;
    edge_weight = &ppa::Simple_alignment::edge_posterior_weight;


    int left_length = left_sequence->sites_length();
    int right_length = right_sequence->sites_length();

    if(Settings::noise>1)  /*DEBUG*/
        cout<<"Simple_alignment: lengths: "<<left_length<<" "<<right_length<<endl;

    align_array M(boost::extents[left_length-1][right_length-1]);
    align_array X(boost::extents[left_length-1][right_length-1]);
    align_array Y(boost::extents[left_length-1][right_length-1]);

    if(Settings::noise>1)
        cout<<"Simple_alignment: matrix created\n";

    match = &M;
    xgap = &X;
    ygap = &Y;

    this->initialise_array_corner();

    int j_max = match->shape()[1];
    int i_max = match->shape()[0];

    // Dynamic programming loop
    //
    for(int j=0;j<j_max;j++)
    {
        for(int i=0;i<i_max;i++)
        {
            this->find_best_transition(i,j);
        }
    }

    if(Settings::noise>1)
        cout<<"Simple_alignment: matrix filled\n";

    if(compute_full_score)
    {
        this->initialise_array_corner_bwd();

        // do the same loop backwards to compute the full probability
        for(int j=j_max-1;j>=0;j--)
        {
            for(int i=i_max-1;i>=0;i--)
            {
                this->compute_bwd_full_score(i,j);
            }
        }
    }

    // It is more convenient to build the sequence forward;
    // thus, backtrace the path into a vector and use that in the next step
    //
    Site * left_site = left->get_site_at(i_max);
    Site * right_site = right->get_site_at(j_max);

    Matrix_pointer max_end;
    this->iterate_bwd_edges_for_end_corner(left_site,right_site,&max_end);

    if(Settings::noise>1)
        cout<<"Simple_alignment: corner found\n";

    Matrix_pointer max_start;
    if(compute_full_score)
    {
        max_start = (*match)[0][0];

        if(Settings::noise>1)
        {
            cout<<" bwd full probability: "<<setprecision(8)<<log(max_start.bwd_score)<<" ["<<max_start.bwd_score<<"]"<<setprecision(4)<<endl; /*DEBUG*/
        }
    
        double ratio = max_end.full_score/max_start.bwd_score;

        if(Settings::noise>0 && (ratio<0.99 || ratio>1.01))
        {
            cout<<"fwd: "<<max_end.full_score<<", bwd: "<<max_start.bwd_score<<endl;
        }
    }

    Path_pointer pp(max_end,true);

    this->backtrack_new_path(&path,pp);

    if(Settings::noise>1)
        cout<<"Simple_alignment: path found\n";


    // Now build the sequence forward following the path saved in a vector;
    // This doesn't copy the edges in the child sequences; they may need
    // extra information (e.g., brach length) and need to be done later

    this->build_ancestral_sequence(&path);

    if(Settings::noise>1)
        cout<<"Simple_alignment: sequence built\n";

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

void Simple_alignment::find_best_transition(int i,int j)
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
//if(Settings::noise>6) cout<<i<<" "<<j<<"\n ";
//cout<<i<<" "<<j;
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

//if(Settings::noise>6) cout<<i<<" "<<j<<" "<<scientific<<max_m->bwd_score<<" before\n ";
if(Settings::noise>6) cout<<i<<" "<<j<<" "<<log(max_m->bwd_score)<<" before\n ";

    if(left_index<i_end)
    {
//    cout<<" x";
        left_site = left->get_site_at(left_index);

        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][j] ];

        this->iterate_fwd_edges_for_gap(left_site,&x_slice,max_x,max_y,max_m);

    }

    if(right_index<j_end)
    {
//    cout<<" y";
        right_site = right->get_site_at(right_index);

        align_slice y_slice = (*ygap)[ indices[i][ range( 0,ygap->shape()[1] ) ] ];

        this->iterate_fwd_edges_for_gap(right_site,&y_slice,max_y,max_x,max_m);

    }

    if(left_index<i_end && right_index<j_end)
    {
//    cout<<" m";
        left_site = left->get_site_at(left_index);
        right_site = right->get_site_at(right_index);

        this->iterate_fwd_edges_for_match(left_site,right_site,max_x,max_y,max_m);

    }
//if(Settings::noise>6) cout<<i<<" "<<j<<" "<<scientific<<max_m->bwd_score<<" after\n ";
if(Settings::noise>6) cout<<i<<" "<<j<<" "<<log(max_m->bwd_score)<<" after\n ";
//cout<<endl;
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
//            cout<<i<<" "<<j<<": ";
            if(vit_mat == Simple_alignment::m_mat)
            {

                vit_mat = (*match)[i][j].matrix;
                x_ind = (*match)[i][j].x_ind;
                y_ind = (*match)[i][j].y_ind;

                left_edges->at((*match)[i][j].x_edge_ind).is_used(true);
                right_edges->at((*match)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*match)[i][j], true );
//if(Settings::noise>2){ cout<<"m "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }

                i--; j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);
//if(Settings::noise>2){ cout<<"m "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }
                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::x_mat)
            {
                vit_mat = (*xgap)[i][j].matrix;
                x_ind = (*xgap)[i][j].x_ind;
                y_ind = (*xgap)[i][j].y_ind;

                left_edges->at((*xgap)[i][j].x_edge_ind).is_used(true);

                Path_pointer pp( (*xgap)[i][j], true );
//if(Settings::noise>2){ cout<<"x "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }

                i--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

//if(Settings::noise>2){ cout<<"x "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }
                this->insert_new_path_pointer(path,&i,&j,pp);
            }
            else if(vit_mat == Simple_alignment::y_mat)
            {
                vit_mat = (*ygap)[i][j].matrix;
                x_ind = (*ygap)[i][j].x_ind;
                y_ind = (*ygap)[i][j].y_ind;

                right_edges->at((*ygap)[i][j].y_edge_ind).is_used(true);

                Path_pointer pp( (*ygap)[i][j], true );
//if(Settings::noise>2){ cout<<"y "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }

                j--;

                // Pre-existing gaps in the middle skipped over
                //
                this->insert_preexisting_gap(path, &i, &j, x_ind, y_ind);

//if(Settings::noise>2){ cout<<"y "<<i<<" "<<j<<": "<<pp.mp.x_ind<<" "<<pp.mp.y_ind<<" "<<pp.mp.matrix<<": "<<pp.mp.score<<endl; }
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

void Simple_alignment::build_ancestral_sequence(vector<Path_pointer> *path)
{

    // The path is given as input.
    // This will create sites with correct child sites.
    this->create_ancestral_sequence(path);
//cout<<"seq created\n";

    // This will add the edges connecting sites.
    this->create_ancestral_edges();
//cout<<"edges created\n";

    // This will do a check-up and delete edges if needed:
    // " To mimic PRANK+F: one has to scan the sequence again to find boundaries 'Match/Skipped' and 'Skipped/Matched'
    //   and record them on the edges. If they are above a limit, the edges (and all between them) are deleted. "
    this->check_skipped_boundaries();
//cout<<"skips checked\n";


    if(Settings::noise>4)
    {
        cout<<"ANCESTRAL SEQUENCE:\n";
        ancestral_sequence->print_sequence();
    }

    if(Settings::noise>6)
    {
        this->print_path(path);
    }

}

void Simple_alignment::create_ancestral_sequence(vector<Path_pointer> *path)
{

    ancestral_sequence = new Sequence(path->size(),model->get_full_alphabet());

    vector<Edge> *edges = ancestral_sequence->get_edges();

    Site first_site( edges, Site::start_site, Site::ends_site );
    first_site.set_state( -1 );
    first_site.set_children(0,0);

    ancestral_sequence->push_back_site(first_site);

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

        ancestral_sequence->push_back_site(site);

        /*DEBUG*/
        if(Settings::noise>6)
            cout<<i<<": m "<<path->at(i).mp.matrix<<" si "<<ancestral_sequence->get_current_site_index()<<": l "<<l_pos<<", r "<<r_pos<<" (st "<<site.get_state()<<")"<<endl;
        /*DEBUG*/

    }

    Site last_site( edges, Site::stop_site, Site::ends_site );
    last_site.set_state( -1 );
    last_site.set_children(left->sites_length()-1,right->sites_length()-1);
    ancestral_sequence->push_back_site(last_site);

}

void Simple_alignment::create_ancestral_edges()
{

    vector<Site> *sites = ancestral_sequence->get_sites();

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
                this->transfer_child_edge(child, &left_child_index, left_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(child, &left_child_index, left_branch_length );
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
                    this->transfer_child_edge(edge, child, left_branch_length );
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
//                    Edge edge( left_child_index.at(prev.match_site_index), i );
                    this->transfer_child_edge(edge, child, left_branch_length );
                }

            }

            if(pstate == Site::xskipped)
                prev.left_skip_site_index = offspring->left_index;
            else
                prev.left_real_site_index = offspring->left_index;

            if(pstate == Site::matched)
                prev.match_site_index = i;
//                prev.match_site_index = offspring->left_index;
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right->get_site_at(offspring->right_index);

            if( tsite->has_bwd_edge() )
            {
                Edge *child = tsite->get_first_bwd_edge();
                this->transfer_child_edge(child, &right_child_index, right_branch_length );

                while( tsite->has_next_bwd_edge() )
                {
                    child = tsite->get_next_bwd_edge();
                    this->transfer_child_edge(child, &right_child_index, right_branch_length );
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
                    this->transfer_child_edge(edge, child, right_branch_length );
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
//                    Edge edge( right_child_index.at(prev.match_site_index), i );
                    Edge edge( prev.match_site_index, i );
                    this->transfer_child_edge(edge, child, right_branch_length );
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

void Simple_alignment::check_skipped_boundaries()
{

    vector<Site> *sites = ancestral_sequence->get_sites();

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

            if( psite->get_path_state()==Site::matched
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
                && nsite->get_path_state()==Site::matched )
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
//                Edge *edge = tsite->get_first_bwd_edge();
//                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
//                    skip_start = i;
//
//                while( tsite->has_next_bwd_edge() )
//                {
//                    edge = tsite->get_next_bwd_edge();
//                    if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
//                        skip_start = i;
//                }
                Edge *edge = tsite->get_first_bwd_edge();
                while( tsite->has_next_bwd_edge() )
                {
                    Edge *another = tsite->get_next_bwd_edge();
                    if( another->get_start_site_index() > edge->get_start_site_index() )
                        edge = another;
                }
                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                    skip_start = i;

            }

            non_skipped = false;
        }
        if(!non_skipped && skip_start>=0 && tstate == Site::matched)
        {

            int edge_ind = -1;
            if( tsite->has_bwd_edge() )
            {
//                Edge *edge = tsite->get_first_bwd_edge();
//                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
//                    edge_ind = edge->get_index();
//
//                while( tsite->has_next_bwd_edge() )
//                {
//                    edge = tsite->get_next_bwd_edge();
//                    if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
//                        edge_ind = edge->get_index();
//                }

                Edge *edge = tsite->get_first_bwd_edge();

                while( tsite->has_next_bwd_edge() )
                {
                    Edge *another = tsite->get_next_bwd_edge();
                    if( another->get_start_site_index() < edge->get_start_site_index() )
                        edge = another;
                }

                if(edge->get_branch_count_as_skipped_edge()>max_allowed_match_skip_branches)
                    edge_ind = edge->get_index();

            }

            if(edge_ind>=0)
            {
                cout<<"delete: "<<edge_ind<<" "<<skip_start<<" "<<i<<endl;
                this->delete_edge_range(edge_ind,skip_start);
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

void Simple_alignment::delete_edge_range(int edge_ind,int skip_start_site)
{
    vector<Edge> *edges = ancestral_sequence->get_edges();

    Edge *edge = &edges->at(edge_ind);

    int this_site_index = edge->get_start_site_index();

    // if not the start of the range, delete further
    if(this_site_index >= skip_start_site)
    {
        ancestral_sequence->delete_all_bwd_edges_at_site(this_site_index);
        ancestral_sequence->delete_all_fwd_edges_at_site(this_site_index);
        --this_site_index;
    }

}

void Simple_alignment::transfer_child_edge(Edge *child, vector<int> *child_index, float branch_length, bool adjust_posterior_weight, float branch_weight)
{
    Edge edge( child_index->at( child->get_start_site_index() ), child_index->at( child->get_end_site_index() ) );

    this->transfer_child_edge(edge, child, branch_length, adjust_posterior_weight, branch_weight);
}

void Simple_alignment::transfer_child_edge(Edge edge, Edge *child, float branch_length, bool adjust_posterior_weight, float branch_weight)
{

    // No identical copies
    if(ancestral_sequence->get_site_at( edge.get_end_site_index() )->contains_bwd_edge( &edge ) )
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
    float dist_start = ancestral_sequence->get_site_at(edge.get_start_site_index())->get_branch_distance_since_last_used();
    float dist_end   = ancestral_sequence->get_site_at(edge.get_end_site_index()  )->get_branch_distance_since_last_used();

    int count_start = ancestral_sequence->get_site_at(edge.get_start_site_index())->get_branch_count_since_last_used();
    int count_end   = ancestral_sequence->get_site_at(edge.get_end_site_index()  )->get_branch_count_since_last_used();

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

    if(!ancestral_sequence->contains_this_bwd_edge_at_site(edge.get_end_site_index(),&edge))
    {

        edge.set_branch_count_as_skipped_edge( child->get_branch_count_as_skipped_edge() );
        ancestral_sequence->push_back_edge(edge);

        ancestral_sequence->get_site_at( edge.get_start_site_index() )->add_new_fwd_edge_index( ancestral_sequence->get_current_edge_index() );
        ancestral_sequence->get_site_at( edge.get_end_site_index()   )->add_new_bwd_edge_index( ancestral_sequence->get_current_edge_index() );
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
//cout<< setprecision (8)<<"f0 "<<max->full_score<<endl;
        this->score_m_match(left_edge,right_edge,m_log_match,max,m_match);
        double best_score = max->score;
//cout<<"f1 "<<max->full_score<<endl;
        align_slice x_slice = (*xgap)[ indices[ range( 0,xgap->shape()[0] ) ][xgap->shape()[1]-1] ];
        this->score_gap_close(left_edge,&x_slice,max,true);
//cout<<"f2 "<<max->full_score<<endl;
        if(this->first_is_bigger(max->score,best_score) )
        {
            best_score = max->score;
            max->y_ind = xgap->shape()[1]-1;
        }

        align_slice y_slice = (*ygap)[ indices[ygap->shape()[0]-1][ range( 0,ygap->shape()[1] ) ] ];
        this->score_gap_close(right_edge,&y_slice,max,false);
//cout<<"f3 "<<max->full_score<< setprecision (4)<<endl;
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

    if(Settings::noise>5)
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
//if(Settings::noise>6) cout<<"new\n ";

        Edge * edge = site->get_first_fwd_edge();
        int slice_end = (int)g_slice->shape()[0];

        if(edge->get_end_site_index() < slice_end )
        {
            this->score_gap_ext_bwd(edge,g_slice,max_s);
//if(Settings::noise>6) cout<<"1 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_s->bwd_score<<endl;
            this->score_gap_double_bwd(edge,g_slice,max_d);
//if(Settings::noise>6) cout<<"2 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_d->bwd_score<<endl;
            this->score_gap_open_bwd(edge,g_slice,max_m);
//if(Settings::noise>6) cout<<"3 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_m->bwd_score<<endl;
        }
        while(site->has_next_fwd_edge())
        {
            edge = site->get_next_fwd_edge();

            if(edge->get_end_site_index() < slice_end )
            {
                this->score_gap_ext_bwd(edge,g_slice,max_s);
//if(Settings::noise>6) cout<<"4 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_s->bwd_score<<endl;
                this->score_gap_double_bwd(edge,g_slice,max_d);
//if(Settings::noise>6) cout<<"5 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_d->bwd_score<<endl;
                this->score_gap_open_bwd(edge,g_slice,max_m);
//if(Settings::noise>6) cout<<"6 "<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<"; "<<scientific<<max_m->bwd_score<<endl;
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
//if(Settings::noise>6) cout<<"A "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<scientific<<max_m->bwd_score<<endl;
//if(Settings::noise>6) cout<<"A "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<log(max_m->bwd_score)<<endl;
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
//if(Settings::noise>6) cout<<"D "<</*edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<*/"; "<<scientific<<max_m->bwd_score<<endl;
//if(Settings::noise>6) cout<<"D "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<log(max_m->bwd_score)<<endl;
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
//if(Settings::noise>6) cout<<"B "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<scientific<<max_m->bwd_score<<endl;
//if(Settings::noise>6) cout<<"B "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<log(max_m->bwd_score)<<endl;
            }

            while(right_site->has_next_fwd_edge())
            {

                right_edge = right_site->get_next_fwd_edge();

                if(left_edge->get_end_site_index() < left_end
                   && right_edge->get_end_site_index() < right_end )
                {
                    this->score_match_bwd(left_edge,right_edge,max_x,max_y,max_m);
//if(Settings::noise>6) cout<<"C "<</*edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<*/"; "<<scientific<<max_m->bwd_score<<endl;
//if(Settings::noise>6) cout<<"C "<<left_edge->get_start_site_index()<<" "<<left_edge->get_end_site_index()<<"; "<<right_edge->get_start_site_index()<<" "<<right_edge->get_end_site_index()<<"; "<<log(max_m->bwd_score)<<endl;
                }
            }
        }
    }
}


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

//void Simple_alignment::score_gap_close_bwd(Edge *edge,align_slice *z_slice,Matrix_pointer *max)
//{
//    int prev_index = edge->get_end_site_index();
//
//    double this_full_score =  (*z_slice)[prev_index].bwd_score * model->gap_close() * this->get_edge_weight(edge);
//    max->bwd_score += this_full_score;
//}


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
        Site *lsite;
        Site *rsite;


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
