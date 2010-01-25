#include "node.h"
#include <iostream>

using namespace std;
using namespace ppa;

//Node::Node()
//{
//    isleaf = true;
//}

Node::~Node()
{
    if(Settings::noise>4)
        cout<<"deleting node "<<this->get_name()<<endl;

    if(this->has_left_child())
        delete left_child;

    if(this->has_right_child())
        delete right_child;

    if(this->node_has_sequence_object)
        delete sequence;
}

void Node::add_sequence( string seq_string, string full_char_alphabet, bool gapped)
{
    if(Settings::noise>4)
        cout<<"Node::add_sequence "<<name<<"\n";
    sequence = new Sequence(seq_string, full_char_alphabet, gapped);
    this->node_has_sequence_object= true;
}

void Node::get_alignment(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes)
{
    vector<Node*> nodes;
    if(include_internal_nodes)
        this->get_all_nodes(&nodes);
    else
        this->get_leaf_nodes(&nodes);

    for(unsigned int i=0;i<nodes.size();i++)
    {
        Fasta_entry entry;
        entry.name = nodes.at(i)->get_name();
        entry.comment = nodes.at(i)->get_name_comment();

        aligned_sequences->push_back(entry);
    }

//    if(Settings_handle::st.get("sample-additional-paths").as<int>()==0)
    if(true)
    {
        Sequence *root = this->get_sequence();
        int root_length = root->sites_length();
//cout<<"print root\n";
//root->print_sequence();
        for(int j=1;j<root_length-1;j++)
        {
            vector<char> column;
            this->get_alignment_column_at(j,&column,include_internal_nodes);

            for(unsigned int i=0;i<aligned_sequences->size();i++)
            {
                aligned_sequences->at(i).sequence.push_back(column.at(i));
            }
        }
    }
    else
    {
        Sequence *root = this->get_sequence();
        int root_length = root->sites_length();

        cout<<"root length: "<<root_length<<endl;

        vector<int> sampled_sites;

        int site_index=0;
        sampled_sites.push_back(site_index);

        int l_index = 0;
        int r_index = 0;

        while(site_index < root_length-1)
        {
            cout<<"current: "<<site_index<<" ("<<l_index<<","<<r_index<<")"<<endl;

            Site *tsite = root->get_site_at( site_index );

            vector<int> next_sites;

            Edge *fwd_edge = tsite->get_first_fwd_edge();

            Site *nsite = root->get_site_at( fwd_edge->get_end_site_index() );
cout<<"next? "<<nsite->children.left_index<<" "<<nsite->children.right_index<<endl;
            if( ( l_index+1 == nsite->children.left_index && r_index+1 == nsite->children.right_index ) ||
                ( l_index+1 == nsite->children.left_index && nsite->children.right_index == -1 ) ||
                ( nsite->children.left_index == -1 && r_index+1 == nsite->children.right_index ) )
            {
                next_sites.push_back(fwd_edge->get_end_site_index());
            }

            while(tsite->has_next_fwd_edge())
            {
                fwd_edge = tsite->get_next_fwd_edge();

                nsite = root->get_site_at( fwd_edge->get_end_site_index() );
cout<<"next? "<<nsite->children.left_index<<" "<<nsite->children.right_index<<endl;

                if( ( l_index+1 == nsite->children.left_index && r_index+1 == nsite->children.right_index ) ||
                    ( l_index+1 == nsite->children.left_index && nsite->children.right_index == -1 ) ||
                    ( nsite->children.left_index == -1 && r_index+1 == nsite->children.right_index ) )
                {
                    next_sites.push_back(fwd_edge->get_end_site_index());
                }
            }

            cout<<"choices: ";
            for(int i=0;i<(int)next_sites.size();i++)
                cout<<next_sites.at(i)<<" ";
            cout<<endl;
            site_index = (int) ( next_sites.size() * (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
            site_index = next_sites.at(site_index);
            sampled_sites.push_back( site_index );
cout<<"site_index: "<<site_index<<endl;

            nsite = root->get_site_at( site_index );

            if(nsite->children.left_index > 0)
                l_index = nsite->children.left_index;
            if(nsite->children.right_index > 0)
                r_index = nsite->children.right_index;

        }

        cout<<"path sampled\n";
        for(int k=1;k<(int)sampled_sites.size()-1;k++)
        {
            int j = sampled_sites.at(k);
cout<<"picking site "<<j<<endl;
            vector<char> column;
            this->get_alignment_column_at(j,&column,include_internal_nodes);

            for(unsigned int i=0;i<aligned_sequences->size();i++)
            {
                aligned_sequences->at(i).sequence.push_back(column.at(i));
            }
        }

    }
}

void Node::get_alignment_column_at(int j,vector<char> *column,bool include_internal_nodes)
{
    if(leaf)
    {
        int state = sequence->get_site_at(j)->get_state();
        column->push_back(sequence->get_full_alphabet().at(state));
    }
    else
    {
        Site_children *offspring = sequence->get_site_at(j)->get_children();
        int lj = offspring->left_index;
        if(lj>=0)
        {
            left_child->get_alignment_column_at(lj,column,include_internal_nodes);
        }
        else
        {
            int nl = left_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = left_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back('-');
        }

        if(include_internal_nodes)
        {
            int cstate = sequence->get_site_at(j)->get_state();
            char c = sequence->get_full_alphabet().at(cstate);

            int pstate = sequence->get_site_at(j)->get_path_state();
            if( pstate == Site::xskipped || pstate == Site::yskipped )
                c = '-';
            column->push_back(c);
        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            right_child->get_alignment_column_at(rj,column,include_internal_nodes);
        }
        else
        {
            int nl = right_child->get_number_of_leaves();
            if(include_internal_nodes)
                nl = right_child->get_number_of_nodes();

            for(int i=0;i<nl;i++)
                column->push_back('-');
        }
    }
}


void Node::write_metapost_sequence_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath c[];\ndefaultscale := 0.5;\n";
    vector<Site> *sites = this->sequence->get_sites();
    string full_alphabet = this->sequence->get_full_alphabet();

    stringstream all_chars;
    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

        char c = 's';
        if(tsite->get_site_type()==Site::real_site)
            c = full_alphabet.at(tsite->get_state());
        else if(tsite->get_site_type()==Site::stop_site)
            c = 'e';

        string color = this->get_node_fill_color(c);
        if(tsite->get_branch_count_since_last_used()>0)
            color = "0.5white";

        *output<<"c"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<c<<"\","<<color<<");\n";
    }

    if(leaf)
        *output<<"label.top(btex $"<<this->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = this->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {
        Site *tsite =  &sites->at(i);

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
                angle = 40;
            else if(start+3==stop)
                angle = 30;
            else if(start+4<=stop)
                angle = 20;

            stringstream label;
            if(tedge->get_branch_count_since_last_used()>0)
                label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

            *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

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
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                label.str("");

                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(c"<<start<<",c"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

            }
        }
    }

    *output<<"endfig;\n";


    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";

    ++*count;

}

void Node::write_metapost_alignment_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
{
    vector<Site> *sites = this->sequence->get_sites();

    vector<int> left_child_index;
    vector<int> right_child_index;

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            left_child_index.push_back(i);
        }
        if(offspring->right_index>=0)
        {
            right_child_index.push_back(i);
        }
    }


    *output<<"beginfig("<<*count<<");\npickup pencircle scaled 1pt;\npath l[]; path r[];\ndefaultscale := 0.5;\n";
    string full_alphabet = this->sequence->get_full_alphabet();

    *output<<"l0 = circle((0cm,1.5cm),\"s\",white);\n";
    *output<<"r0 = circle((0cm,0cm),\"s\",white);\n";


    stringstream all_chars;
    for(unsigned int i=1;i<sites->size();i++)
    {
        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *lsite = this->left_child->get_sequence()->get_site_at(offspring->left_index);

            char lc = 's';
            if(lsite->get_site_type()==Site::real_site)
                lc = full_alphabet.at(lsite->get_state());
            else if(lsite->get_site_type()==Site::stop_site)
                lc = 'e';

            string color = this->get_node_fill_color(lc);
            if(lsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"l"<<i<<" = circle"<<"(("<<0.5*i<<"cm,1.5cm),\""<<lc<<"\","<<color<<");\n";
        }
        if(offspring->right_index>=0)
        {
            Site *rsite = this->right_child->get_sequence()->get_site_at(offspring->right_index);

            char rc = 's';
            if(rsite->get_site_type()==Site::real_site)
                rc = full_alphabet.at(rsite->get_state());
            else if(rsite->get_site_type()==Site::stop_site)
                rc = 'e';

            string color = this->get_node_fill_color(rc);
            if(rsite->get_branch_count_since_last_used()>0)
                color = "0.5white";

            *output<<"r"<<i<<" = circle"<<"(("<<0.5*i<<"cm,0cm),\""<<rc<<"\","<<color<<");\n";
        }
    }


    if(left_child->is_leaf())
        *output<<"label.top(btex $"<<left_child->get_name()<<"$ etex,(0.125cm,1.75cm));\n";
    else
    {
        string n = left_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,1.75cm));\n";
    }

    if(right_child->is_leaf())
        *output<<"label.top(btex $"<<right_child->get_name()<<"$ etex,(0.125cm,0.25cm));\n";
    else
    {
        string n = right_child->get_name();
        n = n.substr(1,n.length()-2);
        *output<<"label.top(btex \\#"<<n<<"\\# etex,(0.125cm,0.25cm));\n";
    }

    *output<<"defaultscale := 0.25;\n";

    for(unsigned int i=1;i<sites->size();i++)
    {

        Site_children *offspring = sites->at(i).get_children();

        if(offspring->left_index>=0)
        {
            Site *tsite = left_child->get_sequence()->get_site_at(offspring->left_index);

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
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = left_child_index.at( start );
                stop  = left_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<" "<<tedge->get_branch_count_as_skipped_edge()<<" "<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

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
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = left_child_index.at( start );
                    stop  = left_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(l"<<start<<",l"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }

        if(offspring->right_index>=0)
        {
            Site *tsite = right_child->get_sequence()->get_site_at(offspring->right_index);

            if(tsite->has_bwd_edge())
            {
                Edge *tedge = tsite->get_first_bwd_edge();
                int start = tedge->get_start_site_index();
                int stop  = tedge->get_end_site_index();

//                int angle = 0;
//                string place = "edgebot";
//                if(start+1==stop)
//                    place = "edgetop";
//                else if(start+2==stop)
//                    angle = 320;
//                else if(start+3==stop)
//                    angle = 330;
//                else if(start+4<=stop)
//                    angle = 340;
                int angle = 0;
                string place = "edgetop";
                if(start+1==stop)
                    place = "edgebot";
                else if(start+2==stop)
                    angle = 40;
                else if(start+3==stop)
                    angle = 30;
                else if(start+4<=stop)
                    angle = 20;

                start = right_child_index.at( start );
                stop  = right_child_index.at( stop );

                stringstream label;
                if(tedge->get_branch_count_since_last_used()>0)
                    label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                while(tsite->has_next_bwd_edge())
                {
                    tedge = tsite->get_next_bwd_edge();
                    start = tedge->get_start_site_index();
                    stop  = tedge->get_end_site_index();

//                    angle = 0;
//                    place = "edgebot";
//                    if(start+1==stop)
//                        place = "edgetop";
//                    else if(start+2==stop)
//                        angle = 320;
//                    else if(start+3==stop)
//                        angle = 330;
//                    else if(start+4<=stop)
//                        angle = 340;
                    angle = 0;
                    place = "edgetop";
                    if(start+1==stop)
                        place = "edgebot";
                    else if(start+2==stop)
                        angle = 40;
                    else if(start+3==stop)
                        angle = 30;
                    else if(start+4<=stop)
                        angle = 20;

                    start = right_child_index.at( start );
                    stop  = right_child_index.at( stop );

                    label.str("");

                    if(tedge->get_branch_count_since_last_used()>0)
                        label<<"["<<tedge->get_branch_count_since_last_used()<<"~"<<tedge->get_branch_count_as_skipped_edge()<<"~"<<tedge->get_branch_distance_since_last_used()<<"]";

                    *output<<place<<"(r"<<start<<",r"<<stop<<","<<angle<<",\""<<label.str()<<"\",0.5);\n";

                }
            }
        }
    }

    *output<<"endfig;\n";

    string file = Settings_handle::st.get("mpost-graph-file").as<string>();
    float width = (float)this->sequence->get_sites()->size()/(float)root_length;

    *output2<<"\\includegraphics[width="<<width<<"\\columnwidth]{"<<file<<"."<<*count<<"}\n\n\\bigskip\n";
    *output2<<"~\n\n\\bigskip\n";

    ++*count;
}

void Node::check_valid_graph() const
{
    vector<Site> *sites = sequence->get_sites();

    for(unsigned int i=0;i<sites->size();i++)
    {
        Site *ssite = &sites->at(i);
        if( ssite->has_fwd_edge() )
        {
            Edge *edge = ssite->get_first_fwd_edge();
            Site *esite = &sites->at(edge->get_end_site_index());

            if(!esite->contains_bwd_edge(edge,true))
            {
                cout<<"site "<<i<<" has fwd edge from "<<edge->get_start_site_index()<<" to "
                        <<edge->get_end_site_index()<<" but no return\n";
            }

            while( ssite->has_next_fwd_edge() )
            {
                edge = ssite->get_next_fwd_edge();
                esite = &sites->at(edge->get_end_site_index());

                if(!esite->contains_bwd_edge(edge,true))
                {
                    cout<<"site "<<i<<" has fwd edge from "<<edge->get_start_site_index()<<" to "
                            <<edge->get_end_site_index()<<" but no return\n";
                }
            }
        }

        if( ssite->has_bwd_edge() )
        {
            Edge *edge = ssite->get_first_bwd_edge();
            Site *esite = &sites->at(edge->get_start_site_index());

            if(!esite->contains_fwd_edge(edge,true))
            {
                cout<<"site "<<i<<" has bwd edge from "<<edge->get_start_site_index()<<" to "
                        <<edge->get_end_site_index()<<" but no return\n";
            }

            while( ssite->has_next_bwd_edge() )
            {
                edge = ssite->get_next_bwd_edge();
                esite = &sites->at(edge->get_start_site_index());

                if(!esite->contains_fwd_edge(edge,true))
                {
                    cout<<"site "<<i<<" has bwd edge from "<<edge->get_start_site_index()<<" to "
                            <<edge->get_end_site_index()<<" but no return\n";
                }
            }
        }

    }
}


void Node::prune_down()
{
//    cout<<"prune down in "<<this->get_name()<<endl;

    if(this->is_leaf())
        return;

    left_child->prune_down();
    right_child->prune_down();

    if(!left_child->has_sequence()){
        delete this->left_child;
        this->has_left_child(false);
    }

    if(!right_child->has_sequence()){
        delete this->right_child;
        this->has_right_child(false);
    }

    if(this->has_left_child() && !left_child->is_leaf())
    {
        if(!left_child->has_left_child() && left_child->has_right_child())
        {
            Node *new_child = left_child->right_child;
            new_child->set_distance_to_parent (left_child->get_distance_to_parent()+
                                        left_child->right_child->get_distance_to_parent());

            left_child->has_right_child(false);
            delete left_child;
            this->add_left_child(new_child);
        }
        else if(left_child->has_left_child() && !left_child->has_right_child())
        {
            Node *new_child = left_child->left_child;
            new_child->set_distance_to_parent (left_child->get_distance_to_parent()+
                                left_child->left_child->get_distance_to_parent());

            left_child->has_left_child(false);
            delete left_child;
            this->add_left_child(new_child);
        }
    }

    if(this->has_right_child() && !right_child->is_leaf())
    {
        if(!right_child->has_left_child() && right_child->has_right_child())
        {
            Node *new_child = right_child->right_child;
            new_child->set_distance_to_parent (right_child->get_distance_to_parent()+
                                        right_child->right_child->get_distance_to_parent() );

            right_child->has_right_child(false);
            delete right_child;
            this->add_right_child(new_child);
        }
        else if(right_child->has_left_child() && !right_child->has_right_child())
        {
            Node *new_child = right_child->left_child;
            new_child->set_distance_to_parent (right_child->get_distance_to_parent()+
                                            right_child->left_child->get_distance_to_parent());
            right_child->has_left_child(false);
            delete right_child;
            this->add_right_child(new_child);
        }
    }

    if(this->has_left_child() && left_child->has_sequence())
        this->has_sequence(true);

    if(this->has_right_child() && right_child->has_sequence())
        this->has_sequence(true);

//    cout<<"prune down out "<<this->get_name()<<endl;
}

void Node::prune_up()
{
//    cout<<"prune up in "<<this->get_name()<<endl;

    if(!this->is_leaf() && !this->has_left_child() && this->has_right_child())
    {
        Node* tmp_child = right_child;
        left_child = tmp_child->left_child;
        right_child = tmp_child->right_child;
        delete tmp_child;
    }

    if(!this->is_leaf() && this->has_left_child() && !this->has_right_child())
    {
        Node* tmp_child = left_child;
        left_child = tmp_child->left_child;
        right_child = tmp_child->right_child;
        delete tmp_child;
    }
//    cout<<"prune up out "<<this->get_name()<<endl;
}
