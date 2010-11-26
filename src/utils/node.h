#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "exceptions.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "main/sequence.h"
#include "main/simple_alignment.h"
#include "utils/model_factory.h"
#include "utils/fasta_entry.h"

extern bool DEBUG;

using namespace std;

namespace ppa
{

struct Insertion_at_node
{
    string node_name_wanted;
    int length;
    bool left_child_wanted;
};

class Node
{
    Node *left_child;
    Node *right_child;
    bool leaf;
    double dist_to_parent;
    string name;
    string name_comment;
    string name_id;
    string nhx_tid;

    Sequence *sequence;
    bool node_has_sequence;

    bool node_has_sequence_object;
    bool node_has_left_child;
    bool node_has_right_child;
    bool adjust_left_node_site_index;
    bool adjust_right_node_site_index;

public:
    Node() : leaf(true), dist_to_parent(0), name("undefined"), nhx_tid(""),
                    node_has_sequence(false), node_has_sequence_object(false),
                    node_has_left_child(false), node_has_right_child(false),
                    adjust_left_node_site_index(false),
                    adjust_right_node_site_index(false){}
    ~Node();

    /**************************************/

    void set_distance_to_parent(double d)
    {
        dist_to_parent = d;

        if( Settings_handle::st.is("scale-branches") || Settings_handle::st.is("truncate-branches")  )
        {
            if( Settings_handle::st.is("scale-branches") &&
                  Settings_handle::st.get("scale-branches").as<float>() > 0 )
            {
                dist_to_parent *= Settings_handle::st.get("scale-branches").as<float>();
            }

            if( Settings_handle::st.is("real-branches") )
                ;
            else if( Settings_handle::st.is("truncate-branches") &&
                  Settings_handle::st.get("truncate-branches").as<float>() > 0 &&
                    dist_to_parent > Settings_handle::st.get("truncate-branches").as<float>() )
            {
                dist_to_parent = Settings_handle::st.get("truncate-branches").as<float>();
            }
        }
        else if( Settings_handle::st.is("fixed-branches") )
        {
            dist_to_parent = Settings_handle::st.get("fixed-branches").as<float>();
        }

        if(Settings::noise>5) cout<<"node.set_distance_to_parent("<<d<<")\n";
    }

    double get_distance_to_parent() { return dist_to_parent; }

    /**************************************/

    void set_name(string nm)
    {
        name = nm;
        if(Settings::noise>5) cout<<"node.set_name("<<name<<")\n";
    }

    void set_name_ids(int *count)
    {
        if(!leaf)
        {
            left_child->set_name_ids(count);
            right_child->set_name_ids(count);
        } else {
            stringstream ss;
            ss << "seq" << (++*count);
            this->name_id = ss.str();
        }
    }

    void set_name_id(int count)
    {
        stringstream ss;
        ss << "seq" << (++count);
        this->name_id = ss.str();
    }

    void set_nhx_tid(string nm)
    {
        nhx_tid = nm;
    }

    string get_name() const { return name; }

    string get_name_id() const { return name_id; }

    string get_nhx_tid() const { return nhx_tid; }

    string get_id_for_name(string query) const
    {
        if(this->name == query)
        {
            return this->name_id;
        }
        else if(!leaf)
        {
            string tmp;
            tmp = left_child->get_id_for_name(query);
            if(tmp!="")
                return tmp;

            tmp = right_child->get_id_for_name(query);
            if(tmp!="")
                return tmp;

            return "";
        }

        return "";
    }

    /**************************************/

    void add_name_comment( string comment ) { name_comment = comment; }

    string get_name_comment() { return name_comment; }

    bool is_leaf() { return leaf; }
    void is_leaf(bool i) { leaf = i; }

    bool has_left_child() { return node_has_left_child; }
    void has_left_child(bool h) { node_has_left_child = h; }

    bool has_right_child() { return node_has_right_child; }
    void has_right_child(bool h) { node_has_right_child = h; }

    bool left_needs_correcting_sequence_site_index() { return adjust_left_node_site_index; }
    void left_needs_correcting_sequence_site_index(bool h) { adjust_left_node_site_index = h; }

    bool right_needs_correcting_sequence_site_index() { return adjust_right_node_site_index; }
    void right_needs_correcting_sequence_site_index(bool h) { adjust_right_node_site_index = h; }

    /**************************************/

    void add_left_child(Node *child)
    {
        left_child = child; leaf = false; this->has_left_child(true);
        if(Settings::noise>5) cout<<"node.add_left_child(Node *)\n";
    }

    void add_right_child(Node *child)
    {
        right_child = child; leaf = false; this->has_right_child(true);
        if(Settings::noise>5) cout<<"node.add_right_child(Node *)\n";
    }
    Node * get_left_child() { return left_child; }

    Node * get_right_child() { return right_child; }

    /**************************************/


    void get_all_nodes(vector<Node*> *nodes)
    {
        if(!leaf)
            left_child->get_all_nodes(nodes);

        nodes->push_back(this);

        if(!leaf)
            right_child->get_all_nodes(nodes);
    }

    void get_leaf_nodes(vector<Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_leaf_nodes(nodes);
            right_child->get_leaf_nodes(nodes);
        }
        else
            nodes->push_back(this);
    }

    void get_leaf_nodes(map<string,Node*> *nodes)
    {
        if(this->is_leaf())
        {
            nodes->insert(pair<string,Node*>(this->get_name(),this));
        }
        else
        {
            left_child->get_leaf_nodes(nodes);
            right_child->get_leaf_nodes(nodes);
        }
    }

    void get_internal_nodes(vector<Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_nodes(nodes);
            nodes->push_back(this);
            right_child->get_internal_nodes(nodes);
        }
    }

    void get_internal_nodes(map<string,Node*> *nodes)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_nodes(nodes);
            nodes->insert(pair<string,Node*>(this->get_name(),this));
            right_child->get_internal_nodes(nodes);
        }
    }

    void get_internal_node_names(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_node_names(list);
            right_child->get_internal_node_names(list);

            list->insert(pair<string,string>(this->get_name(),this->get_name()));
//            list->insert(pair<string,string>(this->get_nhx_tid(),this->get_name()));
        }
    }

    void get_internal_node_names_with_tid_tag(multimap<string,string> *list)
    {
        if(!this->is_leaf())
        {
            left_child->get_internal_node_names_with_tid_tag(list);
            right_child->get_internal_node_names_with_tid_tag(list);

            if(this->get_nhx_tid()!="")
            {
                list->insert(pair<string,string>(this->get_nhx_tid(),this->get_name()));
            }
        }
    }

    void name_internal_nodes(int *count)
    {
        if(leaf)
            return;
        else
        {
            left_child->name_internal_nodes(count);
            right_child->name_internal_nodes(count);

            stringstream ss;
            ss<<"#"<<(*count)<<"#";

            this->set_name(ss.str());
            (*count)++;
        }

    }

    bool sequence_site_index_needs_correcting()
    {
        if(this->is_leaf())
            return false;
        else if(adjust_left_node_site_index || adjust_right_node_site_index)
            return true;
        else if(left_child->sequence_site_index_needs_correcting())
            return true;
        else if(right_child->sequence_site_index_needs_correcting())
            return true;
        else
            return false;
    }

    /************************************/

    void has_sequence(bool s) { node_has_sequence = s; }
    bool has_sequence() { return node_has_sequence; }

    void prune_tree() { this->prune_down(); this->prune_up(); }
    void prune_up();
    void prune_down();

    /************************************/

    /*
     * Note on graphical output: posterior plots have to written *during* the alignment
     * as the dynamic programming matrices aren't stored for later use.
     * In contrast, sequence graphs exist till the very end.
     */
    void start_alignment(Model_factory *mf)
    {
        if( Settings_handle::st.is("seqfile") )
        {
            if(Settings_handle::st.is("mpost-posterior-plot-file"))
            {
                this->start_mpost_plot_file();
            }

            this->align_sequences(mf);

            if(Settings_handle::st.is("mpost-posterior-plot-file"))
            {
                this->finish_mpost_plot_file();
            }

            if(Settings_handle::st.is("output-ancestors"))
            {
                this->reconstruct_parsimony_ancestor(mf);
            }
        }
        else if( Settings_handle::st.is("ref-seqfile") )
        {
            this->read_alignment(mf);

            this->reconstruct_parsimony_ancestor(mf);

        }
    }

    void align_sequences(Model_factory *mf)
    {
        if(leaf)
            return;

        left_child->align_sequences(mf);
        right_child->align_sequences(mf);

        this->align_sequences_this_node(mf);

    }

    void align_sequences_this_node(Model_factory *mf, bool is_reads_sequence=false, bool is_local_alignment=false)
    {

        if(Settings::noise>0)
            cout<<"aligning node "<<this->get_name()<<": "<<left_child->get_name()<<" - "<<right_child->get_name()<<"."<<endl;

        double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();
        Evol_model model = mf->alignment_model(dist,is_local_alignment);

        Simple_alignment sa;
        sa.align(left_child->get_sequence(),right_child->get_sequence(),&model,
                 left_child->get_distance_to_parent(),right_child->get_distance_to_parent(), is_reads_sequence);

        this->add_ancestral_sequence( sa.get_simple_sequence() );

        if( Settings_handle::st.is("check-valid-graphs") )
            this->check_valid_graph();
    }

    void read_alignment(Model_factory *mf)
    {
        if(leaf)
            return;

        left_child->read_alignment(mf);
        right_child->read_alignment(mf);

        this->read_alignment_this_node(mf);
    }

    void read_alignment_this_node(Model_factory *mf)
    {

        if(Settings::noise>0)
            cout<<"reading alignment node "<<this->get_name()<<": "<<left_child->get_name()<<" - "<<right_child->get_name()<<"."<<endl;

        double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();
        Evol_model model = mf->alignment_model(dist);

        Simple_alignment sa;
        sa.read_alignment(left_child->get_sequence(),right_child->get_sequence(),&model,
                 left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

        this->add_ancestral_sequence( sa.get_simple_sequence() );

    }

    void reconstruct_parsimony_ancestor(Model_factory *mf)
    {
        for(int i=1;i<this->get_sequence()->sites_length()-1;i++)
        {
            int state = this->get_sequence()->get_site_at(i)->get_state();
            this->reconstruct_parsimony_ancestor_at_site(mf,i,state,false);
        }
    }

    void reconstruct_parsimony_ancestor_at_site(Model_factory *mf,int pos,int parent_state,bool is_matched)
    {

        if(leaf)
            return;

        Site *site = this->get_sequence()->get_site_at(pos);

        int pstate = site->get_path_state();
//        cout<<pos<<": "<<this->get_name()<<" "<<parent_state<<" "<<site->get_state()<<" - ("<<pstate<<") "<<is_matched<<endl;

        if( pstate == Site::matched )
        {
            int new_state = mf->get_child_parsimony_state(parent_state,site->get_state());

            site->set_state(new_state);

            is_matched = true;

//            cout<<pos<<": "<<this->get_name()<<" "<<parent_state<<" "<<site->get_state()<<" "<<new_state<<" ("<<pstate<<") "<<is_matched<<endl;
        }

        if( !is_matched )
        {
            site->set_site_type(Site::non_real);
        }

        if(site->children.left_index >= 0)
        {
            this->left_child->reconstruct_parsimony_ancestor_at_site(mf,site->children.left_index,site->get_state(),is_matched);
        }
        if(site->children.right_index >= 0)
        {
            this->right_child->reconstruct_parsimony_ancestor_at_site(mf,site->children.right_index,site->get_state(),is_matched);
        }

    }

    bool has_site_at_alignment_column(int j,string node_name)
    {

        if(this->get_name() == node_name)
        {
            return true;
        }
        else if(leaf)
        {
            return false;
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(j)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
            {
                bool l = left_child->has_site_at_alignment_column(lj,node_name);
                if(l)
                    return true;
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                bool r = right_child->has_site_at_alignment_column(rj,node_name);
                if(r)
                    return true;
            }
        }
        return false;
    }

    int get_state_at_alignment_column(int j,string node_name)
    {

        if(this->get_name() == node_name)
        {
            return this->get_sequence()->get_site_at(j)->get_state();
        }
        else
        {
            Site_children *offspring = sequence->get_site_at(j)->get_children();
            int lj = offspring->left_index;
            if(lj>=0)
            {
                int l = left_child->get_state_at_alignment_column(lj,node_name);
                if(l>=0)
                    return l;
            }

            int rj = offspring->right_index;
            if(rj>=0)
            {
                int r = right_child->get_state_at_alignment_column(rj,node_name);
                if(r>=0)
                    return r;
            }
        }
        return -1;
    }

    bool has_additional_sites_before_alignment_column(int j)
    {

        Site_children *offspring = sequence->get_site_at(j)->get_children();

        int lj = offspring->left_index;
        if(lj>=0)
        {
            bool l = left_child->has_additional_sites_before_alignment_column(lj);
            if(l)
                return true;
        }

        int rj = offspring->right_index;
        if(rj>=0)
        {
            bool r = right_child->has_additional_sites_before_alignment_column(rj);
            if(r)
                return true;
        }

        if(!(this->adjust_left_node_site_index && this->adjust_right_node_site_index))
        {
            return false;
        }
        else
        {
            int prev_lj = -1;
            int prev_rj = -1;

            if(j>0)
            {
                Site_children *prev_offspring = sequence->get_site_at(j-1)->get_children();
                prev_lj = prev_offspring->left_index;
                prev_rj = prev_offspring->right_index;
            }

            if(this->adjust_left_node_site_index && lj-prev_lj!=1)
                return true;
            else if(this->adjust_right_node_site_index && rj-prev_rj!=1)
                return true;
            else
                return false;
        }
    }

    void additional_sites_before_alignment_column(int j,vector<Insertion_at_node> *addition);
//    {
//        if(this->is_leaf())
//            return;
//
//        Site_children *offspring = sequence->get_site_at(j)->get_children();
//
//        int lj = offspring->left_index;
//        int rj = offspring->right_index;
//
//        if(lj>=0)
//            left_child->additional_sites_before_alignment_column(lj,addition);
//
//        if(j>0)
//        {
//            Site_children *prev_offspring = sequence->get_site_at(j-1)->get_children();
//
//            int prev_lj = prev_offspring->left_index;
//            int prev_rj = prev_offspring->right_index;
//
//            if(lj>0 && prev_lj>=0 && lj-prev_lj != 1)
//            {
//                Insertion_at_node ins;
//                ins.node_name_wanted = this->get_name();
//                ins.length = lj-prev_lj-1;
//                ins.left_child_wanted = true;
//                addition->push_back(ins);
//            }
//
//            if(rj>0 && prev_rj>=0 && rj-prev_rj != 1)
//            {
//                Insertion_at_node ins;
//                ins.node_name_wanted = this->get_name();
//                ins.length = rj-prev_rj-1;
//                ins.left_child_wanted = false;
//                addition->push_back(ins);
//            }
//        }
//
//
//        if(rj>=0)
//            right_child->additional_sites_before_alignment_column(rj,addition);
//
//    }


    void start_mpost_plot_file()
    {
        string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), (ios::out) );
        if (! output) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        ofstream output2(path2.c_str(), (ios::out) );
        if (! output2) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }


        output<<"vardef circle(expr a,c,col) =\n path p;\n h := 0.6mm;\n"
                " p = a+(-h/2,0){up}..{right}a+(0,h/2){right}..{down}a+(h/2,0){down}..{left}a+(0,-h/2){left}..{up}cycle;\n"
                " fill p withcolor col;\n draw(p);\n pair x;\n x = ((point(0) of p)+(point(2) of p))/2;\n"
                " label(c,x);\n p\nenddef;\n\n";

        output<<"vardef edge(expr px,ax,py,ay,s) =\n pair x,y;\n x = ((point(0) of px)+(point(2) of px))/2;\n"
                " y = ((point(0) of py)+(point(2) of py))/2;\n path p; p = (x){dir(ax)}..{dir(ay)}(y);\n"
                " draw (p) cutbefore px cutafter py withpen pencircle scaled (0.3*s);\n point .5*length p of p\n"
                "enddef;\n\n";

        output<<"def edgetop(expr px,py,a,c,s) =\n label.top(c,edge(px,a,py,-1*a,s));\nenddef;\n"
                "def edgebot(expr px,py,a,c,s) =\n label.bot(c,edge(px,360-a,py,360+1*a,s));\nenddef;\n\n";
        output<<"def edgelft(expr px,py,a,c,s) =\n label.lft(c,edge(px,a,py,180-1*a,s));\nenddef;\n"
                "def edgergt(expr px,py,a,c,s) =\n label.rt(c,edge(px,180-a,py,a,s));\nenddef;\n\n";

        output2<<"\\documentclass{article}\n\\usepackage{geometry,graphicx}\n"
                 "\\geometry{a4paper,tmargin=0.5in,bmargin=0.5in,lmargin=.5in,rmargin=0.5in}\n"
                 "\\begin{document}\n";


        output.close();
        output2.close();
    }

    void finish_mpost_plot_file()
    {
        string file = Settings_handle::st.get("mpost-posterior-plot-file").as<string>();
        cout<<"Plot file: "<<file<<endl;

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), (ios::out|ios::app) );
        if (! output) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        ofstream output2(path2.c_str(), (ios::out|ios::app) );
        if (! output2) { throw IOException ("Node::start_mpost_plot_file. Failed to open file"); }

        output<<"end;\n";
        output2<<"\\end{document}\n";

        output.close();
        output2.close();

        cout<<"\nThe posterior probability plot files can generated using following commands:\n"
              "  mpost "<<file<<".mp\n"
              "  latex "<<file<<".tex\n"
              "  dvipdf "<<file<<".dvi\n\n";

    }

    /************************************/

    void get_alignment(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes=false);

    void get_alignment_column_at(int j,vector<char> *column,bool include_internal_nodes);

    void get_multiple_alignment_columns_before(int j,vector< vector<char> > *columns, string node_name_wanted, bool left_child_wanted,bool include_internal_nodes);

    int get_number_of_leaves()
    {
        if(leaf)
             return 1;
        else
            return left_child->get_number_of_leaves()+right_child->get_number_of_leaves();
    }

    int get_number_of_nodes()
    {
        if(leaf)
             return 1;
        else
            return 1+left_child->get_number_of_nodes()+right_child->get_number_of_nodes();
    }

    /************************************/

    void print_node_info()
    {
        if(!leaf)
            left_child->print_node_info();
        cout<<name<<":"<<dist_to_parent<<"\n";
        if(!leaf)
            right_child->print_node_info();
    }

    /************************************/

    string print_tree(bool int_names=false) {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree(int_names)<<","<<right_child->print_subtree(int_names)<<")";
            if(int_names)
                ss<<this->get_name()<<":0";
            ss<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    /************************************/

    string print_subtree(bool int_names=false) {
//        cout<<"print_subtree: "<<get_name()<<endl;
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree(int_names)<<","<<right_child->print_subtree(int_names)<<")";
            if(int_names)
                ss<<this->get_name();
            ss<<":"<<dist_to_parent;
            return ss.str();
        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent;
            return ss.str();
        }
    }

    void write_nhx_tree(string path, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".nhx_tree").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));

        if (! output) { throw IOException ("Node::write_nhx_tree. Failed to open file"); }

        output<<this->print_nhx_tree();
        output.close();
    }

    string print_nhx_tree() const {
        if(!leaf)
        {
            stringstream tid("");
            if(this->get_nhx_tid()!="")
                tid << "[&&NHX:TID="+this->get_nhx_tid()<<"]";

            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree()<<","<<right_child->print_nhx_subtree()<<"):"<<dist_to_parent<<tid.str()<<";";
            return ss.str();
        } else {
            return "";
        }
    }

    /************************************/

    string print_nhx_subtree() const {

        stringstream tid("");
        if(this->get_nhx_tid()!="")
            tid << "[&&NHX:TID="+this->get_nhx_tid()<<"]";

        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_nhx_subtree()<<","<<right_child->print_nhx_subtree()<<"):"<<dist_to_parent<<tid.str();
            return ss.str();

        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent<<tid.str();
            return ss.str();
        }
    }

    string print_xml_tree() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_xml_subtree()<<","<<right_child->print_xml_subtree()<<")"<<name<<":0;";
            return ss.str();
        } else {
            return "";
        }
    }

    /************************************/

    string print_xml_subtree() const {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_xml_subtree()<<","<<right_child->print_xml_subtree()<<")"<<name<<":"<<dist_to_parent;
            return ss.str();
        } else {
            stringstream ss;
            ss<<name_id<<":"<<dist_to_parent;
            return ss.str();
        }
    }


    void write_sequence_graphs(bool overwrite=true) const throw (Exception)
    {

        string file = Settings_handle::st.get("mpost-graph-file").as<string>();
        cout<<"Graph file: "<<file<<endl;

        string path = file;
        path.append(".mp");

        string path2 = file;
        path2.append(".tex");

        ofstream output(path.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        if (! output) { throw IOException ("Node::write_sequence_graphs. Failed to open file"); }

        output<<"vardef circle(expr a,c,col) =\n path p;\n h := 0.2cm;\n"
                " p = a+(-h/2,0){up}..{right}a+(0,h/2){right}..{down}a+(h/2,0){down}..{left}a+(0,-h/2){left}..{up}cycle;\n"
                " fill p withcolor col;\n draw(p);\n pair x;\n x = ((point(0) of p)+(point(2) of p))/2;\n"
                " label(c,x);\n p\nenddef;\n\n";

        output<<"vardef edge(expr px,ax,py,ay,s) =\n pair x,y;\n x = ((point(0) of px)+(point(2) of px))/2;\n"
                " y = ((point(0) of py)+(point(2) of py))/2;\n path p; p = (x){dir(ax)}..{dir(ay)}(y);\n"
                " draw (p) cutbefore px cutafter py withpen pencircle scaled s;\n point .5*length p of p\n"
                "enddef;\n\n";

        output<<"def edgetop(expr px,py,a,c,s) =\n label.top(c,edge(px,a,py,-1*a,s));\nenddef;"
                "def edgebot(expr px,py,a,c,s) =\n label.bot(c,edge(px,360-a,py,360+1*a,s));\nenddef;\n\n";

        ofstream output2(path2.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        if (! output2) { throw IOException ("Node::write_sequence_graphs. Failed to open file"); }

        output2<<"\\documentclass{article}\n\\usepackage{geometry,graphicx}\n"
                "\\geometry{a4paper,tmargin=0.5in,bmargin=0.5in,lmargin=.5in,rmargin=0.5in}\n"
                "\\begin{document}\n";

        int fig_count = 0;
        int root_length = sequence->get_sites()->size();

        this->write_metapost_graphs(&output,&output2,&fig_count,root_length);

        output<<"end;\n";
        output.close();

        output2<<"\\end{document}\n";
        output2.close();


        cout<<"\nThe sequence graph files can generated using following commands:\n"
              "  mpost "<<file<<".mp\n"
              "  latex "<<file<<".tex\n"
              "  dvipdf "<<file<<".dvi\n\n";

    }

    void write_metapost_graphs(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception)
    {
        if(leaf)
        {
            if(Settings_handle::st.is("output-leaf-graphs"))
                this->write_metapost_sequence_graph(output, output2, count, root_length);
        }
        else
        {
            left_child->write_metapost_graphs(output, output2, count, root_length);

            this->write_metapost_sequence_graph(output, output2, count, root_length);

            if(Settings_handle::st.is("output-alignment-graphs"))
                this->write_metapost_alignment_graph(output, output2, count, root_length);

            right_child->write_metapost_graphs(output, output2, count, root_length);
        }
    }

    void write_metapost_sequence_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception);

    void write_metapost_alignment_graph(ostream *output, ostream *output2, int *count, int root_length) const throw (Exception);

    static string get_node_fill_color(char c)
    {
        string color = "white";
        switch(c)
        {
            case 'A': color = "0.5blue"; break;
            case 'C': color = "0.5green"; break;
            case 'G': color = "(0.85,0.85,0)"; break;
            case 'T': color = "0.5red"; break;
            case 'R': color = "(0,0.5,0.5)"; break;
            case 'Y': color = "(0,0.5,0.5)"; break;
            case 'M': color = "(0,0.5,0.5)"; break;
            case 'K': color = "(0,0.5,0.5)"; break;
            case 'W': color = "(0,0.5,0.5)"; break;
            case 'S': color = "(0,0.5,0.5)"; break;
            case 'B': color = "(0.5,0,0.5)"; break;
            case 'D': color = "(0.5,0,0.5)"; break;
            case 'H': color = "(0.5,0,0.5)"; break;
            case 'V': color = "(0.5,0,0.5)"; break;
        }
        return color;
    }

    void add_sequence( Fasta_entry seq_entry, string full_char_alphabet, bool gapped = false, bool no_trimming = false);

    void add_ancestral_sequence( Sequence* s ) { sequence = s;  node_has_sequence_object = true;}

    Sequence *get_sequence() { return sequence; }

    void check_valid_graph() const;

};

}
#endif // NODE_H
