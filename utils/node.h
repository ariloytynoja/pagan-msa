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

class Node
{
    Node *left_child;
    Node *right_child;
    bool leaf;
    double dist_to_parent;
    string name;
    string name_comment;

    Sequence *sequence;
public:
    Node() : leaf(true), dist_to_parent(0), name(""){}
    ~Node();

    /**************************************/

    void set_distance_to_parent(double d)
    {
        dist_to_parent = d;
        if(Settings::noise>5) cout<<"node.set_distance_to_parent("<<d<<")\n";
    }

    double get_distance_to_parent() { return dist_to_parent; }

    /**************************************/

    void set_name(string nm)
    {
        name = nm;
        if(Settings::noise>5) cout<<"node.set_name("<<name<<")\n";
    }

    string get_name() const { return name; }

    /**************************************/

    void add_name_comment( string comment ) { name_comment = comment; }

    string get_name_comment() { return name_comment; }

    bool is_leaf() { return leaf; }

    /**************************************/

    void add_left_child(Node *child)
    {
        left_child = child;
        leaf = false;
        if(Settings::noise>5) cout<<"node.add_left_child(Node *)\n";
    }

    void add_right_child(Node *child)
    {
        right_child = child; leaf = false;
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
        if(!leaf)
        {
            left_child->get_leaf_nodes(nodes);
            right_child->get_leaf_nodes(nodes);
        }
        else
            nodes->push_back(this);
    }

    void get_internal_nodes(vector<Node*> *nodes)
    {
        if(!leaf)
        {
            left_child->get_internal_nodes(nodes);
            nodes->push_back(this);
            right_child->get_internal_nodes(nodes);
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

    /************************************/

    void align_sequences(Model_factory *mf)
    {
        if(leaf)
            return;

        left_child->align_sequences(mf);
        right_child->align_sequences(mf);

        this->align_sequences_this_node(mf);
    }

    void align_sequences_this_node(Model_factory *mf)
    {

        if(Settings::noise>0)
            cout<<"aligning node "<<this->get_name()<<": "<<left_child->get_name()<<" - "<<right_child->get_name()<<"."<<endl;

        int oldv = Settings::noise;
//        if(name=="#31#")
//            Settings::noise = 7;
        double dist = left_child->get_distance_to_parent()+right_child->get_distance_to_parent();
        Dna_model model = mf->dna_alignment_model(dist);

        Simple_alignment sa;
        sa.align(left_child->get_sequence(),right_child->get_sequence(),&model,left_child->get_distance_to_parent(),right_child->get_distance_to_parent());

        Settings::noise = oldv;

        if(Settings::noise>1)
            cout<<" finished...\n";

        this->add_ancestral_sequence( sa.get_simple_sequence() );

        if(Settings::noise>1)
            cout<<" leaving!\n";

        if( Settings_handle::st.is("check-valid-graphs") )
            this->check_valid_graph();

    }

    /************************************/

    void get_alignment(vector<Fasta_entry> *aligned_sequences,bool include_internal_nodes=false);

    void get_alignment_column_at(int j,vector<char> *column,bool include_internal_nodes);

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

    string print_tree() {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree()<<","<<right_child->print_subtree()<<");";
            return ss.str();
        } else {
            return "";
        }
    }

    /************************************/

    string print_subtree() {
        if(!leaf)
        {
            stringstream ss;
            ss<<"("<<left_child->print_subtree()<<","<<right_child->print_subtree()<<"):"<<dist_to_parent;
            return ss.str();
        } else {
            stringstream ss;
            ss<<name<<":"<<dist_to_parent;
            return ss.str();
        }
    }

    void write_sequence_graphs(bool overwrite=true) const throw (Exception)
    {

        string file = Settings_handle::st.get("mpost-graphfile").as<string>();
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

    string get_node_fill_color(char c) const
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

    void add_sequence( string seq_string, string full_dna_alphabet);

    void add_ancestral_sequence( Sequence* s ) { sequence = s; }

    Sequence *get_sequence() { return sequence; }

    void check_valid_graph() const;

};

}
#endif // NODE_H
