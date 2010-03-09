//
// This file is partly based on the code from the Bio++ library.
// The following copyright information is given.
//

//
// File: Fasta.cpp
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Aug 21 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

Julien.Dutheil@univ-montp2.fr

This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include <algorithm>
#include <iostream>
#include <set>
#include "fasta_reader.h"
#include "text_utils.h"
#include "settings_handle.h"

using namespace ppa;

/****************************************************************************************/

void Fasta_reader::read(istream & input, vector<Fasta_entry> & seqs, bool short_names = false) const throw (Exception)
{
    if (!input) { throw IOException ("Fasta_reader::read. Failed to open file"); }


    char c = input.get();
    while(c==' ' || c=='\n')
    {
        c = input.get();
    }

    if(c=='>')
    {
        input.unget();
        this->read_fasta(input,seqs,short_names);
    }
    else if(c=='@')
    {
        input.unget();
        this->read_fastq(input,seqs);
    }
    else if(c=='#')
    {
        input.unget();
        this->read_graph(input,seqs,short_names);
    }
    else
    {
        cout<<"Input file format unrecognized. Only FASTA and FASTQ formats supported. Exiting.\n\n";
        exit(-1);
    }

    set<string> names;
    for(int i=0;i<(int)seqs.size();i++)
    {
        string name = seqs.at(i).name;
        while(names.find(name) != names.end())
        {
            cout<<"Sequence name "<<name<<" is defined more than once! Adding suffix '.1'.\n";
            seqs.at(i).name += ".1";
            name += ".1";
        }
        names.insert(name);
    }

}

void Fasta_reader::read_fasta(istream & input, vector<Fasta_entry> & seqs, bool short_names = false) const throw (Exception)
{

    string temp, name, comment, sequence = "";  // Initialization

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string


        // If first character is >
        if(temp[0] == '>')
        {
            temp = Text_utils::remove_last_whitespaces(temp);

            // If a name and a sequence were found
            if((name != "") && (sequence != ""))
            {
                Fasta_entry fe;
                fe.name = name;
                fe.comment = comment;
                fe.sequence = sequence;
                seqs.push_back(fe);
                name = "";
                sequence = "";
            }
            // Sequence name isolation
            if (! short_names)
            {
                name = temp;
                comment = "";
            }
            else
            {
                String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
                name = st->next_token();
                comment = "";
                while (st->has_more_token())
                {
                  comment += st->next_token()+" ";
                }
                delete st;
            }
            name.erase(name.begin());  // Character > deletion
        }
        else sequence += temp;  // Sequence isolation
    }

    // Addition of the last sequence in file
    if((name != "") && (sequence != ""))
    {
//        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );

        Fasta_entry fe;
        fe.name = name;
        fe.comment = comment;
        fe.sequence = sequence;
        seqs.push_back(fe);
    }

}

void Fasta_reader::read_fastq(istream & input, vector<Fasta_entry> & seqs) const throw (Exception)
{
    string temp, name, comment, sequence, quality = "";  // Initialization

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string
        // If first character is @
        if(temp[0] == '@')
        {
            temp = Text_utils::remove_last_whitespaces(temp);
            Fasta_entry fe;
            String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
            name = st->next_token();

            name.erase(name.begin());  // Character @ deletion
            fe.name = name;

            comment = "";
            while (st->has_more_token())
            {
                string block = st->next_token();
                block = Text_utils::remove_surrounding_whitespaces(block);
                comment += block+" ";

                if(block.substr(0,4)=="TID=")
                {
                    block = block.substr(4);
                    fe.tid = block;
                }
            }
            delete st;


            getline(input, temp, '\n');  // Copy current line in temporary string
            temp = Text_utils::remove_last_whitespaces(temp);
            sequence = temp;
            fe.sequence = sequence;

            getline(input, temp, '\n');  // Copy current line in temporary string
            temp = Text_utils::remove_last_whitespaces(temp);
            if(temp[0] != '+')
            {
                cout<<"Error in FASTQ comment:"<<temp<<"\nExiting.\n\n";
                exit(-1);
            }

            temp.erase(temp.begin());  // Character + deletion

            if(temp.length()>0)
            {
                comment += " ; ";
                comment += temp;
            }
            fe.comment = comment;

            getline(input, temp, '\n');  // Copy current line in temporary string
            temp = Text_utils::remove_last_whitespaces(temp);
            quality = temp;
            fe.quality = quality;

            seqs.push_back(fe);
        }
        else if(temp != "")
        {
            cout<<"FASTQ file parse error. Expecting a line starting with '@':  \n"<<temp<<endl<<endl<<"Exiting\n\n.";
            exit(-1);
        }
    }
}

void Fasta_reader::read_graph(istream & input, vector<Fasta_entry> & seqs, bool short_names = false) const throw (Exception)
{

    string temp, name, comment, sequence, block = "";  // Initialization
    int prev_site = -1;
    vector<Seq_edge> edges;

    while(!input.eof())
    {
        getline(input, temp, '\n');  // Copy current line in temporary string

        // If first character is >
        if(temp[0] == '#')
        {
            temp = Text_utils::remove_last_whitespaces(temp);

            // If a name and a sequence were found
            if((name != "") && (sequence != ""))
            {
                Fasta_entry fe;
                fe.name = name;
                fe.comment = comment;
                fe.sequence = sequence;
                fe.edges = edges;
                seqs.push_back(fe);
                name = "";
                sequence = "";
            }       
            edges.clear();
            prev_site = -1;

            // Sequence name isolation
            if (! short_names)
            {
                name = temp;
                comment = "";
            }
            else
            {
                String_tokenizer * st = new String_tokenizer(temp, " ", true, false);
                name = st->next_token();
                comment = "";
                while (st->has_more_token())
                {
                  comment += st->next_token()+" ";
                }
                delete st;
            }
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            temp = Text_utils::remove_surrounding_whitespaces(temp);
            if(temp.size()==0)
                continue;

//            cout<<temp<<endl;
            String_tokenizer * st = new String_tokenizer(temp, ";", true, false);
            block = st->next_token();

            String_tokenizer * bt = new String_tokenizer(block, " ", true, false);

            temp = bt->next_token();
            int site = Text_utils::to_int( temp );
            if(site != prev_site+1)
            {
                if(prev_site == -2)
                {
                    cout<<"Error reading the graph input: 'end' site is not the last site.\nExiting.\n\n";
                }
                else
                {
                    cout<<"Error reading the graph input: previous site "<<prev_site<<" and this site "<<site<<".\nExiting.\n\n";
                }
                exit(-1);
            }
            prev_site++;

            temp = bt->next_token();
            temp = Text_utils::remove_surrounding_whitespaces(temp);
            if(temp == "start")
            {
                if(site != 0)
                {
                    cout<<"Error reading the graph input: 'start' is not the site 0.\nExiting.\n\n";
                    exit(-1);
                }
            }
            else if(temp == "end")
            {
                prev_site = -2;
            }
            else
            {
                sequence += temp[0];
            }

            delete bt;

            // Next block
            double sum_weight = 0;

            block = st->next_token();
            block = Text_utils::remove_surrounding_whitespaces(block);

            while(block != "")
            {

                String_tokenizer * bt = new String_tokenizer(block, " ", true, false);

                int start_site = Text_utils::to_int( bt->next_token() );
                int end_site = Text_utils::to_int( bt->next_token() );
                double weight = Text_utils::to_double( bt->next_token() );

                if(start_site < 0 || end_site < start_site || end_site > site)
                {
                    cout<<"Error reading the graph input: edge coordinates at site "<<site<<" appear incorrect.\nExiting.\n\n";
                    exit(-1);
                }

                if(weight < 0 || weight > 1 || weight + sum_weight > 1)
                {
                    cout<<"Warning reading the graph input: edge weight at site "<<site<<" appear incorrect.\n\n";
                }

                if(weight + sum_weight > 1)
                {
                    cout<<"Warning reading the graph input: edge weight at site "<<site<<" appear incorrect.\n\n";
                    weight = 1.0-sum_weight;
                    sum_weight = 1.0;
                }

                Seq_edge edge;
                edge.start_site = start_site;
                edge.end_site = end_site;
                edge.weight = weight;

                edges.push_back(edge);

                delete bt;

                block = st->next_token();
                block = Text_utils::remove_surrounding_whitespaces(block);
            }
        }
    }

    // Addition of the last sequence in file
    if((name != "") && (sequence != ""))
    {
//        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );

        Fasta_entry fe;
        fe.name = name;
        fe.comment = comment;
        fe.sequence = sequence;
        fe.edges = edges;
        seqs.push_back(fe);
    }

}

/****************************************************************************************/

void Fasta_reader::write(ostream & output, const vector<Fasta_entry> & seqs) const throw (Exception)
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Fasta_reader::write. Failed to open file"); }

    string seq, temp = "";  // Initialization

    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        output << ">" << vi->name;
        if(vi->comment != "")
            output << " " << vi->comment;

        output << endl;

        // Sequence cutting to specified characters number per line
        seq = vi->sequence;
        while (seq != "")
        {
            if (seq.size() > chars_by_line)
            {
                temp = string(seq.begin(), seq.begin() + chars_by_line);
                output << temp  << endl;
                seq.erase(seq.begin(), seq.begin() + chars_by_line);
            }
            else
            {
                output << seq << endl;
                seq = "";
            }
        }
    }
}

/****************************************************************************************/

void Fasta_reader::write_graph(ostream & output, Node * root) const throw (Exception)
{

    Sequence *sequence = root->get_sequence();
    output<<"# root node\n";

    string alpha = sequence->get_full_alphabet();

    output<<"0 start;\n";
    for(int i=1; i<sequence->sites_length()-1; i++)
    {
        Site *site = sequence->get_site_at(i);
        output<<i<<" "<<alpha.at(site->get_state())<<"; ";

        Edge *edge = site->get_first_bwd_edge();
        output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";

        while(site->has_next_bwd_edge())
        {
            edge = site->get_next_bwd_edge();
            output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";
        }
        output<<endl;

    }
    output<<sequence->sites_length()-1<<" end; ";

    Site *site = sequence->get_site_at(sequence->sites_length()-1);

    Edge *edge = site->get_first_bwd_edge();
    output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";

    while(site->has_next_bwd_edge())
    {
        edge = site->get_next_bwd_edge();
        output<<edge->get_start_site_index()<<" "<<edge->get_end_site_index()<<" "<<edge->get_posterior_weight()<<";";
    }
    output<<endl;

}

/****************************************************************************************/

bool Fasta_reader::check_alphabet(string alphabet, string full_alphabet, vector<Fasta_entry> & seqs) throw (Exception)
{

    // Check that alphabet is correct but use (faster?) build-in alphabet.
    //
    if(alphabet == "ACGT") {

        bool allow_gaps = Settings_handle::st.is("cds-seqfile");

        dna_pi[0] = dna_pi[1] = dna_pi[2] = dna_pi[3] = 0.0;

        bool chars_ok = true;
        vector<Fasta_entry>::iterator vi = seqs.begin();

        // Main loop : for all sequences in vector container
        for (; vi != seqs.end(); vi++)
        {
            // Convert U -> T and all uppercase
            this->rna_to_DNA(&vi->sequence);

            string::iterator si = vi->sequence.begin();
            for (;si != vi->sequence.end();si++)
            {
                char c = *si;
                switch (c)
                {
                    case 'A':
                        dna_pi[0]++;
                        break;
                    case 'C':
                        dna_pi[1]++;
                        break;
                    case 'G':
                        dna_pi[2]++;
                        break;
                    case 'T':
                        dna_pi[3]++;
                        break;
                    case '-':
                        if(!allow_gaps)
                        {
                            vi->sequence.erase(si);
                            si--;
                            chars_ok = false;
                        }
                        break;
                    default:
                        // Remove characters not in full alphabet
                        if(full_alphabet.find(c) == string::npos) {
                            if(Settings::noise>2)
                                cout<<"deleting "<<*si<<"\n";
                            vi->sequence.erase(si);
                            si--;
                            if(c!=' ')
                            {
                                chars_ok = false;
                            }
                       }
                }
            }
        }

        float tot = dna_pi[0]+dna_pi[1]+dna_pi[2]+dna_pi[3];

        dna_pi[0] /= tot;
        dna_pi[1] /= tot;
        dna_pi[2] /= tot;
        dna_pi[3] /= tot;

        return chars_ok;
    }
    else if(alphabet == "HRKQNEDSTGPACVIMLFYW")
    {
        bool chars_ok = true;
        vector<Fasta_entry>::iterator vi = seqs.begin();

        // Main loop : for all sequences in vector container
        for (; vi != seqs.end(); vi++)
        {
            string::iterator si = vi->sequence.begin();
            for (;si != vi->sequence.end();si++)
            {
                char c = *si;
                // Remove characters not in full alphabet
                if(full_alphabet.find(c) == string::npos) {
                    vi->sequence.erase(si);
                    si--;
                    if(c!=' ')
                    {
                        chars_ok = false;
                    }
                }
            }
        }

        return chars_ok;
    }

    return false;
}

/****************************************************************************************/

int Fasta_reader::check_sequence_data_type(const vector<Fasta_entry> &seqs)
{

    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    int dna = 0;
    int protein = 0;
    string dna_alphabet = "ACGTUN";
    string protein_alphabet = "HRKQNEDSTGPACVIMLFYW";

    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        string::const_iterator si = vi->sequence.begin();
        for (;si != vi->sequence.end();si++)
        {
            char c = *si;
            if(dna_alphabet.find(c) != string::npos)
            {
                dna++;
            }
            if(protein_alphabet.find(c) != string::npos)
            {
                protein++;
            }
        }
    }

    if( ((float)dna )/ (float)protein > 0.9)
        return Model_factory::dna;
    else
        return Model_factory::protein;
}

/****************************************************************************************/

void Fasta_reader::rna_to_DNA(string *sequence) const
{
//    transform( sequence->begin(), sequence->end(), sequence->begin(), (int(*)(int))toupper );
    size_t si = 0;
    si=sequence->find('U',si);
    while(si != string::npos) {
        sequence->replace(si,1,1,'T');
        si=sequence->find('U',si);
    }
}

/****************************************************************************************/

bool Fasta_reader::check_sequence_names(const vector<Fasta_entry> *sequences,const vector<Node*> *leaf_nodes) const
{
    unsigned int names_match = 0;
    vector<Fasta_entry>::const_iterator si = sequences->begin();
    for (; si != sequences->end(); si++)
    {
        string s_name = si->name;
        vector<Node*>::const_iterator ni = leaf_nodes->begin();
        for (; ni != leaf_nodes->end(); ni++)
        {
            if((*ni)->get_name() == s_name)
                names_match++;
        }
    }

    // All the leafs in the guidetree need a sequence. Not all sequences need to be in the tree.

    int overlap = 0;

    if(names_match != leaf_nodes->size())
    {
        if(Settings_handle::st.is("cds-seqfile"))
        {
            cout<<"\nAll leaf node names in the guidetree file '"<<Settings_handle::st.get("cds-treefile").as<string>()<<"' do not match\n"<<
                " with a sequence name in the sequence file '"<<Settings_handle::st.get("cds-seqfile").as<string>()<<"'."<<endl;
        }
        else
        {
            cout<<"\nAll leaf node names in the guidetree file '"<<Settings_handle::st.get("treefile").as<string>()<<"' do not match\n"<<
                " with a sequence name in the sequence file '"<<Settings_handle::st.get("seqfile").as<string>()<<"'."<<endl;
        }

        set<string> snames;
        for (si = sequences->begin(); si != sequences->end(); si++)
        {
            snames.insert(si->name);
        }

        set<string> tnames;
        for (vector<Node*>::const_iterator ni = leaf_nodes->begin(); ni != leaf_nodes->end(); ni++)
        {
            tnames.insert((*ni)->get_name());
        }


        cout<<endl<<"Leaf names not in the sequences:"<<endl;
        for (vector<Node*>::const_iterator ni = leaf_nodes->begin(); ni != leaf_nodes->end(); ni++)
        {
            if(snames.find((*ni)->get_name()) == snames.end() )
            {
                cout<<" "<<(*ni)->get_name()<<endl;
                (*ni)->has_sequence(false);
            }
            else
            {
                (*ni)->has_sequence(true);
                overlap++;
            }
        }
    }
    else
    {
        overlap = names_match;
    }

    // Not all sequences need to be in the tree but give a warning that names don't match.
    if(names_match < 2)
    {
        cout<<"\nNo sequences to align! Exiting.\n\n";
        exit(0);
    }

    if((int)sequences->size() > overlap && overlap == (int)leaf_nodes->size())
    {
        cout<<"\nWarning: "<<leaf_nodes->size()<<" leaf nodes but "<<sequences->size()<<" sequences! Excess sequences will be removed.\n\n";
        return true;
    }
    if(overlap < (int)leaf_nodes->size())
    {
        cout<<"\nWarning: "<<leaf_nodes->size()<<" leaf nodes but "<<overlap<<" matching sequences! Excess branches will be removed.\n\n";
        return false;
    }
    return true;
}

/****************************************************************************************/

void Fasta_reader::place_sequences_to_nodes(const vector<Fasta_entry> *sequences,vector<Node*> *leaf_nodes, string full_char_alphabet, bool gapped)
{

    vector<Fasta_entry>::const_iterator si = sequences->begin();
    for (; si != sequences->end(); si++)
    {
        string s_name = si->name;

        vector<Node*>::iterator ni = leaf_nodes->begin();
        for (; ni != leaf_nodes->end(); ni++)
        {
            if((*ni)->get_name() == s_name) {
                (*ni)->add_name_comment( si->comment );
                (*ni)->add_sequence( *si, full_char_alphabet, gapped);
            }
        }
    }
}
