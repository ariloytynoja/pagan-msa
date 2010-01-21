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
#include "fasta_reader.h"
#include "text_utils.h"

using namespace ppa;

/****************************************************************************************/

void Fasta_reader::read(istream & input, vector<Fasta_entry> & seqs, bool short_names = false) const throw (Exception)
{
    if (!input) { throw IOException ("Fasta_reader::read. Failed to open file"); }

    string temp, name, comment, sequence = "";  // Initialization

    // Main loop : for all file lines
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
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );

        Fasta_entry fe;
        fe.name = name;
        fe.comment = comment;
        fe.sequence = sequence;
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
    string dna_alphabet = "ACGTU";
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

void Fasta_reader::check_sequence_names(const vector<Fasta_entry> *sequences,const vector<Node*> *leaf_nodes,const Settings *st) const
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

    if(names_match != leaf_nodes->size())
    {
        cout<<"All leaf node names in the guidetree file '"<<st->get("treefile").as<string>()<<"'\ndo not match"<<
            " with a sequence name in the sequence file '"<<st->get("seqfile").as<string>()<<"'."<<endl;

        cout<<endl<<"Sequence names:"<<endl;
        for (si = sequences->begin(); si != sequences->end(); si++)
            cout<<" "<<si->name<<endl;

        cout<<endl<<"Leaf names:"<<endl;
        for (vector<Node*>::const_iterator ni = leaf_nodes->begin(); ni != leaf_nodes->end(); ni++)
            cout<<" "<<(*ni)->get_name()<<endl;

        cout<<endl<<"Exiting."<<endl;
        exit(-1);
    }

    // Not all sequences need to be in the tree but give a warning that names don't match.

    if(sequences->size() != leaf_nodes->size())
    {
        cout<<"\nWarning: "<<leaf_nodes->size()<<" leaf nodes but "<<sequences->size()<<" sequences!\n";
    }
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
                (*ni)->add_sequence( si->sequence, full_char_alphabet, gapped);
            }
        }
    }
}
