/***************************************************************************
 *   Copyright (C) 2010 by Ari Loytynoja                                   *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include "utils/exceptions.h"
#include "main/node.h"
#include "utils/fasta_entry.h"

using namespace std;

namespace ppa
{

class Fasta_reader
{
    unsigned int chars_by_line;
    float dna_pi[4];

    std::map<std::string,std::string> codon_to_aa;
    std::map<std::string,std::string> aa_to_codon;

    void rna_to_DNA(string *sequence) const;
    void define_translation_tables();
    string DNA_to_protein(string *sequence) const;
    string protein_to_DNA(string *dna,string *prot) const;

public:
    Fasta_reader() : chars_by_line(60) {
        this->define_translation_tables();
    }

    void set_chars_by_line(int n) { chars_by_line = n; }

    void read(istream & input, vector<Fasta_entry> & seqs, bool short_names) const throw (Exception);
    void read(const string & path, vector<Fasta_entry> & seqs, bool short_names) const throw (Exception)
    {
        ifstream input(path.c_str(), ios::in);
        read(input, seqs, short_names);
        input.close();
    }
    void read_fasta(istream & input, vector<Fasta_entry> & seqs, bool short_names) const throw (Exception);
    void read_fastq(istream & input, vector<Fasta_entry> & seqs) const throw (Exception);
    void read_graph(istream & input, vector<Fasta_entry> & seqs, bool short_names) const throw (Exception);

    void trim_fastq_reads(vector<Fasta_entry> * seqs) const throw (Exception);

    void write(ostream & output, const vector<Fasta_entry> & seqs) const throw (Exception);
    void write(const string & path, const vector<Fasta_entry> & seqs, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".fas").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        write(output, seqs);
        output.close();
    }

    void write_dna(ostream & output, const vector<Fasta_entry> & seqs, const vector<Fasta_entry> & org_seqs, Node *root) const throw (Exception);
    void write_dna(const string & path, const vector<Fasta_entry> & seqs, const vector<Fasta_entry> & org_seqs, Node *root, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".dna.fas").c_str(), overwrite ? (ios::out) : (ios::out|ios::app) );
        write_dna(output, seqs, org_seqs, root);
        output.close();
    }

    void print_fast_entry(ostream & output, const Fasta_entry *entry) const;

    void write_fastq(ostream & output, const vector<Fasta_entry> & seqs) const throw (Exception);
    void write_fastq(const string & path, const vector<Fasta_entry> & seqs, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".fastq").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        write_fastq(output, seqs);
        output.close();
    }

    void write_anctree(const string outfile, Node *root, bool overwrite=true)
    {
        ofstream output( (outfile+".anctree").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        if (! output) { throw IOException ("Fasta_reader::write_anctree. Failed to open file"); }

        output<<root->print_tree(true);
        output.close();
    }

    void write_graph(ostream & output, Node * root) const throw (Exception);
    void write_graph(const string & path, Node * root, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".grp").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        write_graph(output, root);
        output.close();
    }


    bool check_alphabet(vector<Fasta_entry> *sequences, int data_type = -1)  throw (Exception);
    bool check_sequence_names(const vector<Fasta_entry> *sequences,const vector<Node*> *leaf_nodes) const;

    float* base_frequencies() { return dna_pi; }
    int check_sequence_data_type(const vector<Fasta_entry> * sequences) const;

    void place_sequences_to_nodes(const vector<Fasta_entry> *sequences,vector<Node*> *leaf_nodes, bool gapped = false, int data_type = -1);
};
}

#endif // FASTA_READER_H

