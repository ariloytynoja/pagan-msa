/***************************************************************************
 *   Copyright (C) 2010-2012 by Ari Loytynoja                              *
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

#ifndef FASTA_ENTRY_H
#define FASTA_ENTRY_H

#include <string>
#include <vector>

using namespace std;

namespace ppa
{

struct Seq_edge
{
    int start_site;
    int end_site;
    float weight;
};

struct Fasta_entry
{
    string name;
    string comment;
    string sequence;
    string dna_sequence;
    string quality;
    vector<Seq_edge> edges;
    int data_type;
    string tid;
    float node_score;
    string node_to_align;
    int first_read_length; // for pair-end reads
    int node_start_pos1;   // to check overlapping reads
    int node_end_pos1;
    int node_start_pos2;   // for pair-end reads
    int node_end_pos2;
    int trim_start;
    int trim_end;
    int local_qstart;      // local alignment query start
    int local_qend;        // local alignment query end
    int local_tstart;      // local alignment target start
    int local_tend;        // local alignment target end
    bool use_local;
    int cluster_attempts;
    bool reversed;
    int num_duplicates;
};

}

#endif // FASTA_ENTRY_H
