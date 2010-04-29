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
    string quality;
    vector<Seq_edge> edges;
    string tid;
    float node_score;
    string node_to_align;
    int first_read_length; // for pair-end reads
    int node_start_pos1;   // to check overlapping reads
    int node_end_pos1;
    int node_start_pos2;   // for pair-end reads
    int node_end_pos2;
};

}

#endif // FASTA_ENTRY_H
