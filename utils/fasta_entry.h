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
};

}

#endif // FASTA_ENTRY_H
