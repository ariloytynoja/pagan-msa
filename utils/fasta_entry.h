#ifndef FASTA_ENTRY_H
#define FASTA_ENTRY_H

#include <string>

using namespace std;

namespace ppa
{

struct Fasta_entry
{
    string name;
    string comment;
    string sequence;
    string quality;
};
}

#endif // FASTA_ENTRY_H
