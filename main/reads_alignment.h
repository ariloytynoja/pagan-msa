#ifndef READS_ALIGNMENT_H
#define READS_ALIGNMENT_H

#include <vector>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/model_factory.h"
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "utils/node.h"
#include "main/sequence.h"
#include "main/simple_alignment.h"

using namespace std;

namespace ppa
{

class Reads_alignment
{
    Node *global_root;
public:
    Reads_alignment();
    void align(Node *root, Model_factory *mf,int count);

    Node *get_global_root() { return global_root; }
};
}

#endif // READS_ALIGNMENT_H
