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
    void find_nodes_for_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf);
    void remove_overlapping_reads(vector<Fasta_entry> *reads, Model_factory *mf);
    void remove_target_overlapping_identical_reads(vector<Fasta_entry> *reads_for_this, Model_factory *mf);
    void remove_target_overlapping_reads(vector<Fasta_entry> *reads_for_this);
    void align_two_reads(Node *node, Fasta_entry *ri1, Fasta_entry *ri2, Model_factory *mf);
    int reads_pairwise_matching_sites(Node *node);

    double read_match_score(Node *node, Fasta_entry *read, Model_factory *mf, float best_score);
    bool read_alignment_overlaps(Node * node, string read_name, string ref_node_name);
    void add_trimming_comment(vector<Fasta_entry> *reads);
    void merge_paired_reads(vector<Fasta_entry> *reads, Model_factory *mf);
    void find_paired_reads(vector<Fasta_entry> *reads);
    void copy_node_details(Node *reads_node,Fasta_entry *read, string full_alpha);
    bool correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map);

    static bool better_score(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.node_score>b.node_score);
    }

    void sort_reads_vector(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_alignment::better_score);
    }
    
    static bool nodeIsSmaller(const string& l,const string& r)
    {   char a,b;
        int vl,vr;
        stringstream ls(l);
        stringstream rs(r);
        ls>>a>>vl>>b;
        rs>>a>>vr>>b;

        return (vl<vr);
    }
public:
    Reads_alignment();
    void align(Node *root, Model_factory *mf,int count);

    void merge_reads_only();

    Node *get_global_root() { return global_root; }
};
}

#endif // READS_ALIGNMENT_H
