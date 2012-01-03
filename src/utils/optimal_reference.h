#ifndef OPTIMAL_REFERENCE_H
#define OPTIMAL_REFERENCE_H

#include "utils/fasta_entry.h"
#include "main/node.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace ppa{

struct hit {
    std::string query;
    std::string target;
    int score;
    int q_start;
    int q_end;
    char q_strand;
    int t_start;
    int t_end;
    char t_strand;
};

class Optimal_reference
{
    Node *global_root;
    bool find_best_exonerate_hit(Fasta_entry *q, vector<Fasta_entry> *sequences, Fasta_entry *best_read);
    void find_best_hit(Fasta_entry *q, vector<Fasta_entry> *sequences, Model_factory *mf, Fasta_entry *best_read);
    void copy_node_details(Node *reads_node,Fasta_entry *read,bool turn_revcomp = false);
    float read_alignment_overlap(Node * node, string read_name, string ref_node_name);
    float read_alignment_identity(Node * node, string read_name, string ref_node_name);
    void read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity);
    bool split_sugar_string(const std::string& row,hit *h);
    void write_exonerate_input(Fasta_entry *q, vector<Fasta_entry> *seqs,  int r);
    void write_exonerate_input(string q, vector<Fasta_entry> *seqs, int r);
    void delete_files(int r);

public:
    Optimal_reference();
    void find_optimal_reference(vector<Fasta_entry> *sequences);
    void align(Node *root, Model_factory *mf, int count);
    Node *get_global_root() { return global_root; }
};
}

#endif // OPTIMAL_REFERENCE_H
