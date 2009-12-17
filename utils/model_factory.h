#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include <string>
#include <vector>
#include "db_matrix.h"
#include "int_matrix.h"
#include "evol_model.h"
#include "settings.h"

namespace ppa{

struct Char_symbol
{
    int index;
    char symbol;
    int n_residues;
    std::vector<char> residues;
};

class Model_factory
{

    std::string char_alphabet;
    std::string full_char_alphabet;
    int char_as;  // alphabet size
    int char_fas; // full_alphabet size
    int sequence_data_type;

    Db_matrix *charPi;
    float char_ins_rate;
    float char_del_rate;
    float char_ext_prob;

    Int_matrix *parsimony_table;

    Db_matrix * charU;
    Db_matrix * charV;
    Db_matrix * charRoot;

    int i;
    int j;
    int k;
    int l;

    std::vector<Char_symbol> char_letters;

    void define_dna_alphabet();
    void define_protein_alphabet();
    void print_char_alphabet();

    void build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot);
    void print_char_q_matrices(Db_matrix *charQ);
//    void print_char_p_matrices(Evol_model &model);

public:
    Model_factory(int sequence_data_type);
    ~Model_factory();

    enum Data_type {dna,protein,codon};

    void dna_model(float *char_pi,Settings *st);
    void dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob);

    void protein_model(Settings *st);
    void protein_model(float ins_rate,float del_rate, float ext_prob);

    Evol_model alignment_model(double distance);
//    Evol_model char_alignment_model(double distance);

    void print_char_p_matrices(Evol_model &model);

    std::string get_char_alphabet() { return char_alphabet; }
    std::string get_full_char_alphabet() { return full_char_alphabet; }

};

}
#endif // MODEL_FACTORY_H
