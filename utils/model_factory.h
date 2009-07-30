#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include <string>
#include <vector>
#include "db_matrix.h"
#include "int_matrix.h"
#include "dna_model.h"
#include "settings.h"

namespace ppa{

struct Dna_symbol
{
    int index;
    char symbol;
    int n_bases;
    std::vector<char> bases;
};

class Model_factory
{

    std::string dna_alphabet;
    std::string full_dna_alphabet;
    int dna_as;  // alphabet size
    int dna_fas; // full_alphabet size

    Db_matrix *dnaPi;
    float dna_ins_rate;
    float dna_del_rate;
    float dna_ext_prob;

    Int_matrix *parsimony_table;

    Db_matrix * dnaU;
    Db_matrix * dnaV;
    Db_matrix * dnaRoot;

    int i;
    int j;
    int k;
    int l;

    std::vector<Dna_symbol> dna_letters;

    void define_dna_alphabet();
    void print_dna_alphabet();

    void build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot);
    void print_dna_q_matrices(Db_matrix *dnaQ);
//    void print_dna_p_matrices(Dna_model &model);

public:
    Model_factory();
    ~Model_factory();

    void dna_model(float *dna_pi,Settings *st);
    void dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob);

    Dna_model dna_alignment_model(double distance);

    void print_dna_p_matrices(Dna_model &model);

    std::string get_dna_alphabet() { return dna_alphabet; }
    std::string get_full_dna_alphabet() { return full_dna_alphabet; }

};

}
#endif // MODEL_FACTORY_H
