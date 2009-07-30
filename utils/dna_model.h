#ifndef DNA_MODEL_H
#define DNA_MODEL_H

#include "db_matrix.h"
#include "int_matrix.h"
#include <string>

namespace ppa {

class Dna_model
{
    std::string full_dna_alphabet;
public:
    Dna_model(std::string alpha,float dist);
    ~Dna_model();

    Db_matrix *dnaPi;
    Db_matrix *dnaPr;
    Db_matrix *logDnaPi;
    Db_matrix *logDnaPr;
    Int_matrix *parsimony_table;

    float distance;

    float id_prob;
    float ext_prob;
    float match_prob;

    float log_id_prob;
    float log_ext_prob;
    float log_match_prob;

    float ins_rate;
    float del_rate;
    float ins_prob;
    float del_prob;

    float gap_open() { return id_prob; }
    float gap_close() { return id_prob; }
    float gap_ext() { return ext_prob; }
    float non_gap() { return match_prob; }

    float log_gap_open() { return log_id_prob; }
    float log_gap_close() { return log_id_prob; }
    float log_gap_ext() { return log_ext_prob; }
    float log_non_gap() { return log_match_prob; }


    float log_score(int i,int j) { return logDnaPr->g(i,j); }
    float score(int i,int j) { return dnaPr->g(i,j); }

    int parsimony_state(int i,int j) { return parsimony_table->g(i,j); }

    std::string get_full_alphabet() { return full_dna_alphabet; }

};
}
#endif // DNA_MODEL_H
