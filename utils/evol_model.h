#ifndef EVOL_MODEL_H
#define EVOL_MODEL_H

#include "db_matrix.h"
#include "int_matrix.h"
#include <string>

namespace ppa {

class Evol_model
{
    std::string full_char_alphabet;
public:
    Evol_model(std::string alpha,float dist);
    ~Evol_model();

    Db_matrix *charPi;
    Db_matrix *charPr;
    Db_matrix *logCharPi;
    Db_matrix *logCharPr;
    Int_matrix *parsimony_table;

    float distance;

    float id_prob;
    float ext_prob;
    float end_ext_prob;
    float break_ext_prob;
    float match_prob;

    float log_id_prob;
    float log_ext_prob;
    float log_end_ext_prob;
    float log_break_ext_prob;
    float log_match_prob;

    float ins_rate;
    float del_rate;
    float ins_prob;
    float del_prob;

    float gap_open() { return id_prob; }
    float gap_close() { return id_prob; }
    float gap_ext() { return ext_prob; }
    float gap_end_ext() { return end_ext_prob; }
    float gap_break_ext() { return break_ext_prob; }
    float non_gap() { return match_prob; }

    float log_gap_open() { return log_id_prob; }
    float log_gap_close() { return log_id_prob; }
    float log_gap_ext() { return log_ext_prob; }
    float log_gap_end_ext() { return log_end_ext_prob; }
    float log_gap_break_ext() { return log_break_ext_prob; }
    float log_non_gap() { return log_match_prob; }


    float log_score(int i,int j) { return logCharPr->g(i,j); }
    float score(int i,int j) { return charPr->g(i,j); }

    int parsimony_state(int i,int j) { return parsimony_table->g(i,j); }

    std::string get_full_alphabet() { return full_char_alphabet; }

};
}
#endif // EVOL_MODEL_H
