#include "evol_model.h"

using namespace ppa;

Evol_model::Evol_model(std::string alpha,float dist)
{
    full_char_alphabet = alpha;
    int char_fas = full_char_alphabet.length();

    distance = dist;

    charPi = new Db_matrix(char_fas,"pi_char");
    charPr = new Db_matrix(char_fas,char_fas,"P_char");
    charPr->initialise(0);

    logCharPi = new Db_matrix(char_fas,"logpi_char");
    logCharPr = new Db_matrix(char_fas,char_fas,"logP_char");
    logCharPr->initialise(-HUGE_VAL);

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

}

Evol_model::~Evol_model()
{
    delete charPi;
    delete charPr;

    delete logCharPi;
    delete logCharPr;

    delete parsimony_table;
}
