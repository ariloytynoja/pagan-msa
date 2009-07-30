#include "dna_model.h"

using namespace ppa;

Dna_model::Dna_model(std::string alpha,float dist)
{
    full_dna_alphabet = alpha;
    int dna_fas = full_dna_alphabet.length();

    distance = dist;

    dnaPi = new Db_matrix(dna_fas,"pi_dna");
    dnaPr = new Db_matrix(dna_fas,dna_fas,"P_dna");
    dnaPr->initialise(0);

    logDnaPi = new Db_matrix(dna_fas,"logpi_dna");
    logDnaPr = new Db_matrix(dna_fas,dna_fas,"logP_dna");
    logDnaPr->initialise(-HUGE_VAL);

    parsimony_table = new Int_matrix(dna_fas,dna_fas,"parsimony_dna");

}

Dna_model::~Dna_model()
{
    delete dnaPi;
    delete dnaPr;

    delete logDnaPi;
    delete logDnaPr;

    delete parsimony_table;
}
