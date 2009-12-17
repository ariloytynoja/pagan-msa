
#include <iomanip>
#include <iostream>
#include <vector>
#include "model_factory.h"
#include "settings.h"
#include "settings_handle.h"
#include "eigen.h"

using namespace ppa;
using namespace std;

Model_factory::Model_factory(int s)
{

    sequence_data_type = s;

    if(sequence_data_type == Model_factory::dna)
        this->define_dna_alphabet();
    else if(sequence_data_type == Model_factory::protein)
        this->define_protein_alphabet();
    else if(sequence_data_type == Model_factory::codon)
        ;
    else
    {
        cout<<"Model_factory(): invalid sequence data type. Exiting\n";
        exit(-1);
    }
}

Model_factory::~Model_factory()
{
    if(charPi!=0)
        delete charPi;

    if(charU!=0)
        delete charU;

    if(charV!=0)
        delete charV;

    if(charRoot!=0)
        delete charRoot;

    if(parsimony_table!=0)
        delete parsimony_table;
}

/*******************************************/

/*
 * Definition of DNA alpahbet including ambiguity characters.
 */
void Model_factory::define_dna_alphabet()
{
    char_alphabet = "ACGT";
    full_char_alphabet = "ACGTRYMKWSBDHVN";

    char_as = char_alphabet.length();
    char_fas = full_char_alphabet.length();

    int n_residues[] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};
    string ambiguity[] = {"A","C","G","T","AG","CT","AC","GT","AT","CG","CGT","AGT","ACT","ACG","ACGT"};

    for(int i=0;i<15;i++)
    {
        Char_symbol letter;

        letter.index = i;
        letter.symbol = full_char_alphabet.at(i);
        letter.n_residues = n_residues[i];
        for(int j=0;j<letter.n_residues;j++)
            letter.residues.push_back(ambiguity[i].at(j));
        char_letters.push_back(letter);
    }

    /*
     * Table for resolving the parental state from the parsimony ambiguity alphabet.
     */

    int a = 1;
    int c = 2;
    int g = 4;
    int t = 8;

    int r = (a|g);
    int y = (c|t);
    int m = (a|c);
    int k = (g|t);
    int w = (a|t);
    int s = (c|g);

    int b = (c|g|t);
    int d = (a|g|t);
    int h = (a|c|t);
    int v = (a|c|g);

    int n = (a|c|g|t);

    Int_matrix *bin2pos = new Int_matrix(n+1,"bin2pos");

    bin2pos->s(0,a);
    bin2pos->s(1,c);
    bin2pos->s(2,g);
    bin2pos->s(3,t);

    bin2pos->s(4,r);
    bin2pos->s(5,y);
    bin2pos->s(6,m);
    bin2pos->s(7,k);
    bin2pos->s(8,w);
    bin2pos->s(9,s);

    bin2pos->s(10,b);
    bin2pos->s(11,d);
    bin2pos->s(12,h);
    bin2pos->s(13,v);

    bin2pos->s(14,n);

    Int_matrix *pos2bin = new Int_matrix(n+1,"pos2bin");

    pos2bin->s(a,0);
    pos2bin->s(c,1);
    pos2bin->s(g,2);
    pos2bin->s(t,3);

    pos2bin->s(r,4);
    pos2bin->s(y,5);
    pos2bin->s(m,6);
    pos2bin->s(k,7);
    pos2bin->s(w,8);
    pos2bin->s(s,9);

    pos2bin->s(b,10);
    pos2bin->s(d,11);
    pos2bin->s(h,12);
    pos2bin->s(v,13);

    pos2bin->s(n,14);

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            int v = (pos2bin->g(i)&pos2bin->g(j));
            if(v>0)
                parsimony_table->s(bin2pos->g(v),i,j);
            else
                parsimony_table->s(bin2pos->g((pos2bin->g(i)|pos2bin->g(j))),i,j);
        }
    }

    delete bin2pos;
    delete pos2bin;

    if(Settings::noise>5)
    {
        cout<<"\nModel_factory::define_dna_alphabet(). DNA parsimony table.\n\n  ";
        for(int i=0;i<15;i++)
            cout<<full_char_alphabet.at(i)<<" ";
        cout<<endl;

        for(int i=0;i<15;i++)
        {
            cout<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<15;j++)
            {
                cout<<full_char_alphabet.at(parsimony_table->g(i,j))<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
}
/*******************************************/

/*
 * Definition of protein alpahbet including ambiguity characters.
 */
void Model_factory::define_protein_alphabet()
{
    char_alphabet = "HRKQNEDSTGPACVIMLFYW";
    full_char_alphabet = "RKQEDNGHASTPCIVMLFYWZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZX";
    full_char_alphabet = "RKQEDNGHASTPCIVMLFYWabcdefghijklmnopqrstuvxyz12345X";

    char_as = char_alphabet.length();
    char_fas = full_char_alphabet.length();

    int n_residues[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                        2,2,2,2,2,2,2,3,2,4,3,3,5,3,3,3,4,4,7,2,3,2,10,11,12,3,4,16,3,4,20};
    string ambiguity[] = {"R","K","Q","E","D","N","G","H","A","S","T","P","C","I","V","M","L","F","Y","W",
                          "NG","HA","IV","ST","QE","ML","RK","RKQ","ED","NAST","AST","HML","RKHSA","NED",
                          "TIV","MLF","HRKQ","ASTG","HRKQSTA","CV","IML","FY","HRKQNEDSTA","HRKQNEDSTPA",
                          "HRKQNEDSTGPA","CIV","MLFY","HRKQNEDSTGPACVIM","LFY","LFYW","HRKQNEDSTGPACVIMLFYW"};

    for(int i=0;i<51;i++)
    {
        Char_symbol letter;

        letter.index = i;
        letter.symbol = full_char_alphabet.at(i);
        letter.n_residues = n_residues[i];
        for(int j=0;j<letter.n_residues;j++)
            letter.residues.push_back(ambiguity[i].at(j));
        char_letters.push_back(letter);
    }

    /*
     * Table for resolving the parental state from the parsimony ambiguity alphabet.
     */
    int table[] = {
         0, 26, 27, 42, 42, 42, 44, 32, 32, 32, 38, 43, 47, 47, 47, 47, 50, 50, 50, 50,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        26,  1, 27, 42, 42, 42, 44, 32, 32, 32, 38, 43, 47, 47, 47, 47, 50, 50, 50, 50,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
        27, 27,  2, 24, 42, 42, 44, 36, 38, 38, 38, 43, 47, 47, 47, 47, 50, 50, 50, 50,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        42, 42, 24,  3, 28, 33, 44, 42, 42, 42, 42, 43, 47, 47, 47, 47, 50, 50, 50, 50,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
        42, 42, 42, 28,  4, 33, 44, 42, 42, 42, 42, 43, 47, 47, 47, 47, 50, 50, 50, 50,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        42, 42, 42, 33, 33,  5, 20, 42, 29, 29, 29, 43, 47, 47, 47, 47, 50, 50, 50, 50,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
        44, 44, 44, 44, 44, 20,  6, 44, 37, 37, 37, 44, 47, 47, 47, 47, 50, 50, 50, 50,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
        32, 32, 36, 42, 42, 42, 44,  7, 21, 32, 38, 43, 47, 47, 47, 31, 31, 50, 50, 50,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
        32, 32, 38, 42, 42, 29, 37, 21,  8, 29, 29, 43, 47, 47, 47, 47, 50, 50, 50, 50,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
        32, 32, 38, 42, 42, 29, 37, 32, 29,  9, 23, 43, 47, 47, 47, 47, 50, 50, 50, 50,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
        38, 38, 38, 42, 42, 29, 37, 38, 29, 23, 10, 43, 47, 34, 34, 47, 50, 50, 50, 50, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        43, 43, 43, 43, 43, 43, 44, 43, 43, 43, 43, 11, 47, 47, 47, 47, 50, 50, 50, 50, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
        47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 12, 45, 39, 47, 50, 50, 50, 50, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 34, 47, 45, 13, 22, 40, 40, 50, 50, 50, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
        47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 34, 47, 39, 22, 14, 47, 50, 50, 50, 50, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        47, 47, 47, 47, 47, 47, 47, 31, 47, 47, 47, 47, 47, 40, 47, 15, 25, 35, 46, 50, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        50, 50, 50, 50, 50, 50, 50, 31, 50, 50, 50, 50, 50, 40, 50, 25, 16, 35, 46, 49, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 35, 35, 17, 41, 49, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 46, 46, 41, 18, 49, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 49, 49, 49, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 44, 47, 44, 44, 50, 44, 44, 44, 44, 44, 50, 44, 44, 47, 50, 44, 44, 44, 47, 50, 50, 44, 44, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 21, 47, 38, 42, 50, 32, 38, 42, 42, 38, 50, 32, 42, 47, 50, 38, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 47, 47, 22, 47, 47, 50, 47, 47, 47, 47, 47, 50, 47, 47, 34, 50, 47, 47, 47, 45, 50, 50, 47, 47, 47, 45, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 38, 47, 23, 42, 50, 38, 38, 42, 29, 29, 50, 38, 42, 47, 50, 38, 37, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 42, 47, 42, 24, 50, 42, 42, 42, 42, 42, 50, 42, 42, 47, 50, 42, 44, 42, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 25, 50, 50, 50, 50, 50, 31, 50, 50, 50, 35, 50, 50, 50, 50, 40, 46, 50, 50, 50, 50, 46, 50, 46, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 32, 47, 38, 42, 50, 26, 27, 42, 42, 38, 50, 32, 42, 47, 50, 36, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 38, 47, 38, 42, 50, 27, 27, 42, 42, 38, 50, 38, 42, 47, 50, 36, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 42, 47, 42, 42, 50, 42, 42, 28, 42, 42, 50, 42, 33, 47, 50, 42, 44, 42, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 42, 47, 29, 42, 50, 42, 42, 42, 29, 29, 50, 42, 42, 47, 50, 42, 44, 42, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 38, 47, 29, 42, 50, 38, 38, 42, 29, 29, 50, 38, 42, 47, 50, 38, 37, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 31, 50, 50, 50, 50, 50, 31, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 32, 47, 38, 42, 50, 32, 38, 42, 42, 38, 50, 32, 42, 47, 50, 38, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 42, 47, 42, 42, 50, 42, 42, 33, 42, 42, 50, 42, 33, 47, 50, 42, 44, 42, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 47, 47, 34, 47, 47, 50, 47, 47, 47, 47, 47, 50, 47, 47, 34, 50, 47, 47, 47, 47, 50, 50, 47, 47, 47, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 35, 50, 50, 50, 50, 50, 50, 50, 50, 50, 35, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 46, 50, 46, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 38, 47, 38, 42, 50, 36, 36, 42, 42, 38, 50, 38, 42, 47, 50, 36, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 44, 47, 37, 44, 50, 44, 44, 44, 44, 37, 50, 44, 44, 47, 50, 44, 37, 44, 47, 50, 50, 44, 44, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 38, 47, 38, 42, 50, 38, 38, 42, 42, 38, 50, 38, 42, 47, 50, 38, 44, 38, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 47, 47, 45, 47, 47, 50, 47, 47, 47, 47, 47, 50, 47, 47, 47, 50, 47, 47, 47, 39, 50, 50, 47, 47, 47, 45, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 40, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 40, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 41, 50, 50, 50, 50, 46, 50, 46, 49, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 42, 47, 42, 42, 50, 42, 42, 42, 42, 42, 50, 42, 42, 47, 50, 42, 44, 42, 47, 50, 50, 42, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 43, 47, 43, 43, 50, 43, 43, 43, 43, 43, 50, 43, 43, 47, 50, 43, 44, 43, 47, 50, 50, 43, 43, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 44, 44, 47, 44, 44, 50, 44, 44, 44, 44, 44, 50, 44, 44, 47, 50, 44, 44, 44, 47, 50, 50, 44, 44, 44, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 47, 47, 45, 47, 47, 50, 47, 47, 47, 47, 47, 50, 47, 47, 47, 50, 47, 47, 47, 45, 50, 50, 47, 47, 47, 45, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 46, 50, 46, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 47, 47, 47, 47, 47, 50, 47, 47, 47, 47, 47, 50, 47, 47, 47, 50, 47, 47, 47, 47, 50, 50, 47, 47, 47, 47, 50, 47, 50, 50, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 50, 46, 50, 50, 50, 50, 46, 50, 46, 49, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 49, 50, 50, 50, 50, 50, 50, 49, 49, 50,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
     };

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            parsimony_table->s(table[i*char_fas+j],i,j);
        }
    }


    if(Settings::noise>5)
    {
        cout<<"\nModel_factory::define_protein_alphabet(). Protein parsimony table.\n\n  ";
        for(int i=0;i<char_fas;i++)
            cout<<full_char_alphabet.at(i)<<" ";
        cout<<endl;

        for(int i=0;i<char_fas;i++)
        {
            cout<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<char_fas;j++)
            {
                cout<<full_char_alphabet.at(parsimony_table->g(i,j))<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
}
/*******************************************/

/*
 * For debugging.
 */
void Model_factory::print_char_alphabet()
{
    cout<<"\nModel_factory::print_char_alphabet()\n";
    for(unsigned int i=0;i<char_letters.size();i++)
    {
        Char_symbol a = char_letters.at(i);

        cout<<"Index "<<a.index<<"; symbol "<<a.symbol<<"; residues ("<<a.n_residues<<"): ";
        for(int j=0;j<a.n_residues;j++){
            int t = char_alphabet.find(a.residues.at(j));
            cout<<a.residues.at(j)<<"["<<t<<"] ";
        }
        cout<<endl;
    }
}

/*******************************************/


/*
 * Definition of DNA model: Q matrix and indel rate/prob.
 */

void Model_factory::dna_model(float *char_pi,Settings *st)
{
    float char_kappa = 2.0;
    if(st->is("char-kappa"))
        char_kappa =  st->get("char-kappa").as<float>();

    float char_rho = 1.0;
    if(st->is("char-rho"))
        char_rho =  st->get("char-rho").as<float>();

    float ins_rate = 0.05;
    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    float del_rate = 0.05;
    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    float gap_ext = 0.5;
    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    this->dna_model(char_pi,char_kappa,char_rho,ins_rate,del_rate,gap_ext);
}

void Model_factory::dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob)
{

    if (Settings::noise>4){
        cout<<"DNA substitution model: base frequencies "
                <<pi[0]<<", "<<pi[1]<<", "<<pi[2]<<" and "<<pi[3]
                <<"; kappa "<<kappa<<" and rho "<<rho<<"; ins rate "<<ins_rate<<", del rate "<<del_rate<<" and ext. prob "<<ext_prob<<"."<<endl;
    }

    char_ins_rate = ins_rate;
    char_del_rate = del_rate;
    char_ext_prob = ext_prob;


    if(Settings::noise > 4)
        print_char_alphabet();

    // pi and Q
    // Allocate space for P matrices
    //
    Db_matrix *charQ  = new Db_matrix(char_as,char_as,"Q_char");

    charPi = new Db_matrix(char_fas,"pi_char");

    FOR(j,char_as) {
        charPi->s(pi[j],j);
    }

    float ka = kappa/2.0;

    float piR = pi[0]+pi[2];
    float piY = pi[1]+pi[3];

    float beta = 1/(2*piR*piY*(1+ka));

    float alfaY = ( piR*piY*ka - pi[0]*pi[2] - pi[1]*pi[3]) /
                  ( (2+2*ka) * ( piY*pi[0]*pi[2] * rho+piR*pi[1]*pi[3] ) );

    float alfaR = rho*alfaY;

    /////////////////////////////
    /* filling up the Q matrix */
    double t;

    /*1st row*/
    t = beta*pi[1];
    charQ->s(t,0,1);

    t = alfaR*pi[2] / piR+beta*pi[2];
    charQ->s(t,0,2);

    t = beta*pi[3];
    charQ->s(t,0,3);

    t = 0-charQ->g(0,1)-charQ->g(0,2)-charQ->g(0,3);
    charQ->s(t,0,0);

    /*2nd row*/
    t = beta*pi[0];
    charQ->s(t,1,0);

    t = beta*pi[2];
    charQ->s(t,1,2);

    t = alfaY*pi[3]/piY+beta*pi[3];
    charQ->s(t,1,3);

    t = 0-charQ->g(1,0)-charQ->g(1,2)-charQ->g(1,3);
    charQ->s(t,1,1);

    /*3rd row*/
    t = alfaR*pi[0]/piR+beta*pi[0];
    charQ->s(t,2,0);

    t = beta*pi[1];
    charQ->s(t,2,1);

    t = beta*pi[3];
    charQ->s(t,2,3);

    t = 0-charQ->g(2,0)-charQ->g(2,1)-charQ->g(2,3);
    charQ->s(t,2,2);

    /*4th row*/
    t = beta*pi[0];
    charQ->s(t,3,0);

    t = alfaY*pi[1]/piY+beta*pi[1];
    charQ->s(t,3,1);

    t = beta*pi[2];
    charQ->s(t,3,2);

    t = 0-charQ->g(3,0)-charQ->g(3,1)-charQ->g(3,2);
    charQ->s(t,3,3);

    /* filling up the Q matrix */
    /////////////////////////////

    // Find eigenvalues and eigenvectors.
    //
    charU = new Db_matrix(char_as,char_as,"eigenvectors_1");
    charU->initialise();
    charV = new Db_matrix(char_as,char_as,"eigenvectors_2");
    charV->initialise();
    charRoot = new Db_matrix(char_as,"eigenvalues");
    charRoot->initialise();

    build_model(char_as,charPi,charQ,charU,charV,charRoot);

    if (Settings::noise > 4) {
        print_char_q_matrices(charQ);
    }

    delete charQ;
}


void Model_factory::protein_model(Settings *st)
{
    float ins_rate = 0.05;
    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    float del_rate = 0.05;
    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    float gap_ext = 0.5;
    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    this->protein_model(ins_rate,del_rate,gap_ext);
}

/*******************************************/

void Model_factory::protein_model(float ins_rate,float del_rate, float ext_prob)
{

    if (Settings::noise>4){
        cout<<"Protein substitution model: WAG\n";
    }

    char_ins_rate = ins_rate;
    char_del_rate = del_rate;
    char_ext_prob = ext_prob;


    if(Settings::noise > 4)
        print_char_alphabet();

    // pi and Q
    // Allocate space for P matrices
    //
    double tmp_pi[20] = {0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
                         0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956};

    charPi = new Db_matrix(char_fas,"pi_char");

    FOR(j,char_as) {
        charPi->s(tmp_pi[j],j);
    }

    double tmp_q[400] = {-1.0644447077525, 0.024253680012, 0.0199296524112, 0.0421562148098, 0.019829882912, 0.0333710782038, 0.091898529865, 0.117944490096, 0.0077435982602, 0.00937017411, 0.034303854235, 0.056214349179, 0.0174255844392, 0.0080896843586, 0.065832507505, 0.234330242141, 0.129414648097, 0.0016275200247, 0.008491734537, 0.142217282556,
                         0.0477814374309, -0.9283507809291, 0.0248352939324, 0.0084029714104, 0.0101982061898, 0.11148814755, 0.0254969723473, 0.048674413647, 0.052213352795, 0.009062124214, 0.042903719239, 0.331941090612, 0.0133235035374, 0.0039473788809, 0.0310955230559, 0.085103118001, 0.0338262340451, 0.016744036728, 0.0134582713486, 0.0178549859644,
                         0.0441670615592, 0.027937434312, -1.38391347672512, 0.309721806842, 0.0051215097968, 0.056694964284, 0.0549932739622, 0.093704896008, 0.096657307877, 0.026861601976, 0.011338897352, 0.186830763486, 0.0038658446967, 0.00369569221099, 0.0089275113111, 0.276280123717, 0.123859441762, 0.00103458645453, 0.0383077812, 0.0139129779176,
                         0.0640178448442, 0.006477251488, 0.212232770148, -0.94297329251968, 0.00058492787022, 0.0226532677023, 0.358464938024, 0.0720614260512, 0.0227376245588, 0.001911353642, 0.0073109283823, 0.029764733853, 0.0020234831358, 0.00179593805976, 0.0194028221904, 0.074506504504, 0.0228715867982, 0.0018668150853, 0.0114891949562, 0.010799881226,
                         0.088970318416, 0.023225614652, 0.0103686978864, 0.00172817559999, -0.46436390434772, 0.00362939371299, 0.0012396736328, 0.0255311625132, 0.0060827096236, 0.00824576291, 0.033128997983, 0.00459221916954, 0.0076154533014, 0.015296664838, 0.0050066661924, 0.097857567114, 0.0312985388968, 0.010315697313, 0.0191832740086, 0.071047316584,
                         0.0787099366842, 0.133477006, 0.060339961416, 0.0351844479133, 0.00190795624962, -1.31468853172694, 0.317551411783, 0.0274774230936, 0.104910689643, 0.005521101322, 0.074957777201, 0.24159519414, 0.030136742202, 0.00384014619352, 0.0427139961732, 0.071524881773, 0.0523445036856, 0.0031035709083, 0.008032288082, 0.0213594972636,
                         0.137118971515, 0.019310611604, 0.0370254015012, 0.352205574616, 0.0004122601456, 0.200883241107, -1.17853838814971, 0.0472634621406, 0.0139264517825, 0.00617432607, 0.013298858967, 0.160308574698, 0.0061457688348, 0.00311812993141, 0.0312266801005, 0.0490058789081, 0.0501991141155, 0.0022522133463, 0.0069244312826, 0.0417384374836,
                         0.122727478488, 0.02570888938, 0.043997465064, 0.0493773258384, 0.0059212002572, 0.0121221828612, 0.0329610245313, -0.4741383934945, 0.006093410533, 0.0014757945466, 0.0052849306733, 0.0231712797588, 0.00339542007, 0.0019189431989, 0.011146518267, 0.093280508578, 0.0137786810791, 0.0048478037397, 0.0036545482168, 0.0132749884132,
                         0.0274570594166, 0.0939747598, 0.154649002326, 0.0530905054876, 0.0048071015816, 0.157714501491, 0.0330950244725, 0.020763831438, -0.9455247927498, 0.00669751654, 0.043058119558, 0.0552322503552, 0.0078818406807, 0.0261095183349, 0.0318601786938, 0.0514549945251, 0.0288777379989, 0.0037772913771, 0.136632497248, 0.0083910614248,
                         0.0167482050465, 0.008221840588, 0.0216647526984, 0.0022496876087, 0.003284932553, 0.0041839549677, 0.0073964135655, 0.00253502563518, 0.003376161347, -1.17498318692916, 0.27336615273, 0.0200868455952, 0.083031965142, 0.040717445093, 0.00457305166728, 0.022206797976, 0.088966278632, 0.0030567591897, 0.014821160614, 0.55449575628,
                         0.0344705408285, 0.021883589212, 0.0051413506032, 0.00483769259197, 0.0074197365386, 0.0319346789409, 0.0089563400907, 0.00510364337166, 0.0122025059606, 0.15368423202, -0.69175870029473, 0.015975776073, 0.094666495854, 0.081290001923, 0.0190303105564, 0.0239655313281, 0.0199280900994, 0.0095710687431, 0.0140609310556, 0.127636184504,
                         0.0785078337935, 0.23531264024, 0.117737663694, 0.0273733764605, 0.00142943173442, 0.14305227669, 0.150049162927, 0.0310993759044, 0.0217544113216, 0.015694841712, 0.022203558995, -1.07152398375532, 0.0182209045452, 0.0034141362684, 0.0254852873376, 0.067232846627, 0.084623394646, 0.0019781331795, 0.0047007809888, 0.0216539266904,
                         0.0774016821384, 0.030039999464, 0.0077483399574, 0.0059186573054, 0.0075393483596, 0.056754463806, 0.0182957528036, 0.01449413838, 0.0098736900133, 0.20634205636, 0.41846021018, 0.0579518322936, -1.2597234186325, 0.045758173097, 0.0078405461599, 0.0343352383995, 0.092502574724, 0.0074188949454, 0.0151127724254, 0.14593504782,
                         0.0182346531826, 0.004516408092, 0.00375891879174, 0.00266574034104, 0.007684890556, 0.00366990113448, 0.00471054498671, 0.0041568456258, 0.0165979167123, 0.05134827302, 0.18234669053, 0.0055103727096, 0.023220499701, -0.67999937105987, 0.0073881779164, 0.0379519766649, 0.0104882661681, 0.022005248076, 0.227669563576, 0.0460744832752,
                         0.124618565545, 0.029878490308, 0.0076255992414, 0.0241862096784, 0.0021123505512, 0.0342809801532, 0.0396167807095, 0.020277640926, 0.0170090221974, 0.0048431492208, 0.035849495396, 0.0345434792256, 0.0033413780883, 0.0062045996636, -0.5770185229931, 0.112151837712, 0.0485285253768, 0.0020054663895, 0.0076208498132, 0.0223241027972,
                         0.292004459041, 0.05383008268, 0.155350266162, 0.061138656376, 0.027178817748, 0.037788440247, 0.0409279829071, 0.111708930276, 0.0180832908897, 0.01548197904, 0.029719604451, 0.059989719918, 0.0096324810435, 0.0209811655989, 0.073828693968, -1.326554610767, 0.267114820854, 0.0075345000378, 0.0277605484806, 0.0165001710484,
                         0.183747304969, 0.024378648436, 0.079353827364, 0.0213842684566, 0.0099045924752, 0.0315100653768, 0.0477688308585, 0.0188010037494, 0.0115635053091, 0.07067118256, 0.028157755998, 0.086032427628, 0.029568433524, 0.0066065589057, 0.0363992375304, 0.304350756558, -1.1004826896859, 0.0015948784176, 0.0102700127816, 0.098419398788,
                         0.0098004742107, 0.05117989024, 0.00281118065298, 0.0074025714917, 0.013845044146, 0.0079236101097, 0.0090895272073, 0.0280544413194, 0.0064149020097, 0.010298201078, 0.057355623581, 0.008529242643, 0.0100576594062, 0.058786971516, 0.0063796049555, 0.0364094439818, 0.0067641119728, -0.44467569893618, 0.087670143938, 0.0259030544764,
                         0.0208543675065, 0.016776769076, 0.0424510884, 0.0185802165661, 0.0105002187974, 0.008363355651, 0.0113971362467, 0.0086252194872, 0.094633174672, 0.02036395922, 0.034364459162, 0.0082661793504, 0.0083556782799, 0.248050243532, 0.0098869347026, 0.0547101006747, 0.0177637255796, 0.035754572001, -0.6920103710931, 0.022312972188,
                         0.173776433679, 0.011074304228, 0.0076711383924, 0.0086899653085, 0.019349118692, 0.0110654786961, 0.0341810742559, 0.0155886497946, 0.0028916398054, 0.3790671258, 0.15520551106, 0.0189456434124, 0.040145332815, 0.0249765843548, 0.0144102052697, 0.0161795265281, 0.084699660521, 0.0052561618971, 0.011101848966, -1.034275403476 };

    Db_matrix *charQ  = new Db_matrix(char_as,char_as,"Q_char");

    FOR(j,char_as) {
        FOR(i,char_as) {
            charQ->s(tmp_q[j*char_as+i],j,i);
        }
    }

    // Find eigenvalues and eigenvectors.
    //
    charU = new Db_matrix(char_as,char_as,"eigenvectors_1");
    charU->initialise();
    charV = new Db_matrix(char_as,char_as,"eigenvectors_2");
    charV->initialise();
    charRoot = new Db_matrix(char_as,"eigenvalues");
    charRoot->initialise();

    build_model(char_as,charPi,charQ,charU,charV,charRoot);

    if (Settings::noise > 4) {
        print_char_q_matrices(charQ);
    }

    delete charQ;
}

/*******************************************/

void Model_factory::build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot)
{


    // Find eigenvalues and eigenvectors.
    //

    Eigen* e = new Eigen();

    double* tpi = new double[s];
    double* tsq = new double[s];
    double* tcq = new double[s*s];
    double* twr = new double[s];
    double* twu = new double[s*s];
    double* twv = new double[s*s];
    int npi0 = 0;

    FOR(i,s) {

        tpi[i] = pi->g(i);

        FOR(j,s) {
            tcq[i*s+j] = q->g(i,j);
        }
    }

    if (e->getpi_sqrt (tpi, tsq, s, &npi0) != 0) {
        cout << "\nError in eigen square roots!!\n\n";
        exit (-1);
    }


    if (e->eigenQREV (tcq, tpi, tsq, s, npi0, twr, twu, twv) != 0) {
        cout << "\nError in eigen QREV!!\n\n";
        exit (-1);
    }

    FOR(i,s) {

        wRoot->s(twr[i],i);

        FOR(j,s) {
            wU->s(twu[i*s+j],i,j);
            wV->s(twv[i*s+j],i,j);
        }
    }

    delete []tpi;
    delete []tsq;
    delete []tcq;
    delete []twr;
    delete []twu;
    delete []twv;

    delete e;

}

/*******************************************/

Evol_model Model_factory::alignment_model(double distance)
{

    // Compute the P matrix for regular DNA alphabet (four bases).
    //
    Eigen* e = new Eigen();

    double* tmr = new double[char_as*char_as];
    double* twr = new double[char_as];
    double* twu = new double[char_as*char_as];
    double* twv = new double[char_as*char_as];

    FOR(i,char_as) {

        twr[i] = charRoot->g(i);

        FOR(j,char_as) {
            twu[i*char_as+j] = charU->g(i,j);
            twv[i*char_as+j] = charV->g(i,j);
        }
    }

    e->computePMatrix(char_as,tmr,twu,twv,twr,distance);

    Evol_model model(full_char_alphabet, distance);

    model.ins_rate = char_ins_rate;
    model.del_rate = char_del_rate;

    model.ins_prob = (1.0-exp(-1.0*char_ins_rate*distance));
    model.del_prob = (1.0-exp(-1.0*char_del_rate*distance));

    double t = (1.0-exp(-1.0*(char_ins_rate+char_del_rate)*distance));

    t /= 2.0;
    model.log_id_prob = log(t);
    model.log_match_prob = log(1.0-2*t);
    model.log_ext_prob = log(char_ext_prob);

    model.id_prob = t;
    model.match_prob = 1.0-2*t;
    model.ext_prob = char_ext_prob;

cout<<"model.ins_prob "<<model.ins_prob<<"\n";
cout<<"model.del_prob "<<model.del_prob<<"\n";
cout<<"model.match_prob "<<model.match_prob<<"\n";
cout<<"sum "<<model.ins_prob+model.del_prob+model.match_prob<<"\n";
cout<<"model.ext_prob "<<model.ext_prob<<"\n";

    FOR(i,char_as) {

        model.charPi->s(charPi->g(i),i);
        model.logCharPi->s(log( charPi->g(i) ),i);

        FOR(j,char_as) {

            model.charPr->s(tmr[i*char_as+j],i,j);
            model.logCharPr->s(log( tmr[i*char_as+j] ),i,j);

        }
    }

    delete[] tmr;
    delete[] twr;
    delete[] twu;
    delete[] twv;

    delete e;


    if(sequence_data_type == Model_factory::dna)
    {
        // Extend the P matrix for ambiguity DNA alphabet (16 bases).
        //
        Db_matrix *char_ambiguity = new Db_matrix(4,char_fas);
        char_ambiguity->initialise(0);

        float ambiguity_factor = 1.0;
        if( Settings_handle::st.is("ambiguity-factor") )
            ambiguity_factor = Settings_handle::st.get("ambiguity-factor").as<float>();

        if( ambiguity_factor > 1.0 || ambiguity_factor < 0 )
            ambiguity_factor = 1.0;

        for(unsigned int ai=0;ai<char_letters.size();ai++)
        {
            Char_symbol a = char_letters.at(ai);

            float probability = pow(ambiguity_factor, a.n_residues);

            for(int aj=0;aj<a.n_residues;aj++){
                int at = char_alphabet.find(a.residues.at(aj));
                char_ambiguity->s(probability, at, ai);
            }
        }


        int n,m;

        FOR(i,char_fas) {
            FOR(j,char_fas) {

                model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                if(i<char_as && j<char_as)
                    continue;

                double max = 0;

                FOR(n,char_as) {
                    FOR(m,char_as) {
                        double t = model.charPr->g(n,m)*char_ambiguity->g(m,j)*char_ambiguity->g(n,i);
                        if(max < t) max = t;
                    }
                }

                model.charPr->s(max,i,j);
                model.logCharPr->s(log(max),i,j);
            }
        }

        delete char_ambiguity;
    }

    else if(sequence_data_type == Model_factory::protein)
    {

        // Extend the P matrix for ambiguity protein alphabet (currently 51 residues).
        //
        Db_matrix *char_ambiguity = new Db_matrix(20,char_fas);
        char_ambiguity->initialise(0);

        float ambiguity_factor = 1.0;
        if( Settings_handle::st.is("ambiguity-factor") )
            ambiguity_factor = Settings_handle::st.get("ambiguity-factor").as<float>();

        if( ambiguity_factor > 1.0 || ambiguity_factor < 0 )
            ambiguity_factor = 1.0;

        for(unsigned int ai=0;ai<char_letters.size();ai++)
        {
            Char_symbol a = char_letters.at(ai);

            float probability = pow(ambiguity_factor, a.n_residues);

            for(int aj=0;aj<a.n_residues;aj++){
                int at = char_alphabet.find(a.residues.at(aj));
                char_ambiguity->s(probability, at, ai);
            }
        }


        int n,m;

        FOR(i,char_fas) {
            FOR(j,char_fas) {

                model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                if(i<char_as && j<char_as)
                    continue;

                double max = 0;

                FOR(n,char_as) {
                    FOR(m,char_as) {
                        double t = model.charPr->g(n,m)*char_ambiguity->g(m,j)*char_ambiguity->g(n,i);
                        if(max < t) max = t;
                    }
                }

                model.charPr->s(max,i,j);
                model.logCharPr->s(log(max),i,j);
            }
        }

        delete char_ambiguity;

    }

    else if(sequence_data_type == Model_factory::codon)
    {

    }

    if (Settings::noise > 4) {
        print_char_p_matrices(model);
    }

    return model;
}

/*******************************************/

void Model_factory::print_char_p_matrices(Evol_model &model)
{
    // Print out the model
    cout<<"\nModel_factory::print_char_p_matrices()\n\n";
    cout<<"alphabet "<<char_alphabet<<endl;
    cout<<"distance "<<model.distance<<endl;

    cout<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,char_as)
        cout << setw(4) <<char_alphabet.at(i)<<"   ";
    cout<<endl;

    FOR(j,char_as) {
        cout<<" "<<model.charPi->g(j);
    }
    cout<<endl;

    cout<<"\nsubstitution matrix"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,char_as)
        cout << setw(6) <<char_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,char_as) {
        cout<<char_alphabet.at(i)<<" ";
        FOR(j,char_as) {
            cout<<" "<<model.charPr->g(i,j);
        }
        cout<<endl;
    }
    cout<<"\nlog substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,char_as)
        cout << setw(7) <<char_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,char_as) {
        cout<<char_alphabet.at(i)<<" ";
        FOR(j,char_as) {
            cout<<" "<<model.logCharPr->g(i,j);
        }
        cout<<endl;
    }


    cout << noshowpos <<"\nfull alphabet "<<full_char_alphabet<<endl;

    cout<<"\nsubstitution matrix"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,char_fas)
        cout << setw(6) <<full_char_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,char_fas) {
        cout<<full_char_alphabet.at(i)<<" ";
        FOR(j,char_fas) {
            cout<<" "<<model.charPr->g(i,j);
        }
        cout<<endl;
    }


    cout<<"\nlog substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,char_fas)
        cout << setw(7) <<full_char_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,char_fas) {
        cout<<full_char_alphabet.at(i)<<" ";
        FOR(j,char_fas) {
            cout<<" "<<model.logCharPr->g(i,j);
        }
        cout<<endl;
    }

    cout<<endl;
    cout<<"indel prob:     "<<model.id_prob<<", "<<model.log_id_prob<<endl;
    cout<<"extension prob: "<<model.ext_prob<<", "<<model.log_ext_prob<<endl;
    cout<<"match prob:     "<<model.match_prob<<", "<<model.log_match_prob<<endl;;
    cout<<endl;

}

void Model_factory::print_char_q_matrices(Db_matrix *charQ)
{
    // Print out the model
    cout<<"\nModel_factory::print_char_q_matrices()\n\n";
    cout<<"alphabet "<<char_alphabet<<endl;

    cout<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,char_as)
        cout << setw(4) <<char_alphabet.at(i)<<"   ";
    cout<<endl;

    FOR(j,char_as) {
        cout<<" "<<charPi->g(j);
    }
    cout<<endl;

    cout<<"\noriginal substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,char_as)
        cout << setw(7) <<char_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,char_as) {
     cout<<full_char_alphabet.at(i)<<" ";
     FOR(j,char_as) {
            cout<<" "<<charQ->g(i,j);
        }
        cout<<endl;
    }

    cout<<"\neigen values & vectors"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(j,char_as) {
        cout<<" "<<charRoot->g(j);
    }
    cout<<endl<<endl;

    FOR(i,char_as) {
        FOR(j,char_as) {
            cout<<" "<<charU->g(i,j);
        }
        cout<<endl;
    }
    cout<<endl;

    FOR(i,char_as) {
        FOR(j,char_as) {
            cout<<" "<<charV->g(i,j);
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"ins rate:     "<<char_ins_rate<<endl;
    cout<<"del rate:     "<<char_del_rate<<endl;
    cout<<"extension prob: "<<char_ext_prob<<endl;
    cout<<endl;

}



