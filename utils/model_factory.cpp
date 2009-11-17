
#include <iomanip>
#include <iostream>
#include <vector>
#include "model_factory.h"
#include "settings.h"
#include "settings_handle.h"
#include "eigen.h"

using namespace ppa;
using namespace std;

Model_factory::Model_factory()
{

    define_dna_alphabet();

}

Model_factory::~Model_factory()
{
    if(dnaPi!=0)
        delete dnaPi;

    if(dnaU!=0)
        delete dnaU;

    if(dnaV!=0)
        delete dnaV;

    if(dnaRoot!=0)
        delete dnaRoot;

    if(parsimony_table!=0)
        delete parsimony_table;
}

/*******************************************/

/*
 * Definition of DNA alpahbet including ambiguity characters.
 */
void Model_factory::define_dna_alphabet()
{
    dna_alphabet = "ACGT";
    full_dna_alphabet = "ACGTRYMKWSBDHVN";

    dna_as = dna_alphabet.length();
    dna_fas = full_dna_alphabet.length();

    int n_bases[] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};
    string ambiguity[] = {"A","C","G","T","AG","CT","AC","GT","AT","CG","CGT","AGT","ACT","ACG","ACGT"};

    for(int i=0;i<15;i++)
    {
        Dna_symbol letter;

        letter.index = i;
        letter.symbol = full_dna_alphabet.at(i);
        letter.n_bases = n_bases[i];
        for(int j=0;j<letter.n_bases;j++)
            letter.bases.push_back(ambiguity[i].at(j));
        dna_letters.push_back(letter);
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

    parsimony_table = new Int_matrix(dna_fas,dna_fas,"parsimony_dna");

    for(int i=0;i<dna_fas;i++)
    {
        for(int j=0;j<dna_fas;j++)
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
            cout<<full_dna_alphabet.at(i)<<" ";
        cout<<endl;

        for(int i=0;i<15;i++)
        {
            cout<<full_dna_alphabet.at(i)<<" ";
            for(int j=0;j<15;j++)
            {
                cout<<full_dna_alphabet.at(parsimony_table->g(i,j))<<" ";
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
void Model_factory::print_dna_alphabet()
{
    cout<<"\nModel_factory::print_dna_alphabet()\n";
    for(unsigned int i=0;i<dna_letters.size();i++)
    {
        Dna_symbol a = dna_letters.at(i);

        cout<<"Index "<<a.index<<"; symbol "<<a.symbol<<"; bases ("<<a.n_bases<<"): ";
        for(int j=0;j<a.n_bases;j++){
            int t = dna_alphabet.find(a.bases.at(j));
            cout<<a.bases.at(j)<<"["<<t<<"] ";
        }
        cout<<endl;
    }
}

/*******************************************/

/*
 * Definition of DNA model: Q matrix and indel rate/prob.
 */

void Model_factory::dna_model(float *dna_pi,Settings *st)
{
    float dna_kappa = 2.0;
    if(st->is("dna-kappa"))
        dna_kappa =  st->get("dna-kappa").as<float>();

    float dna_rho = 1.0;
    if(st->is("dna-rho"))
        dna_rho =  st->get("dna-rho").as<float>();

    float ins_rate = 0.05;
    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    float del_rate = 0.05;
    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    float gap_ext = 0.5;
    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    this->dna_model(dna_pi,dna_kappa,dna_rho,ins_rate,del_rate,gap_ext);
}

void Model_factory::dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob)
{

    if (Settings::noise>4){
        cout<<"DNA substitution model: base frequencies "
                <<pi[0]<<", "<<pi[1]<<", "<<pi[2]<<" and "<<pi[3]
                <<"; kappa "<<kappa<<" and rho "<<rho<<"; ins rate "<<ins_rate<<", del rate "<<del_rate<<" and ext. prob "<<ext_prob<<"."<<endl;
    }

    dna_ins_rate = ins_rate;
    dna_del_rate = del_rate;
    dna_ext_prob = ext_prob;


    if(Settings::noise > 4)
        print_dna_alphabet();

    // pi and Q
    // Allocate space for P matrices
    //
    Db_matrix *dnaQ  = new Db_matrix(dna_as,dna_as,"Q_dna");

    dnaPi = new Db_matrix(dna_fas,"pi_dna");

    FOR(j,dna_as) {
        dnaPi->s(pi[j],j);
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
    dnaQ->s(t,0,1);

    t = alfaR*pi[2] / piR+beta*pi[2];
    dnaQ->s(t,0,2);

    t = beta*pi[3];
    dnaQ->s(t,0,3);

    t = 0-dnaQ->g(0,1)-dnaQ->g(0,2)-dnaQ->g(0,3);
    dnaQ->s(t,0,0);

    /*2nd row*/
    t = beta*pi[0];
    dnaQ->s(t,1,0);

    t = beta*pi[2];
    dnaQ->s(t,1,2);

    t = alfaY*pi[3]/piY+beta*pi[3];
    dnaQ->s(t,1,3);

    t = 0-dnaQ->g(1,0)-dnaQ->g(1,2)-dnaQ->g(1,3);
    dnaQ->s(t,1,1);

    /*3rd row*/
    t = alfaR*pi[0]/piR+beta*pi[0];
    dnaQ->s(t,2,0);

    t = beta*pi[1];
    dnaQ->s(t,2,1);

    t = beta*pi[3];
    dnaQ->s(t,2,3);

    t = 0-dnaQ->g(2,0)-dnaQ->g(2,1)-dnaQ->g(2,3);
    dnaQ->s(t,2,2);

    /*4th row*/
    t = beta*pi[0];
    dnaQ->s(t,3,0);

    t = alfaY*pi[1]/piY+beta*pi[1];
    dnaQ->s(t,3,1);

    t = beta*pi[2];
    dnaQ->s(t,3,2);

    t = 0-dnaQ->g(3,0)-dnaQ->g(3,1)-dnaQ->g(3,2);
    dnaQ->s(t,3,3);

    /* filling up the Q matrix */
    /////////////////////////////

    // Find eigenvalues and eigenvectors.
    //
    dnaU = new Db_matrix(dna_as,dna_as,"eigenvectors_1");
    dnaU->initialise();
    dnaV = new Db_matrix(dna_as,dna_as,"eigenvectors_2");
    dnaV->initialise();
    dnaRoot = new Db_matrix(dna_as,"eigenvalues");
    dnaRoot->initialise();

    build_model(dna_as,dnaPi,dnaQ,dnaU,dnaV,dnaRoot);

    if (Settings::noise > 4) {
        print_dna_q_matrices(dnaQ);
    }

    delete dnaQ;
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

Dna_model Model_factory::dna_alignment_model(double distance)
{

    // Compute the P matrix for regular DNA alphabet (four bases).
    //
    Eigen* e = new Eigen();

    double* tmr = new double[dna_as*dna_as];
    double* twr = new double[dna_as];
    double* twu = new double[dna_as*dna_as];
    double* twv = new double[dna_as*dna_as];

    FOR(i,dna_as) {

        twr[i] = dnaRoot->g(i);

        FOR(j,dna_as) {
            twu[i*dna_as+j] = dnaU->g(i,j);
            twv[i*dna_as+j] = dnaV->g(i,j);
        }
    }

    e->computePMatrix(dna_as,tmr,twu,twv,twr,distance);

    Dna_model model(full_dna_alphabet, distance);

    model.ins_rate = dna_ins_rate;
    model.del_rate = dna_del_rate;

    model.ins_prob = (1.0-exp(-1.0*dna_ins_rate*distance));
    model.del_prob = (1.0-exp(-1.0*dna_del_rate*distance));

    double t = (1.0-exp(-1.0*(dna_ins_rate+dna_del_rate)*distance));

    t /= 2.0;
    model.log_id_prob = log(t);
    model.log_match_prob = log(1.0-2*t);
    model.log_ext_prob = log(dna_ext_prob);

    model.id_prob = t;
    model.match_prob = 1.0-2*t;
    model.ext_prob = dna_ext_prob;

    FOR(i,dna_as) {

        model.dnaPi->s(dnaPi->g(i),i);
        model.logDnaPi->s(log( dnaPi->g(i) ),i);

        FOR(j,dna_as) {

            model.dnaPr->s(tmr[i*dna_as+j],i,j);
            model.logDnaPr->s(log( tmr[i*dna_as+j] ),i,j);

        }
    }

    delete[] tmr;
    delete[] twr;
    delete[] twu;
    delete[] twv;

    delete e;

    // Extend the P matrix for ambiguity DNA alphabet (16 bases).
    //
    Db_matrix *dna_ambiguity = new Db_matrix(4,dna_fas);
    dna_ambiguity->initialise(0);

    float ambiguity_factor = 1.0;
    if( Settings_handle::st.is("ambiguity-factor") )
        ambiguity_factor = Settings_handle::st.get("ambiguity-factor").as<float>();
    if( ambiguity_factor > 1.0 || ambiguity_factor < 0 )
        ambiguity_factor = 1.0;

    for(unsigned int ai=0;ai<dna_letters.size();ai++)
    {
        Dna_symbol a = dna_letters.at(ai);

        float probability = pow(ambiguity_factor, a.n_bases);

        for(int aj=0;aj<a.n_bases;aj++){
            int at = dna_alphabet.find(a.bases.at(aj));
            dna_ambiguity->s(probability, at, ai);
        }
    }


    int n,m;

    FOR(i,dna_fas) {
        FOR(j,dna_fas) {

            model.parsimony_table->s( parsimony_table->g(i,j), i,j);

            if(i<dna_as && j<dna_as)
                continue;

            double max = 0;

            FOR(n,dna_as) {
                FOR(m,dna_as) {
                    double t = model.dnaPr->g(n,m)*dna_ambiguity->g(m,j)*dna_ambiguity->g(n,i);
                    if(max < t) max = t;
                }
            }

            model.dnaPr->s(max,i,j);
            model.logDnaPr->s(log(max),i,j);

        }
    }

    delete dna_ambiguity;

    if (Settings::noise > 4) {
        print_dna_p_matrices(model);
    }

    return model;
}

/*******************************************/

void Model_factory::print_dna_p_matrices(Dna_model &model)
{
    // Print out the model
    cout<<"\nModel_factory::print_dna_p_matrices()\n\n";
    cout<<"alphabet "<<dna_alphabet<<endl;
    cout<<"distance "<<model.distance<<endl;

    cout<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,dna_as)
        cout << setw(4) <<dna_alphabet.at(i)<<"   ";
    cout<<endl;

    FOR(j,dna_as) {
        cout<<" "<<model.dnaPi->g(j);
    }
    cout<<endl;

    cout<<"\nsubstitution matrix"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,dna_as)
        cout << setw(6) <<dna_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,dna_as) {
        cout<<dna_alphabet.at(i)<<" ";
        FOR(j,dna_as) {
            cout<<" "<<model.dnaPr->g(i,j);
        }
        cout<<endl;
    }
    cout<<"\nlog substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,dna_as)
        cout << setw(7) <<dna_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,dna_as) {
        cout<<dna_alphabet.at(i)<<" ";
        FOR(j,dna_as) {
            cout<<" "<<model.logDnaPr->g(i,j);
        }
        cout<<endl;
    }


    cout << noshowpos <<"\nfull alphabet "<<full_dna_alphabet<<endl;

    cout<<"\nsubstitution matrix"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,dna_fas)
        cout << setw(6) <<full_dna_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,dna_fas) {
        cout<<full_dna_alphabet.at(i)<<" ";
        FOR(j,dna_fas) {
            cout<<" "<<model.dnaPr->g(i,j);
        }
        cout<<endl;
    }


    cout<<"\nlog substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,dna_fas)
        cout << setw(7) <<full_dna_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,dna_fas) {
        cout<<full_dna_alphabet.at(i)<<" ";
        FOR(j,dna_fas) {
            cout<<" "<<model.logDnaPr->g(i,j);
        }
        cout<<endl;
    }

    cout<<endl;
    cout<<"indel prob:     "<<model.id_prob<<", "<<model.log_id_prob<<endl;
    cout<<"extension prob: "<<model.ext_prob<<", "<<model.log_ext_prob<<endl;
    cout<<"match prob:     "<<model.match_prob<<", "<<model.log_match_prob<<endl;;
    cout<<endl;

}

void Model_factory::print_dna_q_matrices(Db_matrix *dnaQ)
{
    // Print out the model
    cout<<"\nModel_factory::print_dna_q_matrices()\n\n";
    cout<<"alphabet "<<dna_alphabet<<endl;

    cout<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    cout << fixed << noshowpos << setprecision (4);

    FOR(i,dna_as)
        cout << setw(4) <<dna_alphabet.at(i)<<"   ";
    cout<<endl;

    FOR(j,dna_as) {
        cout<<" "<<dnaPi->g(j);
    }
    cout<<endl;

    cout<<"\noriginal substitution matrix"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(i,dna_as)
        cout << setw(7) <<dna_alphabet.at(i)<<" ";
    cout<<endl;

    FOR(i,dna_as) {
     cout<<full_dna_alphabet.at(i)<<" ";
     FOR(j,dna_as) {
            cout<<" "<<dnaQ->g(i,j);
        }
        cout<<endl;
    }

    cout<<"\neigen values & vectors"<<endl;
    cout << fixed << showpos << setprecision (4);

    FOR(j,dna_as) {
        cout<<" "<<dnaRoot->g(j);
    }
    cout<<endl<<endl;

    FOR(i,dna_as) {
        FOR(j,dna_as) {
            cout<<" "<<dnaU->g(i,j);
        }
        cout<<endl;
    }
    cout<<endl;

    FOR(i,dna_as) {
        FOR(j,dna_as) {
            cout<<" "<<dnaV->g(i,j);
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"ins rate:     "<<dna_ins_rate<<endl;
    cout<<"del rate:     "<<dna_del_rate<<endl;
    cout<<"extension prob: "<<dna_ext_prob<<endl;
    cout<<endl;

}



