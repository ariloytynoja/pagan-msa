/***************************************************************************
 *   Copyright (C) 2013 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>
#include "utils/bppancestors.h"
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "utils/settings_handle.h"

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;
using namespace ppa;

BppAncestors::BppAncestors() {}

bool BppAncestors::test_executable()
{

    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    return WEXITSTATUS(status) == 0;

    # else

    char path[200];
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);
    //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    epath = string(path).substr(0,length);
    if (epath.find("/")!=std::string::npos)
        epath = epath.substr(0,epath.rfind("/")+1);

    #endif

    bppdistpath = epath;
    epath = epath+"bppancestor >/dev/null 2>/dev/null";
    int status = system(epath.c_str());

    if(WEXITSTATUS(status) == 0)
        return true;

    bppdistpath = "";
    status = system("bppancestor >/dev/null 2>/dev/null");

    return WEXITSTATUS(status) == 0;

    #endif
}

void BppAncestors::infer_ancestors(vector<Fasta_entry> *aligned_sequences,string tree,bool isDna)
{

    string tmp_dir = this->get_temp_dir();

    stringstream f_name;
    stringstream t_name;
    stringstream o_name;

    int r = rand();
    while(true)
    {

        f_name.str("");
        t_name.str("");
        o_name.str("");

        f_name <<tmp_dir<<"f"<<r<<".fas";
        ifstream f_file(f_name.str().c_str());

        t_name <<tmp_dir<<"t"<<r<<".tre";
        ifstream t_file(t_name.str().c_str());

        o_name <<tmp_dir<<"o"<<r<<".fas";
        ifstream o_file(t_name.str().c_str());

        if(!f_file && !t_file && !o_file)
        {
            ofstream f_tmp;
            f_tmp.open(f_name.str().c_str(), (ios::out) );
            ofstream t_tmp;
            t_tmp.open(t_name.str().c_str(), (ios::out) );
            ofstream o_tmp;
            o_tmp.open(o_name.str().c_str(), (ios::out) );

            break;
        }
        r = rand();
    }



    ////////////



    ofstream f_output;
    f_output.open( f_name.str().c_str(), (ios::out) );

    vector<Fasta_entry>::iterator si = aligned_sequences->begin();
    for(;si!=aligned_sequences->end();si++)
        f_output<<">"<<si->name<<"\n"<<si->sequence<<"\n";
    f_output.close();



    ofstream t_output;
    t_output.open( t_name.str().c_str(), (ios::out) );
    t_output<<tree<<endl;
    t_output.close();


    stringstream command;
    command << bppdistpath<<"bppancestor input.sequence.file="<<f_name.str()<<" input.sequence.format=Fasta input.sequence.sites_to_use=all input.tree.file="<<t_name.str()<<
            " input.tree.format=NHX input.sequence.max_gap_allowed=100% initFreqs=observed output.sequence.file="<<o_name.str()<<" output.sequence.format=Phylip";
    if(!isDna)
        command << " alphabet=Protein model=WAG01";
    else
    {
        if(Settings_handle::st.is("codons"))
            command << " alphabet=Codon\\(letter=DNA,type=Standard\\) model=YN98\\(kappa=2,omega=0.5\\)";
        else
            command << " alphabet=DNA model=HKY85";
    }

    Log_output::write_out("BppAncestors: command: "+command.str()+"\n",2);


    FILE *fpipe;
    char line[256];

    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with bppancestor pipe.\nExiting.\n",0);
        exit(1);
    }
    while ( fgets( line, sizeof line, fpipe))
    {
        Log_output::write_out("BppAncestors: command: "+string(line)+"\n",3);
    }
    pclose(fpipe);

    map<string,string> bppa_sequences;

    Fasta_reader fr;
    fr.read_bpp_phylip(o_name.str().c_str(),&bppa_sequences);

    si = aligned_sequences->begin();
    for(;si!=aligned_sequences->end();si++)
    {
        map<string,string>::iterator fi = bppa_sequences.find(si->name);
        if(fi != bppa_sequences.end())
        {
            for(int i=0;i<si->sequence.length();i++)
            {
                if( si->sequence.at(i)!='-' && si->sequence.at(i)!='.' )
                    si->sequence.at(i) = fi->second.at(i);
            }
        }
    }

    this->delete_files(r);

}

void BppAncestors::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream t_name;
    t_name <<tmp_dir<<"t"<<r<<".tre";

    stringstream f_name;
    f_name <<tmp_dir<<"f"<<r<<".fas";

    stringstream o_name;
    o_name <<tmp_dir<<"o"<<r<<".fas";

    if ( remove( t_name.str().c_str() ) != 0 )
        perror( "Error deleting file" );
    if ( remove( f_name.str().c_str() ) != 0 )
        perror( "Error deleting file");
    if ( remove( o_name.str().c_str() ) != 0 )
        perror( "Error deleting file");

}
