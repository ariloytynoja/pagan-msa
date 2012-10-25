/***************************************************************************
 *   Copyright (C) 2010-2012 by Ari Loytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
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

#include "utils/mafft_alignment.h"
#include "utils/log_output.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>

#if defined (__APPLE__)
#include <mach-o/dyld.h>
#endif

using namespace std;
using namespace ppa;

Mafft_alignment::Mafft_alignment()
{
}

bool Mafft_alignment::test_executable()
{
    #if defined (__CYGWIN__)
    char path[200];
    int length = readlink("/proc/self/exe",path,200-1);

    string epath = string(path).substr(0,length);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));
    mafftpath = epath;
    epath = epath+"sh.exe "+epath+"mafft -h >/dev/null 2>/dev/null";
    int status = system(epath.c_str());
    return WEXITSTATUS(status) == 1;

    # else
    int status = system("mafft -h >/dev/null 2>/dev/null");

    if(WEXITSTATUS(status) == 1)
        return true;

    char path[200];
    string epath;

    #if defined (__APPLE__)
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    epath = string(path);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));
    //epath = "DYLD_LIBRARY_PATH="+epath+" "+epath;

    #else
    int length = readlink("/proc/self/exe",path,200-1);
    epath = string(path).substr(0,length);
    epath.replace(epath.rfind("pagan"),string("pagan").size(),string(""));

    #endif

    mafftpath = epath;
    epath = epath+"mafft -h >/dev/null 2>/dev/null";
    status = system(epath.c_str());
    return WEXITSTATUS(status) == 1;

    #endif
}

void Mafft_alignment::align_sequences(vector<Fasta_entry> *sequences)
{
    ofstream m_output;
    string tmp_dir = this->get_temp_dir();

    int r = rand();
    while(true)
    {
        stringstream m_name;
        m_name <<tmp_dir<<"m"<<r<<".fas";
        ifstream m_file(m_name.str().c_str());

        if(!m_file)
        {
            m_output.open( m_name.str().c_str(), (ios::out) );
            break;
        }
        r = rand();
    }

    map<string,string> dna_seqs;
    bool has_dna_seqs = (sequences->at(0).sequence.length() > 0);

    vector<Fasta_entry>::iterator si = sequences->begin();
    for(;si!=sequences->end();si++)
    {
        m_output<<">"<<si->name<<endl<<si->sequence<<endl;
        if(has_dna_seqs)
            dna_seqs.insert(pair<string,string>(si->name,si->dna_sequence));
    }
    m_output.close();
    sequences->clear();

    stringstream command;
    command << mafftpath<<"mafft "+tmp_dir+"m"<<r<<".fas  2>/dev/null";


    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with mafft pipe.\nExiting.\n",0);
        exit(1);
    }


    // read mafft output
    string name, sequence = "";  // Initialization
    char temp[256];

    while ( fgets( temp, sizeof temp, fpipe))
    {
    	string line(temp);

        if (line[0] == '>')
        {
            line = this->remove_last_whitespaces(line);

            // If a name and a sequence were found
            if ((name != "") && (sequence != ""))
            {
                Fasta_entry s;
                s.name = name;
                sequence = this->remove_whitespaces(sequence);
                transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
                s.sequence = sequence;

                if(has_dna_seqs)
                {
                    map<string,string>::iterator it = dna_seqs.find(name);
                    if(it!=dna_seqs.end())
                        s.dna_sequence = it->second;
                }

                sequences->push_back(s);
                name = "";
                sequence = "";
            }
            name = line;
            name.erase(name.begin());  // Character > deletion
        }
        else
        {
            sequence += temp;  // Sequence isolation
        }
    }

    // Addition of the last sequence in file
    if ((name != "") && (sequence != ""))
    {
        Fasta_entry s;
        s.name = name;
        sequence = this->remove_whitespaces(sequence);
        transform( sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int))toupper );
        s.sequence = sequence;

        if(has_dna_seqs)
        {
            map<string,string>::iterator it = dna_seqs.find(name);
            if(it!=dna_seqs.end())
                s.dna_sequence = it->second;
        }

        sequences->push_back(s);
    }

    pclose(fpipe);

    if(!Settings_handle::st.is("keep-temp-files"))
        this->delete_files(r);
}

void Mafft_alignment::delete_files(int r)
{

    string tmp_dir = this->get_temp_dir();

    stringstream m_name;
    m_name <<tmp_dir<<"m"<<r<<".fas";

    if ( remove( m_name.str().c_str() ) != 0 )
        Log_output::write_out( "Error deleting file", 1);
}
