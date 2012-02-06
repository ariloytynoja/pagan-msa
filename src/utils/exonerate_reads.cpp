#include "exonerate_reads.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <boost/regex.hpp>

using namespace std;
using namespace ppa;

Exonerate_reads::Exonerate_reads()
{

}

bool Exonerate_reads::test_executable()
{
    int status = system("exonerate  >/dev/null");
    return WEXITSTATUS(status) == 1;
}

bool Exonerate_reads::split_sugar_string(const string& row,hit *h)
{

    const boost::regex pattern("sugar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );
    }

    return valid;
}

bool Exonerate_reads::split_vulgar_string(const string& row,hit *h)
{

    const boost::regex pattern("vulgar:\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S)\\s+(\\d+)\\s+(.+)\n");
    boost::match_results<string::const_iterator> result;
    bool valid = boost::regex_match(row, result, pattern);

    if(valid)
    {
        h->query    = result[1];
        h->q_start  = atoi( string(result[2]).c_str() );
        h->q_end    = atoi( string(result[3]).c_str() );
        h->q_strand = string(result[4]).at(0);

        h->node     = result[5];
        h->t_start  = atoi( string(result[6]).c_str() );
        h->t_end    = atoi( string(result[7]).c_str() );
        h->t_strand = string(result[8]).at(0);

        h->score    = atoi( string(result[9]).c_str() );

    }



    return valid;
}

void Exonerate_reads::write_exonerate_input(Node *root, vector<Fasta_entry> *reads, map<string,string> *names, int r)
{
    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,true);

    // create exonerate input
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    vector<Fasta_entry>::iterator it = reads->begin();
    for(;it!=reads->end();it++)
    {
        q_output<<">"<<it->name<<endl<<it->sequence<<endl;
    }
    q_output.close();


    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    it = aligned_sequences.begin();
    for(;it!=aligned_sequences.end();it++)
    {
        if(names->find(it->name) != names->end())
        {
            string seq = it->sequence;
            for (string::iterator si = seq.begin();si != seq.end();)
                if(*si == '-')
                    seq.erase(si);
                else
                    si++;
            t_output<<">"<<it->name<<endl<<seq<<endl;
        }
    }
    t_output.close();

    return;
}

void Exonerate_reads::write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int r)
{
    vector<Fasta_entry> aligned_sequences;
    root->get_alignment(&aligned_sequences,true);

    // create exonerate input
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    ofstream q_output( q_name.str().c_str(), (ios::out));
    q_output<<">"<<read->name<<endl<<read->sequence<<endl;
    q_output.close();

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    ofstream t_output( t_name.str().c_str(), (ios::out));
    vector<Fasta_entry>::iterator it = aligned_sequences.begin();
    for(;it!=aligned_sequences.end();it++)
    {
        if(names->find(it->name) != names->end())
        {
            string seq = it->sequence;
            for (string::iterator si = seq.begin();si != seq.end();)
                if(*si == '-')
                    seq.erase(si);
                else
                    si++;
            t_output<<">"<<it->name<<endl<<seq<<endl;
        }
    }
    t_output.close();

    return;
}

void Exonerate_reads::all_local_alignments(Node *root, vector<Fasta_entry> *reads, std::multimap<std::string,std::string> *tid_nodes, std::map<std::string,std::multimap<std::string,hit> > *hits, bool is_local)
{

    int r = rand();

    bool ignore_tid_tags = Settings_handle::st.is("test-every-internal-node") || Settings_handle::st.is("test-every-node") ;

    set<string> tid_tags;
    vector<Fasta_entry>::iterator ri = reads->begin();
    for(;ri!=reads->end();ri++)
    {
        tid_tags.insert(ri->tid);
    }

    map<string,string> names;
    multimap<string,string>::iterator it = tid_nodes->begin();
    for(;it!=tid_nodes->end();it++)
    {
        if( ignore_tid_tags )
        {
            names.insert(pair<string,string>(it->second,"empty"));
        }
        else
        {
            if(tid_tags.find(it->first)!=tid_tags.end())
            {
                names.insert(pair<string,string>(it->second,it->first));
            }
        }
    }

    this->write_exonerate_input(root,reads,&names,r);

    // exonerate command for local alignment

    stringstream command;
    if(is_local)
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";
    else
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -m affine:local -E 2>&1";

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }

    // read exonerate output, summing the multiple hit scores

    char line[256];
    map<string,multimap<string,hit> > all_hits;

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {

            map<string,multimap<string,hit> >::iterator iter = all_hits.find(h.query);

            if( iter != all_hits.end() )
            {

                multimap<string,hit>::iterator iter2 = iter->second.find(h.node);

                if( iter2 != iter->second.end() )
                {
                    if(iter2->second.t_strand == h.t_strand && iter2->second.q_strand == h.q_strand)
                    {
                        iter2->second.score += h.score;

                        if(iter2->second.q_start > h.q_start)
                            iter2->second.q_start = h.q_start;
                        if(iter2->second.q_end < h.q_end)
                            iter2->second.q_end = h.q_end;
                        if(iter2->second.t_start > h.t_start)
                            iter2->second.t_start = h.t_start;
                        if(iter2->second.t_end < h.t_end)
                            iter2->second.t_end = h.t_end;
                    }
                    else if(iter2->second.score < h.score)
                    {
                        iter2->second = h;
                    }
                }
                else
                {
                    iter->second.insert( make_pair(h.node, h) );
                }
            }
            else
            {
                multimap<string,hit> new_hit;
                new_hit.insert( make_pair(h.node, h) );
                all_hits.insert( make_pair(h.query, new_hit ) );
            }
        }
    }
    pclose(fpipe);


    tid_nodes->clear();
    hits->clear();

    ri = reads->begin();
    for( ;ri != reads->end(); ri++ )
    {

        map<string,multimap<string,hit> >::iterator iter = all_hits.find(ri->name);
        if( iter != all_hits.end() )
        {

            vector<hit> best_hits;

            multimap<string,hit>::iterator iter2 = iter->second.begin();
            for( ;iter2 != iter->second.end(); iter2++ )
            {
                bool tid_matches = false;

                if(names.find(iter2->first) != names.end() && names.find(iter2->first)->second == ri->tid)
                    tid_matches = true;

                if(ignore_tid_tags || tid_matches)
                {
                    best_hits.push_back(iter2->second);
                }
            }

            sort (best_hits.begin(), best_hits.end(), Exonerate_reads::better);

            multimap<string,hit> best_hits_for_this;

            if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                    Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

                ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                    Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
            {
                int lim = best_hits.at(0).score;
                if(is_local)
                    lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
                else
                    lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );

                for(int i=0; i<(int)best_hits.size(); i++)
                {
                    if(best_hits.at(i).score > lim)
                    {
                        string tid = names.find(best_hits.at(i).node)->second;
                        tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                        best_hits_for_this.insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                        Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
                    }
                }
            }

            // keep a fixed number of hits

            else
            {
                int lim = best_hits.size();
                if( is_local && Settings_handle::st.is("exonerate-local-keep-best") &&
                        Settings_handle::st.get("exonerate-local-keep-best").as<int>()>0 )
                    lim = Settings_handle::st.get("exonerate-local-keep-best").as<int>();

                if( !is_local && Settings_handle::st.is("exonerate-gapped-keep-best") &&
                        Settings_handle::st.get("exonerate-gapped-keep-best").as<int>()>0 )
                    lim = Settings_handle::st.get("exonerate-gapped-keep-best").as<int>();

                for(int i=0; i<lim && i<(int)best_hits.size(); i++)
                {
                    string tid = names.find(best_hits.at(i).node)->second;
                    tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                    best_hits_for_this.insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                    Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
                }
            }

            hits->insert( make_pair<string,multimap<string,hit> >(ri->name,best_hits_for_this) );
        }

        else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
        {
            ri->node_to_align = "discarded_read";
        }

    }

    ri = reads->begin();
    for( ;ri != reads->end(); ri++ )
    {

        map<string,multimap<string,hit> >::iterator iter = hits->find(ri->name);
        if( iter != hits->end() )
        {
            Log_output::write_out("Exonerate_reads: "+iter->first+" has "+Log_output::itos(iter->second.size())+" hits\n",0);

            multimap<string,hit>::iterator iter2 = iter->second.begin();

            for( ;iter2 != iter->second.end(); iter2++ )
            {
                Log_output::write_out("  "+ri->name+" matches "+iter2->first+" with score "+Log_output::itos(iter2->second.score)+"\n",0);

            }
        }
    }


    if(!Settings_handle::st.is("keep-exonerate-files"))
        this->delete_files(r);


}

void Exonerate_reads::local_alignment(Node *root, Fasta_entry *read, multimap<string,string> *tid_nodes, map<string,hit> *hits, bool is_local, bool all_nodes)
{

    int r = rand();

    map<string,string> names;
    multimap<string,string>::iterator it = tid_nodes->begin();
    for(;it!=tid_nodes->end();it++)
    {
        if(all_nodes)
        {
            names.insert(pair<string,string>(it->second,"empty"));
        }
        else
        {
            if(it->first==read->tid)
            {
                names.insert(pair<string,string>(it->second,it->first));
            }
        }
    }

    this->write_exonerate_input(root,read,&names,r);

    // exonerate command for local alignment

    stringstream command;
    if(is_local)
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no 2>&1";
    else
        command << "exonerate -q q"<<r<<".fas -t t"<<r<<".fas --showalignment no --showsugar yes --showvulgar no -m affine:local -E 2>&1";

    FILE *fpipe;
    if ( !(fpipe = (FILE*)popen(command.str().c_str(),"r")) )
    {
        Log_output::write_out("Problems with exonerate pipe.\nExiting.\n",0);
        exit(1);
    }

    // read exonerate output, summing the multiple hit scores

    char line[256];
    map<string,hit> all_hits;
    vector<string> hit_names;

    while ( fgets( line, sizeof line, fpipe))
    {
        hit h;
        bool valid = split_sugar_string(string(line),&h);

        if(valid)
        {
            map<string,hit>::iterator iter = all_hits.find(h.node);
            if( iter != all_hits.end() )
            {
                if(iter->second.t_strand == h.t_strand && iter->second.q_strand == h.q_strand)
                {
                    iter->second.score += h.score;

                    if(iter->second.q_start > h.q_start)
                        iter->second.q_start = h.q_start;
                    if(iter->second.q_end < h.q_end)
                        iter->second.q_end = h.q_end;
                    if(iter->second.t_start > h.t_start)
                        iter->second.t_start = h.t_start;
                    if(iter->second.t_end < h.t_end)
                        iter->second.t_end = h.t_end;
                }
                else if(iter->second.score < h.score)
                {
                    iter->second = h;
                }
            }
            else
            {
                all_hits.insert( make_pair(h.node, h) );
                hit_names.push_back(h.node);
            }
        }
    }
    pclose(fpipe);


    Log_output::write_out("Exonerate_reads: "+read->name+" has "+Log_output::itos(hit_names.size())+" hits\n",2);

    if(hit_names.size()>0)
    {
        tid_nodes->clear();
        hits->clear();

        vector<hit> best_hits;
        vector<string>::iterator iter = hit_names.begin();

        for(;iter!=hit_names.end();iter++)
        {
             map<string,hit>::iterator iter2 = all_hits.find(*iter);
             if( iter2 != all_hits.end() )
             {
               best_hits.push_back(iter2->second);
             }
        }

        sort (best_hits.begin(), best_hits.end(), Exonerate_reads::better);


        // keep hits that are above a relative threshold

        if( ( is_local && Settings_handle::st.is("exonerate-local-keep-above") &&
                Settings_handle::st.get("exonerate-local-keep-above").as<float>()>0 ) ||

            ( !is_local && Settings_handle::st.is("exonerate-gapped-keep-above") &&
                Settings_handle::st.get("exonerate-gapped-keep-above").as<float>()>0 ) )
        {
            int lim = best_hits.at(0).score;
            if(is_local)
                lim = int (lim * Settings_handle::st.get("exonerate-local-keep-above").as<float>() );
            else
                lim = int (lim * Settings_handle::st.get("exonerate-gapped-keep-above").as<float>() );


            for(int i=0; i<(int)hit_names.size(); i++)
            {
                if(best_hits.at(i).score > lim)
                {
                    string tid = names.find(best_hits.at(i).node)->second;
                    tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                    hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                    Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
                }
            }
        }


        // keep a fixed number of hits

        else
        {
            int lim = hit_names.size();
            if( is_local && Settings_handle::st.is("exonerate-local-keep-best") &&
                    Settings_handle::st.get("exonerate-local-keep-best").as<int>()>0 )
                lim = Settings_handle::st.get("exonerate-local-keep-best").as<int>();

            if( !is_local && Settings_handle::st.is("exonerate-gapped-keep-best") &&
                    Settings_handle::st.get("exonerate-gapped-keep-best").as<int>()>0 )
                lim = Settings_handle::st.get("exonerate-gapped-keep-best").as<int>();

            for(int i=0; i<lim && i<(int)hit_names.size(); i++)
            {
                string tid = names.find(best_hits.at(i).node)->second;
                tid_nodes->insert(pair<string,string>(tid,best_hits.at(i).node));
                hits->insert(pair<string,hit>(best_hits.at(i).node,best_hits.at(i)));

                Log_output::write_out("Exonerate_reads: adding "+best_hits.at(i).node+" "+Log_output::itos(best_hits.at(i).score)+"\n",3);
            }
        }
    }
    else if(!Settings_handle::st.is("keep-despite-exonerate-fails"))
    {
        tid_nodes->clear();
        read->node_to_align = "discarded_read";
    }

    if(!Settings_handle::st.is("keep-exonerate-files"))
        this->delete_files(r);
}

void Exonerate_reads::delete_files(int r)
{
    stringstream q_name;
    q_name <<"q"<<r<<".fas";

    stringstream t_name;
    t_name <<"t"<<r<<".fas";

    if( remove( q_name.str().c_str() ) != 0 )
       Log_output::write_out( "Error deleting file", 1);
    if( remove( t_name.str().c_str() ) != 0 )
       Log_output::write_out( "Error deleting file", 1 );
}
