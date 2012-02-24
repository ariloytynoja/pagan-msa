#ifndef EXONERATE_READS_H
#define EXONERATE_READS_H

#include "utils/fasta_entry.h"
#include "main/node.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>

namespace ppa{

struct hit {
    std::string query;
    std::string node;
    int score;
    int q_start;
    int q_end;
    char q_strand;
    int t_start;
    int t_end;
    char t_strand;
};

class Exonerate_reads
{
    static bool better (hit i,hit j) { return (i.score>j.score); }

    bool split_sugar_string(const std::string& row,hit *h);
    bool split_vulgar_string(const std::string& row,hit *h);
    void write_exonerate_input(Node *root, vector<Fasta_entry> *reads, map<string,string> *names, int r);
    void write_exonerate_input(Node *root, Fasta_entry *read, map<string,string> *names, int r);
    void delete_files(int r);

    string get_temp_dir()
    {
        string tmp_dir = "/tmp/";
        if(Settings_handle::st.is("temp-folder"))
            tmp_dir = Settings_handle::st.get("temp-folder").as<string>()+"/";

        struct stat st;
        if(stat(tmp_dir.c_str(),&st) != 0)
            tmp_dir = "";

        return tmp_dir;
    }

public:
    Exonerate_reads();
    bool test_executable();

    void local_alignment(Node *root, Fasta_entry *read, std::multimap<std::string,std::string> *good_hits, std::map<std::string,hit> *hits, bool is_local, bool all_nodes=false);
    void all_local_alignments(Node *root, vector<Fasta_entry> *reads, std::multimap<std::string,std::string> *tid_nodes, std::map<std::string,std::multimap<std::string,hit> > *hits, bool is_local);

};

}

#endif // EXONERATE_READS_H
