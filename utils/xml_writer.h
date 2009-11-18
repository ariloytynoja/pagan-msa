#ifndef XML_WRITER_H
#define XML_WRITER_H

#include <string>
#include <vector>
#include <fstream>
#include "exceptions.h"
#include "node.h"
#include "settings.h"
#include "utils/fasta_entry.h"

using namespace std;

namespace ppa
{

class Xml_writer
{
public:
    Xml_writer();

    void write(ostream & output, const Node *root, const vector<Fasta_entry> & seqs) const throw (Exception);
    void write(const string & path, const Node *root, const vector<Fasta_entry> & seqs, bool overwrite=true) const throw (Exception)
    {
        ofstream output( (path+".xml").c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
        write(output, root, seqs);
        output.close();
    }

};

}
#endif // XML_WRITER_H
