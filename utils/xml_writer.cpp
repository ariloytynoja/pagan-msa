#include <iostream>

#include "xml_writer.h"

using namespace ppa;

Xml_writer::Xml_writer()
{

}

/****************************************************************************************/

void Xml_writer::write(ostream & output, const Node *root, const vector<Fasta_entry> & seqs) const throw (Exception)
{
    // Checking the existence of specified file, and possibility to open it in write mode
    if (! output) { throw IOException ("Xml_writer::write. Failed to open file"); }


    output << "<ms_alignment>\n<newick>" << root->print_xml_tree() << "</newick>\n<nodes>\n";

    string seq, temp = "";  // Initialization

    vector<Fasta_entry>::const_iterator vi = seqs.begin();

    // Main loop : for all sequences in vector container
    for (; vi != seqs.end(); vi++)
    {
        string id = root->get_id_for_name(vi->name);

        output << "<leaf id=\"" << id <<"\" name=\"" << vi->name << "\">\n";
        output << "  <sequence>\n    " << vi->sequence << "\n  </sequence>\n</leaf>\n";
    }

    output << "</nodes>\n</ms_alignment>\n";
}

/****************************************************************************************/
