/***************************************************************************
 *   Copyright (C) 2012 by Ari Loytynoja                                   *
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

#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <vector>
#include "utils/fasta_entry.h"
#include "main/node.h"
#include "utils/model_factory.h"

namespace ppa{

class Input_output_parser
{
public:
    Input_output_parser();

    void parse_input_sequences(std::vector<Fasta_entry> *sequences,bool *reference_alignment);
    Node * parse_input_tree(vector<Fasta_entry> *sequences,bool reference_alignment);
    void match_sequences_and_tree(std::vector<Fasta_entry> *sequences,Node *root, bool reference_alignment, int *data_type);
    void define_alignment_model(Model_factory *mf,int data_type);
    void output_aligned_sequences(std::vector<Fasta_entry> *sequences,Node *root);
};
}

#endif // INPUT_PARSER_H
