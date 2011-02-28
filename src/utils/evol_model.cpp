/***************************************************************************
 *   Copyright (C) 2010 by Ari Loytynoja                                   *
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

#include "utils/evol_model.h"

using namespace ppa;

Evol_model::Evol_model(std::string alpha,float dist)
{
    full_char_alphabet = alpha;
    int char_fas = full_char_alphabet.length();

    distance = dist;

    charPi = new Db_matrix(char_fas,"pi_char");
    charPr = new Db_matrix(char_fas,char_fas,"P_char");
    charPr->initialise(0);

    logCharPi = new Db_matrix(char_fas,"logpi_char");
    logCharPr = new Db_matrix(char_fas,char_fas,"logP_char");
    logCharPr->initialise(-HUGE_VAL);

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

}

Evol_model::~Evol_model()
{
    delete charPi;
    delete charPr;

    delete logCharPi;
    delete logCharPr;

    delete parsimony_table;
}
