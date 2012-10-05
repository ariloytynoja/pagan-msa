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

#ifndef FIND_SUBSTRINGS_H
#define FIND_SUBSTRINGS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "substring_hit.h"

/*
Adapted from http://www.drdobbs.com/architecture-and-design/algorithm-alley/184404588
*/

namespace ppa{

class Find_substrings
{
    static int pstrcmp(const void *p, const void *q)
    {
        const char **cp = (const char **)p;
        const char **cq = (const char **)q;
        return strcmp(*cp, *cq);
    }

    static bool sort_by_length(Substring_hit p, Substring_hit q)
    {
        return (p.length > q.length);
    }

    static bool sort_by_start_site_1(Substring_hit p, Substring_hit q)
    {
        return (p.start_site_1 < q.start_site_1);
    }

    int identical_prefix_length(char *p, char *q)
    {
        int i = 0;
        while (*p && (*p++ == *q++))
            i++;
        return i;
    }

    bool different_strings(char *p, char *q)
    {
        if(p-beg1>=0 && p-beg1<len1 && q-beg2>=0 && q-beg2<len2)
            return true;
        if(q-beg1>=0 && q-beg1<len1 && p-beg2>=0 && p-beg2<len2)
            return true;

        return false;
    }

    int position_string1(char *p, char *q)
    {
        if(p-beg1>=0 && p-beg1<len1)
            return p-beg1;
        if(q-beg1>=0 && q-beg1<len1)
            return q-beg1;

        return -1;
    }

    int position_string2(char *p, char *q)
    {
        if(p-beg2>=0 && p-beg2<len2)
            return p-beg2;
        if(q-beg2>=0 && q-beg2<len2)
            return q-beg2;

        return -1;
    }

    bool hits_overlap(Substring_hit h1,Substring_hit h2)
    {
        if(h1.start_site_1 <= h2.start_site_1 && h1.start_site_1+h1.length >= h2.start_site_1 )
            return true;
        if(h1.start_site_2 <= h2.start_site_2 && h1.start_site_2+h1.length >= h2.start_site_2 )
            return true;

        return false;
    }


    char *beg1;
    char *beg2;
    int len1;
    int len2;


public:
    Find_substrings();
    void find_long_substrings(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits,int min_length);
    void define_tunnel(std::vector<Substring_hit> *hits,std::vector<int> *upper_bound,std::vector<int> *lower_bound,int length1, int length2);
};

}
#endif // FIND_SUBSTRINGS_H
