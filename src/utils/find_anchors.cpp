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

#include "utils/find_anchors.h"
#include "utils/settings_handle.h"
#include <iostream>
#include <algorithm>

using namespace ppa;
using namespace std;

Find_anchors::Find_anchors()
{
}


void Find_anchors::find_long_substrings(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits,int min_length)
{
//    cout<<"find_long_substrings\n";

    len1 = seq1->length();
    len2 = seq2->length();

    char c1[len1];
    char c2[len2];
    char *a[len1+len2];

    int m=0; int n=0;
    for(; n<len1; n++,m++) {
        c1[n] = seq1->at(n);
        a[m] = &c1[n];
    }
    c1[n] = 0;

    int l = m;

    n=0;
    for(; n<len2; n++,m++) {
        c2[n] = seq2->at(n);
        a[m] = &c2[n];
    }
    c2[n] = 0;

    beg1 = a[0];
    beg2 = a[l];

    qsort(a, m, sizeof(char *), Find_anchors::pstrcmp);

    for (int i = 0; i < m-1; i++)
    {
        if(different_strings(a[i], a[i+1]))
        {
            int length = this->identical_prefix_length(a[i], a[i+1]);
            if (length >= min_length) {
                int pos1 = position_string1(a[i],a[i+1]);
                int pos2 = position_string2(a[i],a[i+1]);

                Substring_hit s;
                s.start_site_1 = pos1;
                s.start_site_2 = pos2;
                s.length = length;
                hits->push_back(s);
            }
        }
    }

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_length);
//    cout<<"sorted "<<hits->size()<<"\n";

    vector<bool> hit_site1;
    hit_site1.reserve(len1);
    for(int i=0;i<len1;i++)
        hit_site1.push_back(false);

    vector<bool> hit_site2;
    hit_site2.reserve(len2);
    for(int i=0;i<len2;i++)
        hit_site2.push_back(false);

    vector<Substring_hit>::iterator it1 = hits->begin();
    for(;it1!=hits->end();)
    {
//        cout<<"hit "<<it1->start_site_1+it1->length<<endl;

        bool overlap = false;
        for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length,j<it1->start_site_2+it1->length;i++,j++)
        {
            if(hit_site1.at(i) || hit_site2.at(j))
            {
                overlap = true;
                break;
            }
        }


        if(overlap)
        {
            hits->erase(it1);
        }
        else
        {
            for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length,j<it1->start_site_2+it1->length;i++,j++)
            {
                hit_site1.at(i)=true;
                hit_site2.at(j)=true;
            }
            it1++;
        }
    }
}

void Find_anchors::check_hits_order_conflict(vector<Substring_hit> *hits)
{
    sort(hits->begin(),hits->end(),Find_anchors::sort_by_start_site_1);

    bool order_conflicts = false;

    if(hits->size()>1)
    {
        vector<Substring_hit>::iterator it1 = hits->begin();
        vector<Substring_hit>::iterator it2 = hits->begin();
        it2++;
        for(;it2!=hits->end();it1++,it2++)
        {
            if(it1->start_site_1>it2->start_site_1 || it1->start_site_2>it2->start_site_2)
                order_conflicts = true;
        }
    }

    if(order_conflicts)
    {
        hits->clear();
    }
}


void Find_anchors::define_tunnel(std::vector<Substring_hit> *hits,std::vector<int> *upper_bound,std::vector<int> *lower_bound,string str1,string str2)
{
//    cout<<"define_tunnel\n";

    int length1 = str1.length();
    int length2 = str2.length();

    // Has to consider skpped-over gaps as those are in the sequence objects
    //
    vector<int> index1;
    vector<int> index2;
    index1.push_back(0);
    index2.push_back(0);

    for(int i=0;i<length1;i++)
        if(str1.at(i)!='-')
            index1.push_back(i);

    for(int i=0;i<length2;i++)
        if(str2.at(i)!='-')
            index2.push_back(i);

    upper_bound->reserve(length1);
    lower_bound->reserve(length1);

    vector<int> diagonals;
    diagonals.reserve(length1);
    for(int i=0;i<length1;i++)
        diagonals.push_back(-1);

    for(vector<Substring_hit>::iterator it = hits->begin();it!=hits->end();it++)
    {
//        cout<<it->start_site_1<<" "<<it->start_site_2<<" "<<it->length<<" "<<length1<<" "<<length2<<"\n";
        for(int i=0;i<it->length;i++)
            diagonals.at(index1.at(it->start_site_1+i)) = index2.at(it->start_site_2+i);
    }

    int width = Settings_handle::st.get("anchors-tunnel-offset").as<int>();

    int y1 = 0;
    int y2 = 0;
    int y;

    int prev_y = 0;
    int m_count = 0;

    upper_bound->push_back(0);

    for(int i=0;i<length1;i++)
    {
        if(i>=width && diagonals.at(i-width)>=0)
            y1 = diagonals.at(i-width);

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)-width;
            m_count++;
        }
        else
        {
            m_count = 0;
        }

        y = min(y1,y2);
        y = max(y,0);

        if(diagonals.at(i)>=0 && m_count>=width)
        {
            prev_y = y;
        }

        y = min(y,prev_y);

        upper_bound->push_back(y);
    }

    y1 = length2;
    y2 = length2;

    prev_y = length2;
    m_count = 0;

    for(int i=length1-1;i>=0;i--)
    {
        if(i<=length1-1-width && diagonals.at(i+width)>=0)
            y1 = diagonals.at(i+width);

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)+width;
            m_count++;
        }
        else
        {
            m_count = 0;
        }

        y = max(y1,y2);
        y = min(y,length2);

        if(diagonals.at(i)>=0 && m_count>=width)
        {
            prev_y = y;
        }

        y = max(y,prev_y);

        lower_bound->push_back(y);
    }
    lower_bound->push_back(y);
    reverse(lower_bound->begin(),lower_bound->end());
}
