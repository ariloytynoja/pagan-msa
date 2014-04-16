/***************************************************************************
 *   Copyright (C) 2010-2014 by Ari Loytynoja                              *
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
#include "utils/log_output.h"
#include <iostream>
#include <algorithm>

using namespace ppa;
using namespace std;

Find_anchors::Find_anchors()
{
}


void Find_anchors::find_long_substrings(std::string *seq1,std::string *seq2,std::vector<Substring_hit> *hits,int min_length)
{
    len1 = seq1->length();
    len2 = seq2->length();

//    cout<<*seq1<<endl<<*seq2<<endl<<endl;
//    cout<<">root\n"<<*seq1<<endl;

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
                s.score = length;
                hits->push_back(s);
            }
        }
    }

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_length);

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
        bool overlap = false;
        for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
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
            for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
            {
                hit_site1.at(i)=true;
                hit_site2.at(j)=true;
            }
            it1++;
        }
    }
}


void Find_anchors::check_hits_order_conflict(std::string *seq1,std::string *seq2,vector<Substring_hit> *hits)
{
    len1 = seq1->length();
    len2 = seq2->length();

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_score);

    vector<bool> hit_site1;
    hit_site1.reserve(len1);
    for(int i=0;i<len1;i++)
        hit_site1.push_back(false);

    vector<bool> hit_site2;
    hit_site2.reserve(len2);
    for(int i=0;i<len2;i++)
        hit_site2.push_back(false);

//    cout<<"\nHits "<<hits->size()<<endl;

    vector<Substring_hit>::iterator it1 = hits->begin();
    for(;it1!=hits->end();)
    {
        bool overlap = false;
        for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
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
            for(int i=it1->start_site_1,j=it1->start_site_2;i<it1->start_site_1+it1->length && j<it1->start_site_2+it1->length;i++,j++)
            {
                hit_site1.at(i)=true;
                hit_site2.at(j)=true;
            }
            it1++;
        }
    }

    sort(hits->begin(),hits->end(),Find_anchors::sort_by_start_site_1);

    it1 = hits->begin();
    vector<Substring_hit>::iterator it2 = it1;
    it2++;
    for(;it1!=hits->end() && it2!=hits->end();)
    {
        if(it1->start_site_2 > it2->start_site_2)
        {
            if(it1->score < it2->score)
            {
                hits->erase(it1);
                it2 = it1;
                it2++;
            }
            else
            {
                hits->erase(it2);
                it2 = it1;
                it2++;
            }
            continue;
        }
        it1++;it2++;
    }


//    cout<<"Hits "<<hits->size()<<endl;

//    if(hits->size()>0)
//    {
//        vector<Substring_hit>::iterator it1 = hits->begin();
//        for(;it1!=hits->end();it1++)
//        {
//            cout<<"REMAINS "<<it1->start_site_1<<".."<<it1->start_site_1+it1->length<<" "<<it1->start_site_2<<".."<<it1->start_site_2+it1->length<<" ("<<it1->score<<")\n";
//        }
//    }
}


void Find_anchors::define_tunnel(std::vector<Substring_hit> *hits,std::vector<int> *upper_bound,std::vector<int> *lower_bound,string str1,string str2)
{
    int length1 = str1.length();
    int length2 = str2.length();

//    cout<<"\n>S1\n"<<str1<<"\n>S2\n"<<str2<<"\n";

    // Has to consider skpped-over gaps as those are in the sequence objects
    //
    vector<int> index1;
    vector<int> index2;

    for(int i=0;i<length1;i++)
        if(str1.at(i)!='-')
            index1.push_back(i+1);

    for(int i=0;i<length2;i++)
        if(str2.at(i)!='-')
            index2.push_back(i+1);


    upper_bound->reserve(length1+1);
    lower_bound->reserve(length1+1);

    vector<int> diagonals;
    diagonals.reserve(length1+1);
    for(int i=0;i<length1+1;i++)
        diagonals.push_back(-1);

    for(vector<Substring_hit>::iterator it = hits->begin();it!=hits->end();it++)
    {
//        cout<<it->start_site_1<<" "<<it->start_site_2<<" "<<it->length<<" "<<length1<<" "<<length2<<"\n";
        int i=0;
        for(;i<it->length;i++)
        {
            diagonals.at(index1.at(it->start_site_1+i)) = index2.at(it->start_site_2+i);
        }

        if(it->start_site_1+i < (int)index1.size() && index1.at(it->start_site_1+i) < (int)diagonals.size())
            diagonals.at(index1.at(it->start_site_1+i)) = -2;

    }


    int width = Settings_handle::st.get("anchors-offset").as<int>();

    int y1 = 0;
    int y2 = 0;
    int y;

    int prev_y = 0;
    int m_count = 0;

    for(int i=0;i<=length1;i++)
    {
        if(i>=width && diagonals.at(i-width)>=0)
            y1 = diagonals.at(i-width)+0;

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)-width+0;
        }

        if(diagonals.at(i)>=0 && i>0 && diagonals.at(i-1)+1 == diagonals.at(i))
        {
            m_count++;
        }
        else if(diagonals.at(i)==-2)
        {
            m_count = 0;
        }

        y = min(y1,y2);
        y = max(y,0);

        if(diagonals.at(i)>=0 && i>0 && diagonals.at(i-1)+1 == diagonals.at(i) && m_count>=width)
        {
            prev_y = y;
        }

//        if(i>0)
//            cout<<i<<" "<<y1<<" "<<y2<<" "<<prev_y<<"; "<<diagonals.at(i)<<" "<<diagonals.at(i-1)<<" "<<m_count<<endl;

        y = min(y,prev_y);
        y = max(y,0);

        upper_bound->push_back(y);

    }

    y1 = length2;
    y2 = length2;

    prev_y = length2;
    m_count = 0;

    for(int i=length1;i>=0;i--)
    {
        if(i<=length1-width && diagonals.at(i+width)>=0)
            y1 = diagonals.at(i+width)+0;

        if(diagonals.at(i)>=0)
        {
            y2 = diagonals.at(i)+width+0;
        }

        if(diagonals.at(i)>=0 && i<length1 && diagonals.at(i+1)-1 == diagonals.at(i))
        {
            m_count++;
        }
        else if(diagonals.at(i)==-2)
        {
            m_count = 0;
        }

        y = max(y1,y2);
        y = min(y,length2);

        if(diagonals.at(i)>=0 && i<length1 && diagonals.at(i+1)-1 == diagonals.at(i)&& m_count>=width)
        {
            prev_y = y;
        }

        y = max(y,prev_y);
        y = min(y,length2);

        lower_bound->insert(lower_bound->begin(),y);
    }

    /*
    if(Settings::noise>0)
    {
        int sum = 0;
        for(int i=0;i<length1;i++)
            sum += lower_bound->at(i)-upper_bound->at(i);

        stringstream s;
        s.precision(4);
        s<<"Anchoring: Computing "<<((float)sum/(length1*length2))*100.0<<"% of DP matrix.\n";
        Log_output::write_out(s.str(),1);
    }
    */

    if(Settings_handle::st.is("plot-anchors-for-R")){
        cout<<"\n\nl1="<<str1.length()<<"; l2="<<str2.length()<<"; plot(1,1,type=\"n\",xlim=c(0,l1),ylim=c(0,l2))\n";
        cout<<"segments(0,0,l1,0,col=\"red\"); segments(0,l2,l1,l2,col=\"red\"); segments(0,0,0,l2,col=\"red\"); segments(l1,0,l1,l2,col=\"red\")\n";
        stringstream xh;
        stringstream yh;
        int lastx = 0,lasty = 0;
        for(int i=0;i<(int)diagonals.size();i++)
        {
            if(diagonals.at(i)>=0)
                {
                xh<<i<<",";
                yh<<diagonals.at(i)<<",";
                lastx=i;lasty=diagonals.at(i);
            }
        }
        xh<<lastx;yh<<lasty;
        cout<<"xhit=c("<<xh.str()<<");yhit=c("<<yh.str()<<");points(xhit,yhit,pch=\".\")\n";
        cout<<"upper=c("<<upper_bound->at(0);
        for(int i=1;i<(int)upper_bound->size();i++)
            cout<<","<<upper_bound->at(i);
        cout<<")\nlower=c("<<lower_bound->at(0);
        for(int i=1;i<(int)lower_bound->size();i++)
            cout<<","<<lower_bound->at(i);
        cout<<")\nlines(0:l1,upper,col=\"green\"); lines(0:l1,lower,col=\"blue\")\n\n";
    }

}
