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

#include <cstdlib>
#include <iostream>
#include "utils/int_matrix.h"
#include "utils/settings.h"


using namespace std;
using namespace ppa;

Int_matrix::Int_matrix(int xa, std::string n)
{
    assert(xa>0);
    x = xa;
    y = z = w = 1;
    name = n;

    allocate();
	xar=yar=zar=war=false;
}

Int_matrix::Int_matrix(int xa, int ya, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    z = w = 1;
    name = n;

    allocate();
	xar=yar=zar=war=false;
}

Int_matrix::Int_matrix(int xa, int ya, int za, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    assert(za>0);
    z = za;
    w = 1;
    name = n;

    allocate();
	xar=yar=zar=war=false;
}

Int_matrix::Int_matrix(int xa, int ya, int za, int wa, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    assert(za>0);
    z = za;
    assert(wa>0);
    w = wa;
    name = n;

    allocate();
	xar=yar=zar=war=false;
}

Int_matrix::~Int_matrix()
{
//    cout<<"int delete "<<name<<endl;
    delete []data;
    x=y=z=w=0;
}

void Int_matrix::allocate()
{
    data = new int[x*y*z*w];
}

void Int_matrix::initialise(int v)
{
    FOR(i,x) {
        FOR(j,y) {
            FOR(k,z) {
                FOR(l,w) {
                    data[i + j*x + k*x*y + l*x*y*z] = v;
                }
            }
        }
    }
}

void Int_matrix::s(int v, int xa, int ya, int za, int wa)
{
	assert(xa>=0);
	assert(ya>=0);
	assert(za>=0);
	assert(wa>=0);

	if(xa>=x && xar){
		resize(1);
		this->s(v,xa,ya,za,wa);
	}
	else if(xa>=x){
        cout<<"Int_matrix("<<name<<"): x ("<<xa<<") over the limit ("<<x<<")!"<<endl;
		exit(-1);
	}
	if(ya>=y && yar){
		resize(2);
		this->s(v,xa,ya,za,wa);
	}
	else if(ya>=y){
        cout<<"Int_matrix("<<name<<"): y ("<<ya<<") over the limit ("<<y<<")!"<<endl;
		exit(-1);
	}
	if(za>=z && zar){
		resize(3);
		this->s(v,xa,ya,za,wa);
	}
	else if(za>=z){
        cout<<"Int_matrix("<<name<<"): z ("<<za<<") over the limit ("<<z<<")!"<<endl;
		exit(-1);
	}
	if(wa>=w && war){
		resize(4);
		this->s(v,xa,ya,za,wa);
	}
	else if(wa>=w){
        cout<<"Int_matrix("<<name<<"): w ("<<wa<<") over the limit ("<<w<<")!"<<endl;
		exit(-1);
	}

	data[xa + ya*x + za*x*y + wa*x*y*z] = v;
}

void Int_matrix::resize(int i)
{
    assert(Settings::resize_factor>1);

	if(i==1){
        int new_x = (int)(Settings::resize_factor*x);
		if(new_x == x)
			new_x++;
		int *tmp = new int[new_x*y*z*w];
		copyData(tmp,new_x,y,z,w);
		delete[] data;
		data = tmp;
		x = new_x;
	}
	else if(i==2){
        int new_y = (int)(Settings::resize_factor*y);
		if(new_y == y)
			new_y++;
		int *tmp = new int[x*new_y*z*w];
		copyData(tmp,x,new_y,z,w);
		delete[] data;
		data = tmp;
		y = new_y;
	}
	else if(i==3){
        int new_z = (int)(Settings::resize_factor*z);
		if(new_z == z)
			new_z++;
		int *tmp = new int[x*y*new_z*w];
		copyData(tmp,x,y,new_z,w);
		delete[] data;
		data = tmp;
		z = new_z;
	}
	else if(i==4){
        int new_w = (int)(Settings::resize_factor*w);
		if(new_w == w)
			new_w++;
		int *tmp = new int[x*y*z*new_w];
		copyData(tmp,x,y,z,new_w);
		delete[] data;
		data = tmp;
		w = new_w;
	}
}

void Int_matrix::copyData(int *tmp,int new_x,int new_y,int new_z,int )
{
	cout<<"Resizing matrix '"<<name<<"': consider using a greater initial size!"<<endl;
	FOR(i,x) {
		FOR(j,y) {
			FOR(k,z) {
				FOR(l,w) {
					tmp[i + j*new_x + k*new_x*new_y + l*new_x*new_y*new_z] = data[i + j*x + k*x*y + l*x*y*z];
				}
			}
		}
	}
}

void Int_matrix::print()
{
    FOR(i,x) {
        FOR(j,y) {
            FOR(k,z) {
                FOR(l,w) {
                    cout<<data[i + j*x + k*x*y + l*x*y*z]<<" ";
                }
                if (w>1)
                    cout<<endl;
            }
            if (z>1)
                cout<<endl;
        }
        if (y>1)
            cout<<endl;
    }
    if (x>1)
        cout<<endl;
}

void Int_matrix::print(int m)
{
	FOR(i,m) {
		FOR(j,y) {
			FOR(k,z) {
				FOR(l,w) {
					cout<<data[i + j*x + k*x*y + l*x*y*z]<<" ";
				}
				if (w>1)
					cout<<endl;
			}
			if (z>1)
				cout<<endl;
		}
		if (y>1)
			cout<<endl;
	}
	if (x>1)
		cout<<endl;
}

void Int_matrix::allowResize(bool xr, bool yr, bool zr, bool wr)
{
	xar=xr;
	yar=yr;
	zar=zr;
	war=wr;
}
