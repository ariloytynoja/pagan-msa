
#ifndef INT_MATRIX_H
#define INT_MATRIX_H

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#include <string>
#include <iostream>
#include <cassert>

class Int_matrix{
private:
    int x;
    int y;
    int z;
    int w;
	bool xar;
	bool yar;
	bool zar;
	bool war;

	std::string name;
    int* data;
    int i,j,k,l;
public:
    Int_matrix(int x, std::string name="");
    Int_matrix(int x, int y, std::string name="");
    Int_matrix(int x, int y, int z, std::string name="");
    Int_matrix(int x, int y, int z, int w, std::string name="");

    ~Int_matrix();

    void allocate();
    void initialise(int v = 0);

    int g(int xa, int ya=0, int za = 0, int wa = 0) {  assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z); assert(wa>=0); assert(wa<w); return data[xa + ya*x + za*x*y + wa*x*y*z]; }
/*    void s(int v, int xa, int ya=0, int za = 0, int wa = 0) { assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z);  assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] = v; }*/
	void s(int v, int xa, int ya=0, int za = 0, int wa = 0);
	void a(int v, int xa, int ya=0, int za = 0, int wa = 0) { assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z);  assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] += v; }
    void printName() { std::cout<<"Name "<<name<<": x = "<<x<<", y = "<<y<<", z = "<<z<<", w = "<<w<<std::endl; }
    void print();
	void print(int m);

	int X() { return x; }

	void resize(int i);
	void copyData(int *tmp,int new_x,int new_y,int new_z,int new_w);

	void allowResize(bool xr, bool yr=false, bool zr=false, bool wr=false);
};

#endif
