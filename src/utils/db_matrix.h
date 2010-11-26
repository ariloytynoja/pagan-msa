
#ifndef DB_MATRIX_H
#define DB_MATRIX_H

#define NDEBUG

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#include <string>
#include <iostream>
#include <cassert>
#include <cmath>

namespace ppa {

class Db_matrix{
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
    double* data;
    int i,j,k,l;
public:
    Db_matrix(int x, std::string name="");
    Db_matrix(int x, int y, std::string name="");
    Db_matrix(int x, int y, int z, std::string name="");
    Db_matrix(int x, int y, int z, int w, std::string name="");

    ~Db_matrix();

    void allocate();
    void initialise(double v = 0);

    void allow_resize(bool xr, bool yr=false, bool zr=false, bool wr=false);

    void resize(int i);
    void copy_data(double *tmp,int new_x,int new_y,int new_z,int new_w);


    /*set*/
    void s(double v, int xa, int ya=0, int za = 0, int wa = 0);

    /*get*/
    double g(int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0); assert(xa<x);
        assert(ya>=0); assert(ya<y);
        assert(za>=0); assert(za<z);
        assert(wa>=0); assert(wa<w);
        return data[xa + ya*x + za*x*y + wa*x*y*z];
    }

    /*add*/
    void a(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0); assert(xa<x);
        assert(ya>=0); assert(ya<y);
        assert(za>=0); assert(za<z);
        assert(wa>=0); assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] += v;
    }

    /*divide*/
    void d(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0); assert(xa<x);
        assert(ya>=0); assert(ya<y);
        assert(za>=0); assert(za<z);
        assert(wa>=0); assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] /= v;
    }

    /*multiply*/
    void m(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0); assert(xa<x);
        assert(ya>=0); assert(ya<y);
        assert(za>=0); assert(za<z);
        assert(wa>=0); assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] *= v;
    }

//    void alog(double v, int xa, int ya=0, int za = 0, int wa = 0) { assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z);  assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] = sumLogs(data[xa + ya*x + za*x*y + wa*x*y*z], v); }
//    void clog(int xa, int ya=0, int za = 0, int wa = 0) { assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z); assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] = log(data[xa + ya*x + za*x*y + wa*x*y*z]); }
//    double sumLogs(double a, double b);

    void printName()
    {
        std::cout<<"Name "<<name<<": x = "<<x<<", y = "<<y<<", z = "<<z<<", w = "<<w<<std::endl;
    }

    void print();
    void print(int m);

};

}
#endif
