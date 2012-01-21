#ifndef LOG_OUTPUT_H
#define LOG_OUTPUT_H

#include <ostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

namespace ppa
{

class Log_output
{
    static ofstream fs;
    static ostream *os;
    static bool newline;
    static int header_length;
    static int msg_length;
    static int msg2_length;
    static string prev_msg;
    static string prev_msg2;
public:
    Log_output();
    static void open_stream();
    static void write_out(const string str,const int priority);
    static void write_msg(const string str,const int priority);
    static void append_msg(const string str,const int priority);
    static void write_header(const string str,const int priority);
    static void flush() {os->flush();}
    static void write_out(const string str,const string option);
    static void clean_output();
    static string itos(int i)
    {
        stringstream s;
        s << i;
        return s.str();
    }

};

}

#endif // LOG_OUTPUT_H
