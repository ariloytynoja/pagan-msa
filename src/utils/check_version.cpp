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

#include "utils/check_version.h"
#include "utils/log_output.h"

#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <netdb.h>
#include <string>
#include <sstream>
#include <iostream>
#include <unistd.h>

using namespace std;
using namespace ppa;

#define PORT 80

Check_version::Check_version(float version)
{

    stringstream ss;
    ss<<"\nThis is PAGAN v."<<version<<".\nChecking if updates are available at http://http://code.google.com/p/pagan-msa.\n";
    Log_output::write_out(ss.str(),0);

    struct sockaddr_in *remote;
    char buf[BUFSIZ+1];

    int sock = create_tcp_socket();
    char *ip = get_ip("pagan-msa.googlecode.com");

    remote = (struct sockaddr_in *)malloc(sizeof(struct sockaddr_in *));
    remote->sin_family = AF_INET;
    int tmpres = inet_pton(AF_INET, ip, (void *)(&(remote->sin_addr.s_addr)));
    if( tmpres < 0)
    {
        Log_output::write_out("Can't set remote->sin_addr.s_addr",0);
        exit(1);
    }
    else if(tmpres == 0)
    {
        stringstream ss;
        ss<<ip<<" is not a valid IP address\n";
        Log_output::write_out(ss.str(),0);

        exit(1);
    }
    remote->sin_port = htons(PORT);

    if(connect(sock, (struct sockaddr *)remote, sizeof(struct sockaddr)) < 0){
        Log_output::write_out("Could not connect",0);
        exit(1);
    }

    char get[] = "GET /git/VERSION_HISTORY HTTP/1.0\r\nHost: pagan-msa.googlecode.com\r\nUser-Agent: HTMLGET 1.0\r\n\r\n";

    //Send the query to the server
    int sent = 0;
    while(sent < (int)strlen(get))
    {
        tmpres = send(sock, get+sent, strlen(get)-sent, 0);
        if(tmpres == -1)
        {
            Log_output::write_out("Can't send query",0);
            exit(1);
        }
        sent += tmpres;
    }


    stringstream output;

    //now it is time to receive the page
    memset(buf, 0, sizeof(buf));
    int htmlstart = 0;
    char * htmlcontent;
    while((tmpres = recv(sock, buf, BUFSIZ, 0)) > 0){
        if(htmlstart == 0)
        {
            /* Under certain conditions this will not work.
            * If the \r\n\r\n part is splitted into two messages
            * it will fail to detect the beginning of HTML content
            */
            htmlcontent = strstr(buf, "\r\n\r\n");
            if(htmlcontent != NULL)
            {
                htmlstart = 1;
                htmlcontent += 4;
            }
        }
        else
        {
            htmlcontent = buf;
        }
        if(htmlstart)
        {
            output<<htmlcontent;
        }

        memset(buf, 0, tmpres);
    }
    if(tmpres < 0)
    {
      Log_output::write_out("Error receiving data",0);
    }

    bool print_this = true;
    bool has_printed = false;
    string s;

    while( getline(output,s) )
    {
        istringstream ss(s);
        double d;
        char v,p;
        while( ss >> v >> p >> d )
        {
            if(v=='v' && p=='.' && int(d*10000) <= int(version*10000)+10)
            {
               print_this = false;
            }
        }

        if(print_this)
        {
            if(!has_printed)
                Log_output::write_out("\nFound updates. Changes in the more recent versions:\n\n",0);

            has_printed = true;
            Log_output::write_out(s+"\n",0);
        }
        else
        {
            break;
        }
    }

    if(!has_printed)
        Log_output::write_out("\nNo updates are available.\n\n",0);

    free(remote);
    free(ip);
    close(sock);
}

int Check_version::create_tcp_socket()
{
  int sock;
  if((sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0){
    Log_output::write_out("Can't create TCP socket",0);
    exit(1);
  }
  return sock;
}


char *Check_version::get_ip(const char *host)
{
  struct hostent *hent;
  int iplen = 15; //XXX.XXX.XXX.XXX
  char *ip = (char *)malloc(iplen+1);
  memset(ip, 0, iplen+1);
  if((hent = gethostbyname(host)) == NULL)
  {
    Log_output::write_out("Can't get IP",0);
    exit(1);
  }
  if(inet_ntop(AF_INET, (void *)hent->h_addr_list[0], ip, iplen) == NULL)
  {
    Log_output::write_out("Can't resolve host",0);
    exit(1);
  }
  return ip;
}

